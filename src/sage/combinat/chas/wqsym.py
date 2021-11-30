# -*- coding: utf-8 -*-
r"""
Word Quasi-symmetric functions

AUTHORS:

- Travis Scrimshaw (2018-04-09): initial implementation
- Darij Grinberg and Amy Pang (2018-04-12): further bases and methods

TESTS:

We check that the coercion `C \to M` goes through the `X` basis::

    sage: WQSym = algebras.WQSym(QQ)
    sage: Q = WQSym.Q()
    sage: C = WQSym.C()
    sage: M = WQSym.M()
    sage: M.coerce_map_from(C)
    Composite map:
      From: Word Quasi-symmetric functions over Rational Field in the Cone basis
      To:   Word Quasi-symmetric functions over Rational Field in the Monomial basis
      Defn:   Generic morphism:
              From: Word Quasi-symmetric functions over Rational Field in the Cone basis
              To:   Word Quasi-symmetric functions over Rational Field in the Characteristic basis
            then
              Generic morphism:
              From: Word Quasi-symmetric functions over Rational Field in the Characteristic basis
              To:   Word Quasi-symmetric functions over Rational Field in the Monomial basis
"""

# ****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.global_options import GlobalOptions
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.shuffle import ShuffleProduct_overlapping, ShuffleProduct
from sage.rings.integer_ring import ZZ


class WQSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of `WQSym`.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis (must match one
      of the names that the basis can be constructed from `WQSym`)
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: TestSuite(M).run()  # long time
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         OrderedSetPartitions(),
                                         category=WQSymBases(alg, graded),
                                         bracket="", prefix=self._prefix)

    def _repr_term(self, osp):
        r"""
        Return a string representation of an element of WordQuasiSymmetricFunctions
        in the basis ``self``.

        TESTS::

            sage: M = WordQuasiSymmetricFunctions(QQ).M()
            sage: elt = M[[[1,2]]] * M[[[1]]]; elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
            sage: M.options.objects = "words"
            sage: elt
            M[1, 1, 2] + M[1, 1, 1] + M[2, 2, 1]
            sage: M.options._reset()
        """
        return self._prefix + self.options._dispatch(self, '_repr_', 'objects', osp)

    def _repr_compositions(self, osp):
        """
        Return a string representation of ``osp`` indexed by ordered set partitions.

        This method is called by ``self_repr_term``.

        EXAMPLES::

            sage: M = WordQuasiSymmetricFunctions(QQ).M()
            sage: elt = M[[[1,2]]] * M[[[1]]]; elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
            sage: M.options.display = "tight"
            sage: elt
            M[{1,2},{3}] + M[{1,2,3}] + M[{3},{1,2}]
            sage: M.options.display = "compact"
            sage: elt
            M[12.3] + M[123] + M[3.12]
            sage: osp = OrderedSetPartition([[2,4], [1,3,7],[5,6]])
            sage: M._repr_compositions(osp) == '[24.137.56]'
            True
            sage: M.options._reset(); elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
        """
        display = self.options.display
        disp = repr(osp)
        if display == 'tight':
            disp = disp.replace(", ", ",")
            return disp
        elif display == 'compact':
            disp = disp.replace("}, ", ".").replace("}", "").replace("{", "")
            return disp.replace(", ", "")
        else:
            # treat display as 'normal'
            return disp

    def _repr_words(self, osp):
        """
        Return a string representation of ``self`` indexed by packed words.

        This method is called by ``self_repr_term``.

        EXAMPLES::

            sage: M = WordQuasiSymmetricFunctions(QQ).M()
            sage: elt = M[[[1,2]]]*M[[[1]]]; elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
            sage: M.options.objects = "words"
            sage: elt
            M[1, 1, 2] + M[1, 1, 1] + M[2, 2, 1]
            sage: M.options.display = "tight"
            sage: elt
            M[1,1,2] + M[1,1,1] + M[2,2,1]
            sage: M.options.display = "compact"
            sage: elt
            M[112] + M[111] + M[221]
            sage: osp = OrderedSetPartition([[2,4], [1,3,7],[5,6]])
            sage: M._repr_words(osp) == '[2121332]'
            True
            sage: M.options._reset(); elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
        """
        display = self.options.display
        disp = repr(list(osp.to_packed_word()))
        if display == 'tight':
            return disp.replace(", ", ",")
        elif display == 'compact':
            return disp.replace(", ", "")
        else:
            # treat display as 'normal'
            return disp

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - word quasi-symmetric functions over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: M = algebras.WQSym(GF(7)).M(); M
            Word Quasi-symmetric functions over Finite Field of size 7 in the Monomial basis

        Elements of the word quasi-symmetric functions canonically coerce in::

            sage: x, y = M([[1]]), M([[2,1]])
            sage: M.coerce(x+y) == x+y
            True

        The word quasi-symmetric functions over `\ZZ` coerces in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: N = algebras.WQSym(ZZ).M()
            sage: Nx, Ny = N([[1]]), N([[2,1]])
            sage: z = M.coerce(Nx+Ny); z
            M[{1}] + M[{1, 2}]
            sage: z.parent() is M
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so word
        quasi-symmetric functions over `\GF{7}` does not coerce
        to the same algebra over `\ZZ`::

            sage: N.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Word Quasi-symmetric functions
             over Finite Field of size 7 in the Monomial basis to
             Word Quasi-symmetric functions over Integer Ring in the Monomial basis

        TESTS::

            sage: M = algebras.WQSym(ZZ).M()
            sage: N = algebras.WQSym(QQ).M()
            sage: M.has_coerce_map_from(N)
            False
            sage: N.has_coerce_map_from(M)
            True
            sage: M.has_coerce_map_from(QQ)
            False
            sage: N.has_coerce_map_from(QQ)
            True
            sage: M.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # word quasi-symmetric functions in the same variables
        # over any base that coerces in:
        if isinstance(R, WQSymBasis_abstract):
            if R.realization_of() == self.realization_of():
                return True
            if not self.base_ring().has_coerce_map_from(R.base_ring()):
                return False
            if self._basis_name == R._basis_name:  # The same basis
                def coerce_base_ring(self, x):
                    return self._from_dict(x.monomial_coefficients())
                return coerce_base_ring
            # Otherwise lift that basis up and then coerce over
            target = getattr(self.realization_of(), R._basis_name)()
            return self._coerce_map_via([target], R)
        return super(WQSymBasis_abstract, self)._coerce_map_from_(R)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: M.an_element()
            M[{1}] + 2*M[{1}, {2}]
        """
        return self([[1]]) + 2 * self([[1], [2]])

    def some_elements(self):
        """
        Return some elements of the word quasi-symmetric functions.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: M.some_elements()
            [M[], M[{1}], M[{1, 2}],
             M[{1}] + M[{1}, {2}],
             M[] + 1/2*M[{1}]]
        """
        u = self.one()
        o = self([[1]])
        s = self.base_ring().an_element()
        return [u, o, self([[1, 2]]), o + self([[1], [2]]), u + s * o]


class WordQuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The word quasi-symmetric functions.

    The ring of word quasi-symmetric functions can be defined as a
    subring of the ring of all bounded-degree noncommutative power
    series in countably many indeterminates (i.e., elements in
    `R \langle \langle x_1, x_2, x_3, \ldots \rangle \rangle` of bounded
    degree). Namely, consider words over the alphabet `\{1, 2, 3, \ldots\}`;
    every noncommutative power series is an infinite `R`-linear
    combination of these words.
    For each such word `w`, we define the *packing* of `w` to be the
    word `\operatorname{pack}(w)` that is obtained from `w` by replacing
    the smallest letter that appears in `w` by `1`, the second-smallest
    letter that appears in `w` by `2`, etc. (for example,
    `\operatorname{pack}(4112774) = 3112443`).
    A word `w` is said to be *packed* if `\operatorname{pack}(w) = w`.
    For each packed word `u`, we define the noncommutative power series
    `\mathbf{M}_u = \sum w`, where the sum ranges over all words `w`
    satisfying `\operatorname{pack}(w) = u`.
    The span of these power series `\mathbf{M}_u` is a subring of the
    ring of all noncommutative power series; it is called the ring of
    word quasi-symmetric functions, and is denoted by `WQSym`.

    For each nonnegative integer `n`, there is a bijection between
    packed words of length `n` and ordered set partitions of
    `\{1, 2, \ldots, n\}`. Under this bijection, a packed word
    `u = (u_1, u_2, \ldots, u_n)` of length `n` corresponds to the
    ordered set partition `P = (P_1, P_2, \ldots, P_k)` of
    `\{1, 2, \ldots, n\}` whose `i`-th part `P_i` (for each `i`) is the
    set of all `j \in \{1, 2, \ldots, n\}` such that `u_j = i`.

    The basis element `\mathbf{M}_u` is also denoted as `\mathbf{M}_P`
    in this situation. The basis `(\mathbf{M}_P)_P` is called the
    *Monomial basis* and is implemented as
    :class:`~sage.combinat.chas.wqsym.WordQuasiSymmetricFunctions.Monomial`.

    Other bases are the cone basis (aka C basis), the characteristic
    basis (aka X basis), the Q basis and the Phi basis.

    Bases of `WQSym` are implemented (internally) using ordered set
    partitions. However, the user may access specific basis vectors using
    either packed words or ordered set partitions. See the examples below,
    noting especially the section on ambiguities.

    `WQSym` is endowed with a connected graded Hopf algebra structure (see
    Section 2.2 of [NoThWi08]_, Section 1.1 of [FoiMal14]_ and
    Section 4.3.2 of [MeNoTh11]_) given by

    .. MATH::

        \Delta(\mathbf{M}_{(P_1,\ldots,P_{\ell})}) = \sum_{i=0}^{\ell}
            \mathbf{M}_{\operatorname{st}(P_1, \ldots, P_i)} \otimes
            \mathbf{M}_{\operatorname{st}(P_{i+1}, \ldots, P_{\ell})}.

    Here, for any ordered set partition `(Q_1, \ldots, Q_k)` of a
    finite set `Z` of integers, we let `\operatorname{st}(Q_1, \ldots, Q_k)`
    denote the set partition obtained from `Z` by replacing the smallest
    element appearing in it by `1`, the second-smallest element by `2`,
    and so on.

    A rule for multiplying elements of the monomial basis relies on the
    *quasi-shuffle product* of two ordered set partitions.
    The quasi-shuffle product `\Box` is given by
    :class:`~sage.combinat.shuffle.ShuffleProduct_overlapping` with the
    ``+`` operation in the overlapping of the shuffles being the
    union of the sets.  The product `\mathbf{M}_P \mathbf{M}_Q`
    for two ordered set partitions `P` and `Q` of `[n]` and `[m]`
    is then given by

    .. MATH::

        \mathbf{M}_P \mathbf{M}_Q
        = \sum_{R \in P \Box Q^+} \mathbf{M}_R ,

    where `Q^+` means `Q` with all numbers shifted upwards by `n`.

    Sometimes, `WQSym` is also denoted as `NCQSym`.

    REFERENCES:

    - [FoiMal14]_
    - [MeNoTh11]_
    - [NoThWi08]_
    - [BerZab05]_

    EXAMPLES:

    Constructing the algebra and its Monomial basis::

        sage: WQSym = algebras.WQSym(ZZ)
        sage: WQSym
        Word Quasi-symmetric functions over Integer Ring
        sage: M = WQSym.M()
        sage: M
        Word Quasi-symmetric functions over Integer Ring in the Monomial basis
        sage: M[[]]
        M[]

    Calling basis elements using packed words::

        sage: x = M[1,2,1]; x
        M[{1, 3}, {2}]
        sage: x == M[[1,2,1]] == M[Word([1,2,1])]
        True
        sage: y = M[1,1,2] - M[1,2,2]; y
        -M[{1}, {2, 3}] + M[{1, 2}, {3}]

    Calling basis elements using ordered set partitions::

        sage: z = M[[1,2,3],]; z
        M[{1, 2, 3}]
        sage: z == M[[[1,2,3]]] == M[OrderedSetPartition([[1,2,3]])]
        True
        sage: M[[1,2],[3]]
        M[{1, 2}, {3}]

    Note that expressions above are output in terms of ordered set partitions,
    even when input as packed words. Output as packed words can be achieved
    by modifying the global options. (See :meth:`OrderedSetPartitions.options`
    for further details.)::

        sage: M.options.objects = "words"
        sage: y
        -M[1, 2, 2] + M[1, 1, 2]
        sage: M.options.display = "compact"
        sage: y
        -M[122] + M[112]
        sage: z
        M[111]

    The options should be reset to display as ordered set partitions::

        sage: M.options._reset()
        sage: z
        M[{1, 2, 3}]

    Illustration of the Hopf algebra structure::

        sage: M[[2, 3], [5], [6], [4], [1]].coproduct()
        M[] # M[{2, 3}, {5}, {6}, {4}, {1}] + M[{1, 2}] # M[{3}, {4}, {2}, {1}]
         + M[{1, 2}, {3}] # M[{3}, {2}, {1}] + M[{1, 2}, {3}, {4}] # M[{2}, {1}]
         + M[{1, 2}, {4}, {5}, {3}] # M[{1}] + M[{2, 3}, {5}, {6}, {4}, {1}] # M[]
        sage: _ == M[5,1,1,4,2,3].coproduct()
        True
        sage: M[[1,1,1]] * M[[1,1,2]]   # packed words
        M[{1, 2, 3}, {4, 5}, {6}] + M[{1, 2, 3, 4, 5}, {6}]
         + M[{4, 5}, {1, 2, 3}, {6}] + M[{4, 5}, {1, 2, 3, 6}]
         + M[{4, 5}, {6}, {1, 2, 3}]
        sage: M[[1,2,3],].antipode()  # ordered set partition
        -M[{1, 2, 3}]
        sage: M[[1], [2], [3]].antipode()
        -M[{1, 2, 3}] - M[{2, 3}, {1}] - M[{3}, {1, 2}] - M[{3}, {2}, {1}]
        sage: x = M[[1],[2],[3]] + 3*M[[2],[1]]
        sage: x.counit()
        0
        sage: x.antipode()
        3*M[{1}, {2}] + 3*M[{1, 2}] - M[{1, 2, 3}] - M[{2, 3}, {1}]
         - M[{3}, {1, 2}] - M[{3}, {2}, {1}]

    .. rubric:: Ambiguities

    Some ambiguity arises when accessing basis vectors with the dictionary syntax,
    i.e., ``M[...]``. A common example is when referencing an ordered set partition
    with one part. For example, in the expression ``M[[1,2]]``, does ``[[1,2]]``
    refer to an ordered set partition or does ``[1,2]`` refer to a packed word?
    We choose the latter: if the received arguments do not behave like a tuple of
    iterables, then view them as describing a packed word. (In the running example,
    one argument is received, which behaves as a tuple of integers.) Here are a
    variety of ways to get the same basis vector::

        sage: x = M[1,1]; x
        M[{1, 2}]
        sage: x == M[[1,1]]  # treated as word
        True
        sage: x == M[[1,2],] == M[[[1,2]]]  # treated as ordered set partitions
        True

        sage: M[[1,3],[2]]  # treat as ordered set partition
        M[{1, 3}, {2}]
        sage: M[[1,3],[2]] == M[1,2,1]  # treat as word
        True

    TESTS::

        sage: M = WordQuasiSymmetricFunctions(QQ).M()
        sage: a = M[OrderedSetPartition([[1]])]
        sage: b = M[OrderedSetPartitions(1)([[1]])]
        sage: c = M[[1]]
        sage: a == b == c
        True

    .. TODO::

        - Dendriform structure.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.WQSym(QQ)
            sage: TestSuite(A).run()  # long time
        """
        category = HopfAlgebras(R).Graded().Connected()
        Parent.__init__(self, base=R, category=category.WithRealizations())

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.WQSym(QQ)  # indirect doctest
            Word Quasi-symmetric functions over Rational Field
        """
        return "Word Quasi-symmetric functions over {}".format(self.base_ring())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the `M`-basis).

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: WQSym.a_realization()
            Word Quasi-symmetric functions over Rational Field in the Monomial basis
        """
        return self.M()

    _shorthands = tuple(['M', 'X', 'C', 'Q', 'Phi'])

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the global options for bases of WordQuasiSymmetricFunctions.
        If no parameters are set, then the function returns a copy of the options
        dictionary.

        The ``options`` can be accessed as the method
        :obj:`WordQuasiSymmetricFunctions.options` of
        :class:`WordQuasiSymmetricFunctions` or of any associated basis.

        @OPTIONS@

        The ``'words'`` representation of a basis element of
        :class:`WordQuasiSymmetricFunctions`, indexed by an ordered
        set partition `A`, is the packed word associated to `A`.
        See :meth:`OrderedSetPartition.to_packed_word` for details.)

        EXAMPLES::

            sage: WQ = WordQuasiSymmetricFunctions(QQ)
            sage: M = WQ.M()
            sage: elt = M[[[1,2]]]*M[[[1]]]; elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
            sage: M.options.display = "tight"
            sage: elt
            M[{1,2},{3}] + M[{1,2,3}] + M[{3},{1,2}]
            sage: M.options.display = "compact"
            sage: elt
            M[12.3] + M[123] + M[3.12]
            sage: WQ.options._reset()
            sage: M.options.objects = "words"
            sage: elt
            M[1, 1, 2] + M[1, 1, 1] + M[2, 2, 1]
            sage: M.options.display = "tight"
            sage: elt
            M[1,1,2] + M[1,1,1] + M[2,2,1]
            sage: WQ.options.display = "compact"
            sage: elt
            M[112] + M[111] + M[221]
            sage: M.options._reset()
            sage: elt
            M[{1, 2}, {3}] + M[{1, 2, 3}] + M[{3}, {1, 2}]
        """
        NAME = 'WordQuasiSymmetricFunctions element'
        module = 'sage.combinat.chas.wqsym'
        option_class = 'WordQuasiSymmetricFunctions'
        objects = dict(default="compositions",
                       description='Specifies how basis elements of WordQuasiSymmetricFunctions should be indexed',
                       values=dict(compositions="Indexing the basis by ordered set partitions",
                                   words="Indexing the basis by packed words"),
                       case_sensitive=False)
        display = dict(default="normal",
                       description='Specifies how basis elements of WordQuasiSymmetricFunctions should be printed',
                       values=dict(normal="Using the normal representation",
                                   tight="Dropping spaces after commas",
                                   compact="Using a severely compacted representation"),
                       case_sensitive=False)

    class Monomial(WQSymBasis_abstract):
        r"""
        The Monomial basis of `WQSym`.

        The family `(\mathbf{M}_u)`, as defined in
        :class:`~sage.combinat.chas.wqsym.WordQuasiSymmetricFunctions`
        with `u` ranging over all packed words, is a basis for the
        free `R`-module `WQSym` and called the *Monomial basis*.
        Here it is labelled using ordered set partitions.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: M = WQSym.M(); M
            Word Quasi-symmetric functions over Rational Field in the Monomial basis
            sage: sorted(M.basis(2))
            [M[{1, 2}], M[{2}, {1}], M[{1}, {2}]]
        """
        _prefix = "M"
        _basis_name = "Monomial"

        def product_on_basis(self, x, y):
            r"""
            Return the (associative) `*` product of the basis elements
            of ``self`` indexed by the ordered set partitions `x` and
            `y`.

            This is the shifted quasi-shuffle product of `x` and `y`.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
                sage: x = OrderedSetPartition([[1],[2,3]])
                sage: y = OrderedSetPartition([[1,2]])
                sage: z = OrderedSetPartition([[1,2],[3]])
                sage: A.product_on_basis(x, y)
                M[{1}, {2, 3}, {4, 5}] + M[{1}, {2, 3, 4, 5}]
                 + M[{1}, {4, 5}, {2, 3}] + M[{1, 4, 5}, {2, 3}]
                 + M[{4, 5}, {1}, {2, 3}]
                sage: A.product_on_basis(x, z)
                M[{1}, {2, 3}, {4, 5}, {6}] + M[{1}, {2, 3, 4, 5}, {6}]
                 + M[{1}, {4, 5}, {2, 3}, {6}] + M[{1}, {4, 5}, {2, 3, 6}]
                 + M[{1}, {4, 5}, {6}, {2, 3}] + M[{1, 4, 5}, {2, 3}, {6}]
                 + M[{1, 4, 5}, {2, 3, 6}] + M[{1, 4, 5}, {6}, {2, 3}]
                 + M[{4, 5}, {1}, {2, 3}, {6}] + M[{4, 5}, {1}, {2, 3, 6}]
                 + M[{4, 5}, {1}, {6}, {2, 3}] + M[{4, 5}, {1, 6}, {2, 3}]
                 + M[{4, 5}, {6}, {1}, {2, 3}]
                sage: A.product_on_basis(y, y)
                M[{1, 2}, {3, 4}] + M[{1, 2, 3, 4}] + M[{3, 4}, {1, 2}]

            TESTS::

                sage: one = OrderedSetPartition([])
                sage: all(A.product_on_basis(one, z) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
                sage: all(A.product_on_basis(z, one) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
            """
            K = self.basis().keys()
            if not x:
                return self.monomial(y)
            m = max(max(part) for part in x)  # The degree of x
            x = [set(part) for part in x]
            yshift = [[val + m for val in part] for part in y]

            def union(X, Y):
                return X.union(Y)
            return self.sum_of_monomials(ShuffleProduct_overlapping(x, yshift,
                                                                    K, union))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of ``self`` on the basis element
            indexed by the ordered set partition ``x``.

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()

                sage: M.coproduct(M.one())  # indirect doctest
                M[] # M[]
                sage: M.coproduct( M([[1]]) )  # indirect doctest
                M[] # M[{1}] + M[{1}] # M[]
                sage: M.coproduct( M([[1,2]]) )
                M[] # M[{1, 2}] + M[{1, 2}] # M[]
                sage: M.coproduct( M([[1], [2]]) )
                M[] # M[{1}, {2}] + M[{1}] # M[{1}] + M[{1}, {2}] # M[]
            """
            if not len(x):
                return self.one().tensor(self.one())
            K = self.indices()

            def standardize(P):  # standardize an ordered set partition
                base = sorted(sum((list(part) for part in P), []))
                # base is the ground set of P, as a sorted list.
                d = {val: i + 1 for i, val in enumerate(base)}
                # d is the unique order isomorphism from base to
                # {1, 2, ..., |base|} (encoded as dict).
                return K([[d[x] for x in part] for part in P])
            T = self.tensor_square()
            return T.sum_of_monomials((standardize(x[:i]), standardize(x[i:]))
                                      for i in range(len(x) + 1))

    M = Monomial

    class Characteristic(WQSymBasis_abstract):
        r"""
        The Characteristic basis of `WQSym`.

        The *Characteristic basis* is a graded basis `(X_P)` of `WQSym`,
        indexed by ordered set partitions `P`. It is defined by

        .. MATH::

            X_P = (-1)^{\ell(P)} \mathbf{M}_P ,

        where `(\mathbf{M}_P)_P` denotes the Monomial basis,
        and where `\ell(P)` denotes the number of blocks in an ordered
        set partition `P`.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: X = WQSym.X(); X
            Word Quasi-symmetric functions over Rational Field in the Characteristic basis

            sage: X[[[1,2,3]]] * X[[1,2],[3]]
            X[{1, 2, 3}, {4, 5}, {6}] - X[{1, 2, 3, 4, 5}, {6}]
             + X[{4, 5}, {1, 2, 3}, {6}] - X[{4, 5}, {1, 2, 3, 6}]
             + X[{4, 5}, {6}, {1, 2, 3}]

            sage: X[[1, 4], [3], [2]].coproduct()
            X[] # X[{1, 4}, {3}, {2}] + X[{1, 2}] # X[{2}, {1}]
             + X[{1, 3}, {2}] # X[{1}] + X[{1, 4}, {3}, {2}] # X[]

            sage: M = WQSym.M()
            sage: M(X[[1, 2, 3],])
            -M[{1, 2, 3}]
            sage: M(X[[1, 3], [2]])
            M[{1, 3}, {2}]
            sage: X(M[[1, 2, 3],])
            -X[{1, 2, 3}]
            sage: X(M[[1, 3], [2]])
            X[{1, 3}, {2}]
        """
        _prefix = "X"
        _basis_name = "Characteristic"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: X = algebras.WQSym(QQ).X()
                sage: TestSuite(X).run()  # long time
            """
            WQSymBasis_abstract.__init__(self, alg)

            M = self.realization_of().M()
            mone = -self.base_ring().one()

            def sgn(P):
                return mone**len(P)
            self.module_morphism(codomain=M, diagonal=sgn).register_as_coercion()
            M.module_morphism(codomain=self, diagonal=sgn).register_as_coercion()

        class Element(WQSymBasis_abstract.Element):
            def algebraic_complement(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the algebraic complement involution.

                See
                :meth:`WQSymBases.ElementMethods.algebraic_complement`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`coalgebraic_complement`, :meth:`star_involution`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: X = WQSym.X()
                    sage: X[[1,2],[5,6],[3,4]].algebraic_complement()
                    X[{3, 4}, {5, 6}, {1, 2}]
                    sage: X[[3], [1, 2], [4]].algebraic_complement()
                    X[{4}, {1, 2}, {3}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(X[A]).algebraic_complement() == M(X[A].algebraic_complement())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.algebraic_complement doc
                # for the formula we're using here.
                Q = self.parent()
                OSPs = Q.basis().keys()
                return Q._from_dict({OSPs(A.reversed()): c for (A, c) in self},
                                    remove_zeros=False)

            def coalgebraic_complement(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the coalgebraic complement involution.

                See
                :meth:`WQSymBases.ElementMethods.coalgebraic_complement`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`algebraic_complement`, :meth:`star_involution`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: X = WQSym.X()
                    sage: X[[1,2],[5,6],[3,4]].coalgebraic_complement()
                    X[{5, 6}, {1, 2}, {3, 4}]
                    sage: X[[3], [1, 2], [4]].coalgebraic_complement()
                    X[{2}, {3, 4}, {1}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(X[A]).coalgebraic_complement()
                    ....:     == M(X[A].coalgebraic_complement())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.coalgebraic_complement doc
                # for the formula we're using here.
                Q = self.parent()
                OSPs = Q.basis().keys()
                return Q._from_dict({OSPs(A.complement()): c for (A, c) in self},
                                    remove_zeros=False)

            def star_involution(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the star involution.

                See
                :meth:`WQSymBases.ElementMethods.star_involution`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`algebraic_complement`, :meth:`coalgebraic_complement`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: X = WQSym.X()
                    sage: X[[1,2],[5,6],[3,4]].star_involution()
                    X[{3, 4}, {1, 2}, {5, 6}]
                    sage: X[[3], [1, 2], [4]].star_involution()
                    X[{1}, {3, 4}, {2}]

                TESTS:

                    sage: M = WQSym.M()
                    sage: all(M(X[A]).star_involution() == M(X[A].star_involution())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.star_involution doc
                # for the formula we're using here.
                Q = self.parent()
                OSPs = Q.basis().keys()
                return Q._from_dict({OSPs(A.complement().reversed()): c for (A, c) in self},
                                    remove_zeros=False)

    X = Characteristic

    class Cone(WQSymBasis_abstract):
        r"""
        The Cone basis of `WQSym`.

        Let `(X_P)_P` denote the Characteristic basis of `WQSym`.
        Denote the quasi-shuffle of two ordered set partitions `A` and
        `B` by `A \Box B`. For an ordered set partition
        `P = (P_1, \ldots, P_{\ell})`, we form a list of ordered set
        partitions `[P] := (P'_1, \ldots, P'_k)` as follows.
        Define a strictly decreasing sequence of integers
        `\ell + 1 = i_0 > i_1 > \cdots > i_k = 1` recursively by
        requiring that `\min P_{i_j} \leq \min P_a` for all `a < i_{j-1}`.
        Set `P'_j = (P_{i_j}, \ldots, P_{i_{j-1}-1})`.

        The *Cone basis* `(C_P)_P` is defined by

        .. MATH::

            C_P = \sum_Q X_Q,

        where the sum is over all elements `Q` of the quasi-shuffle
        product `P'_1 \Box P'_2 \Box \cdots \Box P'_k` with
        `[P] = (P'_1, \ldots, P'_k)`.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: C = WQSym.C()
            sage: C
            Word Quasi-symmetric functions over Rational Field in the Cone basis

            sage: X = WQSym.X()
            sage: X(C[[2,3],[1,4]])
            X[{1, 2, 3, 4}] + X[{1, 4}, {2, 3}] + X[{2, 3}, {1, 4}]
            sage: X(C[[1,4],[2,3]])
            X[{1, 4}, {2, 3}]
            sage: X(C[[2,3],[1],[4]])
            X[{1}, {2, 3}, {4}] + X[{1}, {2, 3, 4}] + X[{1}, {4}, {2, 3}]
             + X[{1, 2, 3}, {4}] + X[{2, 3}, {1}, {4}]
            sage: X(C[[3], [2, 5], [1, 4]])
            X[{1, 2, 3, 4, 5}] + X[{1, 2, 4, 5}, {3}] + X[{1, 3, 4}, {2, 5}]
             + X[{1, 4}, {2, 3, 5}] + X[{1, 4}, {2, 5}, {3}]
             + X[{1, 4}, {3}, {2, 5}] + X[{2, 3, 5}, {1, 4}]
             + X[{2, 5}, {1, 3, 4}] + X[{2, 5}, {1, 4}, {3}]
             + X[{2, 5}, {3}, {1, 4}] + X[{3}, {1, 2, 4, 5}]
             + X[{3}, {1, 4}, {2, 5}] + X[{3}, {2, 5}, {1, 4}]
            sage: C(X[[2,3],[1,4]])
            -C[{1, 2, 3, 4}] - C[{1, 4}, {2, 3}] + C[{2, 3}, {1, 4}]

        REFERENCES:

        - Section 4 of [Early2017]_

        .. TODO::

            Experiments suggest that :meth:`algebraic_complement`,
            :meth:`coalgebraic_complement`, and :meth:`star_involution`
            should have reasonable formulas on the C basis; at least
            the coefficients of the outputs on any element of the C
            basis seem to be always `0, 1, -1`.
            Is this true? What is the formula?
        """
        _prefix = "C"
        _basis_name = "Cone"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: C = algebras.WQSym(QQ).C()
                sage: TestSuite(C).run()  # long time
            """
            WQSymBasis_abstract.__init__(self, alg)

            X = self.realization_of().X()
            phi = self.module_morphism(self._C_to_X, codomain=X, unitriangular="upper")
            phi.register_as_coercion()
            inv_phi = ~phi
            inv_phi.register_as_coercion()
            # We need to explicitly construct the coercion to/from M via X
            #   as otherwise, when another basis B is created before X, the
            #   coercion is attempted to be built via B, which is cannot do.
            #   So the coercion map returned is the default one that calls
            #   the _element_constructor_.
            #   This is only a problem because X is not the default
            #   a_realization(), which is M, and the coercions are always
            #   first attempted through M to another basis. -- TS
            M = self.realization_of().M()
            M.register_coercion(M.coerce_map_from(X) * phi)
            self.register_coercion(inv_phi * X.coerce_map_from(M))

        def some_elements(self):
            """
            Return some elements of the word quasi-symmetric functions
            in the Cone basis.

            EXAMPLES::

                sage: C = algebras.WQSym(QQ).C()
                sage: C.some_elements()
                [C[], C[{1}], C[{1, 2}], C[] + 1/2*C[{1}]]
            """
            u = self.one()
            o = self([[1]])
            s = self.base_ring().an_element()
            return [u, o, self([[1, 2]]), u + s * o]

        def _C_to_X(self, P):
            """
            Return the image of the basis element of ``self`` indexed
            by ``P`` in the Characteristic basis.

            EXAMPLES::

                sage: C = algebras.WQSym(QQ).C()
                sage: OSP = C.basis().keys()
                sage: C._C_to_X(OSP([[2,3],[1,4]]))
                X[{1, 2, 3, 4}] + X[{1, 4}, {2, 3}] + X[{2, 3}, {1, 4}]
            """
            X = self.realization_of().X()
            if not P:
                return X.one()

            OSP = self.basis().keys()

            # Convert to standard set of ordered set partitions
            temp = list(P)
            data = []
            while temp:
                i = min(min(X) for X in temp)
                for j, A in enumerate(temp):
                    if i in A:
                        data.append(OSP(temp[j:]))
                        temp = temp[:j]
                        break

            # Perform the quasi-shuffle product
            cur = {data[0]: 1}
            for B in data[1:]:
                ret = {}
                for A in cur:
                    for C in ShuffleProduct_overlapping(A, B, element_constructor=OSP):
                        if C in ret:
                            ret[C] += cur[A]
                        else:
                            ret[C] = cur[A]
                cur = ret

            # Return the result in the X basis
            return X._from_dict(cur, coerce=True)

    C = Cone

    class StronglyCoarser(WQSymBasis_abstract):
        r"""
        The Q basis of `WQSym`.

        We define a partial order `\leq` on the set of all ordered
        set partitions as follows: `A \leq B` if and only if
        `A` is strongly finer than `B` (see
        :meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.is_strongly_finer`
        for a definition of this).

        The *Q basis* `(Q_P)_P` is a basis of `WQSym` indexed by ordered
        set partitions, and is defined by

        .. MATH::

            Q_P = \sum \mathbf{M}_W,

        where the sum is over ordered set partitions `W` satisfying
        `P \leq W`.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: M = WQSym.M(); Q = WQSym.Q()
            sage: Q
            Word Quasi-symmetric functions over Rational Field in the Q basis

            sage: Q(M[[2,3],[1,4]])
            Q[{2, 3}, {1, 4}]
            sage: Q(M[[1,2],[3,4]])
            Q[{1, 2}, {3, 4}] - Q[{1, 2, 3, 4}]
            sage: M(Q[[1,2],[3,4]])
            M[{1, 2}, {3, 4}] + M[{1, 2, 3, 4}]
            sage: M(Q[[2,3],[1],[4]])
            M[{2, 3}, {1}, {4}] + M[{2, 3}, {1, 4}]
            sage: M(Q[[3], [2, 5], [1, 4]])
            M[{3}, {2, 5}, {1, 4}]
            sage: M(Q[[1, 4], [2, 3], [5], [6]])
            M[{1, 4}, {2, 3}, {5}, {6}] + M[{1, 4}, {2, 3}, {5, 6}]
             + M[{1, 4}, {2, 3, 5}, {6}] + M[{1, 4}, {2, 3, 5, 6}]

            sage: Q[[1, 3], [2]] * Q[[1], [2]]
            Q[{1, 3}, {2}, {4}, {5}] + Q[{1, 3}, {4}, {2}, {5}]
             + Q[{1, 3}, {4}, {5}, {2}] + Q[{4}, {1, 3}, {2}, {5}]
             + Q[{4}, {1, 3}, {5}, {2}] + Q[{4}, {5}, {1, 3}, {2}]

            sage: Q[[1, 3], [2]].coproduct()
            Q[] # Q[{1, 3}, {2}] + Q[{1, 2}] # Q[{1}] + Q[{1, 3}, {2}] # Q[]

        REFERENCES:

        - Section 6 of [BerZab05]_
        """
        _prefix = "Q"
        _basis_name = "Q"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: Q = algebras.WQSym(QQ).Q()
                sage: TestSuite(Q).run()  # long time
            """
            WQSymBasis_abstract.__init__(self, alg)

            M = self.realization_of().M()
            phi = self.module_morphism(self._Q_to_M, codomain=M, unitriangular="lower")
            phi.register_as_coercion()
            phi_inv = M.module_morphism(self._M_to_Q, codomain=self, unitriangular="lower")
            phi_inv.register_as_coercion()

        def some_elements(self):
            """
            Return some elements of the word quasi-symmetric functions
            in the Q basis.

            EXAMPLES::

                sage: Q = algebras.WQSym(QQ).Q()
                sage: Q.some_elements()
                [Q[], Q[{1}], Q[{1, 2}], Q[] + 1/2*Q[{1}]]
            """
            u = self.one()
            o = self([[1]])
            s = self.base_ring().an_element()
            return [u, o, self([[1, 2]]), u + s * o]

        def _Q_to_M(self, P):
            """
            Return the image of the basis element of ``self`` indexed
            by ``P`` in the Monomial basis.

            EXAMPLES::

                sage: Q = algebras.WQSym(QQ).Q()
                sage: OSP = Q.basis().keys()
                sage: Q._Q_to_M(OSP([[2,3],[1,4]]))
                M[{2, 3}, {1, 4}]
                sage: Q._Q_to_M(OSP([[1,2],[3,4]]))
                M[{1, 2}, {3, 4}] + M[{1, 2, 3, 4}]
            """
            M = self.realization_of().M()
            if not P:
                return M.one()

            OSP = self.basis().keys()
            R = M.base_ring()
            one = R.one()
            return M._from_dict({OSP(G): one for G in P.strongly_fatter()},
                                coerce=False)

        def _M_to_Q(self, P):
            """
            Return the image of the basis element of the monomial
            basis indexed by ``P`` in the Q basis ``self``.

            EXAMPLES::

                sage: Q = algebras.WQSym(QQ).Q()
                sage: M = algebras.WQSym(QQ).M()
                sage: OSP = Q.basis().keys()
                sage: Q._M_to_Q(OSP([[2,3],[1,4]]))
                Q[{2, 3}, {1, 4}]
                sage: Q._M_to_Q(OSP([[1,2],[3,4]]))
                Q[{1, 2}, {3, 4}] - Q[{1, 2, 3, 4}]

            TESTS::

                sage: Q = algebras.WQSym(QQ).Q()
                sage: M = algebras.WQSym(QQ).M()
                sage: OSP4 = OrderedSetPartitions(4)
                sage: all(M(Q(M[P])) == M[P] for P in OSP4) # long time
                True
                sage: all(Q(M(Q[P])) == Q[P] for P in OSP4) # long time
                True
            """
            Q = self
            if not P:
                return Q.one()

            OSP = self.basis().keys()
            R = self.base_ring()
            one = R.one()
            lenP = len(P)

            def sign(R):
                # the coefficient with which another
                # ordered set partition will appear
                if len(R) % 2 == lenP % 2:
                    return one
                return -one
            return Q._from_dict({OSP(G): sign(G) for G in P.strongly_fatter()},
                                coerce=False)

        def product_on_basis(self, x, y):
            r"""
            Return the (associative) `*` product of the basis elements
            of the Q basis ``self`` indexed by the ordered set partitions
            `x` and `y`.

            This is the shifted shuffle product of `x` and `y`.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).Q()
                sage: x = OrderedSetPartition([[1],[2,3]])
                sage: y = OrderedSetPartition([[1,2]])
                sage: z = OrderedSetPartition([[1,2],[3]])
                sage: A.product_on_basis(x, y)
                Q[{1}, {2, 3}, {4, 5}] + Q[{1}, {4, 5}, {2, 3}]
                 + Q[{4, 5}, {1}, {2, 3}]
                sage: A.product_on_basis(x, z)
                Q[{1}, {2, 3}, {4, 5}, {6}] + Q[{1}, {4, 5}, {2, 3}, {6}]
                 + Q[{1}, {4, 5}, {6}, {2, 3}] + Q[{4, 5}, {1}, {2, 3}, {6}]
                 + Q[{4, 5}, {1}, {6}, {2, 3}] + Q[{4, 5}, {6}, {1}, {2, 3}]
                sage: A.product_on_basis(y, y)
                Q[{1, 2}, {3, 4}] + Q[{3, 4}, {1, 2}]

            TESTS::

                sage: one = OrderedSetPartition([])
                sage: all(A.product_on_basis(one, z) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
                sage: all(A.product_on_basis(z, one) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
            """
            K = self.basis().keys()
            if not x:
                return self.monomial(y)
            m = max(max(part) for part in x)  # The degree of x
            x = [set(part) for part in x]
            yshift = [[val + m for val in part] for part in y]
            return self.sum_of_monomials(ShuffleProduct(x, yshift, K))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of ``self`` on the basis element
            indexed by the ordered set partition ``x``.

            EXAMPLES::

                sage: Q = algebras.WQSym(QQ).Q()

                sage: Q.coproduct(Q.one())  # indirect doctest
                Q[] # Q[]
                sage: Q.coproduct( Q([[1]]) )  # indirect doctest
                Q[] # Q[{1}] + Q[{1}] # Q[]
                sage: Q.coproduct( Q([[1,2]]) )
                Q[] # Q[{1, 2}] + Q[{1, 2}] # Q[]
                sage: Q.coproduct( Q([[1], [2]]) )
                Q[] # Q[{1}, {2}] + Q[{1}] # Q[{1}] + Q[{1}, {2}] # Q[]
                sage: Q[[1,2],[3],[4]].coproduct()
                Q[] # Q[{1, 2}, {3}, {4}] + Q[{1, 2}] # Q[{1}, {2}]
                 + Q[{1, 2}, {3}] # Q[{1}] + Q[{1, 2}, {3}, {4}] # Q[]
            """
            # The coproduct on the Q basis satisfies the same formula
            # as on the M basis. This is easily derived from the
            # formula on the M basis.
            if not len(x):
                return self.one().tensor(self.one())
            K = self.indices()

            def standardize(P):  # standardize an ordered set partition
                base = sorted(sum((list(part) for part in P), []))
                # base is the ground set of P, as a sorted list.
                d = {val: i + 1 for i, val in enumerate(base)}
                # d is the unique order isomorphism from base to
                # {1, 2, ..., |base|} (encoded as dict).
                return K([[d[x] for x in part] for part in P])
            T = self.tensor_square()
            return T.sum_of_monomials((standardize(x[:i]), standardize(x[i:]))
                                      for i in range(len(x) + 1))

        class Element(WQSymBasis_abstract.Element):
            def algebraic_complement(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the algebraic complement involution.

                See
                :meth:`WQSymBases.ElementMethods.algebraic_complement`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`coalgebraic_complement`, :meth:`star_involution`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: Q = WQSym.Q()
                    sage: Q[[1,2],[5,6],[3,4]].algebraic_complement()
                    Q[{3, 4}, {1, 2, 5, 6}] + Q[{3, 4}, {5, 6}, {1, 2}]
                     - Q[{3, 4, 5, 6}, {1, 2}]
                    sage: Q[[3], [1, 2], [4]].algebraic_complement()
                    Q[{1, 2, 4}, {3}] + Q[{4}, {1, 2}, {3}] - Q[{4}, {1, 2, 3}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(Q[A]).algebraic_complement()  # long time
                    ....:     == M(Q[A].algebraic_complement())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.algebraic_complement doc
                # for the formula we're using here.
                BR = self.base_ring()
                one = BR.one()
                mine = -one
                Q = self.parent()
                OSPs = Q.basis().keys()
                from sage.data_structures.blas_dict import linear_combination

                def img(A):
                    # The image of the basis element Q[A], written as a
                    # dictionary (of its coordinates in the Q-basis).
                    Rs = [Rr.reversed() for Rr in A.strongly_fatter()]
                    return {OSPs(P): (one if (len(R) % 2 == len(P) % 2)
                                      else mine)
                            for R in Rs for P in R.strongly_fatter()}
                return Q._from_dict(linear_combination((img(A), c) for (A, c) in self))

            def coalgebraic_complement(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the coalgebraic complement involution.

                See
                :meth:`WQSymBases.ElementMethods.coalgebraic_complement`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`algebraic_complement`, :meth:`star_involution`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: Q = WQSym.Q()
                    sage: Q[[1,2],[5,6],[3,4]].coalgebraic_complement()
                    Q[{1, 2, 5, 6}, {3, 4}] + Q[{5, 6}, {1, 2}, {3, 4}] - Q[{5, 6}, {1, 2, 3, 4}]
                    sage: Q[[3], [1, 2], [4]].coalgebraic_complement()
                    Q[{2}, {1, 3, 4}] + Q[{2}, {3, 4}, {1}] - Q[{2, 3, 4}, {1}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(Q[A]).coalgebraic_complement()  # long time
                    ....:     == M(Q[A].coalgebraic_complement())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.coalgebraic_complement doc
                # for the formula we're using here.
                BR = self.base_ring()
                one = BR.one()
                mine = -one
                Q = self.parent()
                OSPs = Q.basis().keys()
                from sage.data_structures.blas_dict import linear_combination

                def img(A):
                    # The image of the basis element Q[A], written as a
                    # dictionary (of its coordinates in the Q-basis).
                    Rs = [Rr.complement() for Rr in A.strongly_fatter()]
                    return {OSPs(P): (one if (len(R) % 2 == len(P) % 2)
                                      else mine)
                            for R in Rs for P in R.strongly_fatter()}
                return Q._from_dict(linear_combination((img(A), c) for (A, c) in self))

            def star_involution(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the star involution.

                See
                :meth:`WQSymBases.ElementMethods.star_involution`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`algebraic_complement`, :meth:`coalgebraic_complement`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: Q = WQSym.Q()
                    sage: Q[[1,2],[5,6],[3,4]].star_involution()
                    Q[{3, 4}, {1, 2}, {5, 6}]
                    sage: Q[[3], [1, 2], [4]].star_involution()
                    Q[{1}, {3, 4}, {2}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(Q[A]).star_involution() == M(Q[A].star_involution())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.star_involution doc
                # for the formula we're using here.
                Q = self.parent()
                OSPs = Q.basis().keys()
                return Q._from_dict({OSPs(A.complement().reversed()): c for (A, c) in self},
                                    remove_zeros=False)

    Q = StronglyCoarser

    class StronglyFiner(WQSymBasis_abstract):
        r"""
        The Phi basis of `WQSym`.

        We define a partial order `\leq` on the set of all ordered
        set partitions as follows: `A \leq B` if and only if
        `A` is strongly finer than `B` (see
        :meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.is_strongly_finer`
        for a definition of this).

        The *Phi basis* `(\Phi_P)_P` is a basis of `WQSym` indexed by ordered
        set partitions, and is defined by

        .. MATH::

            \Phi_P = \sum \mathbf{M}_W,

        where the sum is over ordered set partitions `W` satisfying
        `W \leq P`.

        Novelli and Thibon introduced this basis in [NovThi06]_
        Section 2.7.2, and called it the quasi-ribbon basis.
        It later reappeared in [MeNoTh11]_ Section 4.3.2.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: M = WQSym.M(); Phi = WQSym.Phi()
            sage: Phi
            Word Quasi-symmetric functions over Rational Field in the Phi basis

            sage: Phi(M[[2,3],[1,4]])
            Phi[{2}, {3}, {1}, {4}] - Phi[{2}, {3}, {1, 4}]
             - Phi[{2, 3}, {1}, {4}] + Phi[{2, 3}, {1, 4}]
            sage: Phi(M[[1,2],[3,4]])
            Phi[{1}, {2}, {3}, {4}] - Phi[{1}, {2}, {3, 4}]
             - Phi[{1, 2}, {3}, {4}] + Phi[{1, 2}, {3, 4}]
            sage: M(Phi[[1,2],[3,4]])
            M[{1}, {2}, {3}, {4}] + M[{1}, {2}, {3, 4}]
             + M[{1, 2}, {3}, {4}] + M[{1, 2}, {3, 4}]
            sage: M(Phi[[2,3],[1],[4]])
            M[{2}, {3}, {1}, {4}] + M[{2, 3}, {1}, {4}]
            sage: M(Phi[[3], [2, 5], [1, 4]])
            M[{3}, {2}, {5}, {1}, {4}] + M[{3}, {2}, {5}, {1, 4}]
             + M[{3}, {2, 5}, {1}, {4}] + M[{3}, {2, 5}, {1, 4}]
            sage: M(Phi[[1, 4], [2, 3], [5], [6]])
            M[{1}, {4}, {2}, {3}, {5}, {6}] + M[{1}, {4}, {2, 3}, {5}, {6}]
             + M[{1, 4}, {2}, {3}, {5}, {6}] + M[{1, 4}, {2, 3}, {5}, {6}]

            sage: Phi[[1],] * Phi[[1, 3], [2]]
            Phi[{1, 2, 4}, {3}] + Phi[{2}, {1, 4}, {3}]
             + Phi[{2, 4}, {1, 3}] + Phi[{2, 4}, {3}, {1}]
            sage: Phi[[3, 5], [1, 4], [2]].coproduct()
            Phi[] # Phi[{3, 5}, {1, 4}, {2}]
             + Phi[{1}] # Phi[{4}, {1, 3}, {2}]
             + Phi[{1, 2}] # Phi[{1, 3}, {2}]
             + Phi[{2, 3}, {1}] # Phi[{2}, {1}]
             + Phi[{2, 4}, {1, 3}] # Phi[{1}]
             + Phi[{3, 5}, {1, 4}, {2}] # Phi[]

        REFERENCES:

        - Section 2.7.2 of [NovThi06]_
        """
        _prefix = "Phi"
        _basis_name = "Phi"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: Phi = algebras.WQSym(QQ).Phi()
                sage: TestSuite(Phi).run()  # long time
            """
            WQSymBasis_abstract.__init__(self, alg)

            M = self.realization_of().M()
            phi = self.module_morphism(self._Phi_to_M, codomain=M, unitriangular="lower")
            phi.register_as_coercion()
            phi_inv = M.module_morphism(self._M_to_Phi, codomain=self, unitriangular="lower")
            phi_inv.register_as_coercion()

        def some_elements(self):
            """
            Return some elements of the word quasi-symmetric functions
            in the Phi basis.

            EXAMPLES::

                sage: Phi = algebras.WQSym(QQ).Phi()
                sage: Phi.some_elements()
                [Phi[], Phi[{1}], Phi[{1, 2}], Phi[] + 1/2*Phi[{1}]]
            """
            u = self.one()
            o = self([[1]])
            s = self.base_ring().an_element()
            return [u, o, self([[1, 2]]), u + s * o]

        def _Phi_to_M(self, P):
            """
            Return the image of the basis element of ``self`` indexed
            by ``P`` in the Monomial basis.

            EXAMPLES::

                sage: Phi = algebras.WQSym(QQ).Phi()
                sage: OSP = Phi.basis().keys()
                sage: Phi._Phi_to_M(OSP([[2,3],[1,4]]))
                M[{2}, {3}, {1}, {4}] + M[{2}, {3}, {1, 4}]
                 + M[{2, 3}, {1}, {4}] + M[{2, 3}, {1, 4}]
                sage: Phi._Phi_to_M(OSP([[1,2],[3,4]]))
                M[{1}, {2}, {3}, {4}] + M[{1}, {2}, {3, 4}]
                 + M[{1, 2}, {3}, {4}] + M[{1, 2}, {3, 4}]
            """
            M = self.realization_of().M()
            if not P:
                return M.one()

            OSP = self.basis().keys()
            R = M.base_ring()
            one = R.one()
            return M._from_dict({OSP(G): one for G in P.strongly_finer()},
                                coerce=False)

        def _M_to_Phi(self, P):
            """
            Return the image of the basis element of the monomial
            basis indexed by ``P`` in the Phi basis ``self``.

            EXAMPLES::

                sage: Phi = algebras.WQSym(QQ).Phi()
                sage: M = algebras.WQSym(QQ).M()
                sage: OSP = Phi.basis().keys()
                sage: Phi._M_to_Phi(OSP([[2,3],[1,4]]))
                Phi[{2}, {3}, {1}, {4}] - Phi[{2}, {3}, {1, 4}]
                 - Phi[{2, 3}, {1}, {4}] + Phi[{2, 3}, {1, 4}]
                sage: Phi._M_to_Phi(OSP([[1,2],[3,4]]))
                Phi[{1}, {2}, {3}, {4}] - Phi[{1}, {2}, {3, 4}]
                 - Phi[{1, 2}, {3}, {4}] + Phi[{1, 2}, {3, 4}]

            TESTS::

                sage: Phi = algebras.WQSym(QQ).Phi()
                sage: M = algebras.WQSym(QQ).M()
                sage: OSP4 = OrderedSetPartitions(4)
                sage: all(M(Phi(M[P])) == M[P] for P in OSP4) # long time
                True
                sage: all(Phi(M(Phi[P])) == Phi[P] for P in OSP4) # long time
                True
            """
            Phi = self
            if not P:
                return Phi.one()

            OSP = self.basis().keys()
            R = self.base_ring()
            one = R.one()
            lenP = len(P)

            def sign(R):
                # the coefficient with which another
                # ordered set partition will appear
                if len(R) % 2 == lenP % 2:
                    return one
                return -one
            return Phi._from_dict({OSP(G): sign(G) for G in P.strongly_finer()},
                                  coerce=False)

        def product_on_basis(self, x, y):
            r"""
            Return the (associative) `*` product of the basis elements
            of the Phi basis ``self`` indexed by the ordered set partitions
            `x` and `y`.

            This is obtained by the following algorithm (going back to
            [NovThi06]_):

            Let `x` be an ordered set partition of `[m]`, and `y` an
            ordered set partition of `[n]`.
            Transform `x` into a list `u` of all the `m` elements of `[m]`
            by writing out each block of `x` (in increasing order) and
            putting bars between each two consecutive blocks; this is
            called a barred permutation.
            Do the same for `y`, but also shift each entry of the
            resulting barred permutation by `m`. Let `v` be the barred
            permutation of `[m+n] \setminus [m]` thus obtained.
            Now, shuffle the two barred permutations `u` and `v`
            (ignoring the bars) in all the `\binom{n+m}{n}` possible ways.
            For each shuffle obtained, place bars between some entries
            of the shuffle, according to the following rule:

            * If two consecutive entries of the shuffle both come from
              `u`, then place a bar between them if the corresponding
              entries of `u` had a bar between them.

            * If the first of two consecutive entries of the shuffle
              comes from `v` and the second from `u`, then place a bar
              between them.

            This results in a barred permutation of `[m+n]`.
            Transform it into an ordered set partition of `[m+n]`,
            by treating the bars as dividers separating consecutive
            blocks.

            The product `\Phi_x \Phi_y` is the sum of `\Phi_p` with
            `p` ranging over all ordered set partitions obtained this
            way.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).Phi()
                sage: x = OrderedSetPartition([[1],[2,3]])
                sage: y = OrderedSetPartition([[1,2]])
                sage: z = OrderedSetPartition([[1,2],[3]])
                sage: A.product_on_basis(x, y)
                Phi[{1}, {2, 3, 4, 5}] + Phi[{1}, {2, 4}, {3, 5}]
                 + Phi[{1}, {2, 4, 5}, {3}] + Phi[{1, 4}, {2, 3, 5}]
                 + Phi[{1, 4}, {2, 5}, {3}] + Phi[{1, 4, 5}, {2, 3}]
                 + Phi[{4}, {1}, {2, 3, 5}] + Phi[{4}, {1}, {2, 5}, {3}]
                 + Phi[{4}, {1, 5}, {2, 3}] + Phi[{4, 5}, {1}, {2, 3}]
                sage: A.product_on_basis(x, z)
                Phi[{1}, {2, 3, 4, 5}, {6}] + Phi[{1}, {2, 4}, {3, 5}, {6}]
                 + Phi[{1}, {2, 4, 5}, {3, 6}] + Phi[{1}, {2, 4, 5}, {6}, {3}]
                 + Phi[{1, 4}, {2, 3, 5}, {6}] + Phi[{1, 4}, {2, 5}, {3, 6}]
                 + Phi[{1, 4}, {2, 5}, {6}, {3}] + Phi[{1, 4, 5}, {2, 3, 6}]
                 + Phi[{1, 4, 5}, {2, 6}, {3}] + Phi[{1, 4, 5}, {6}, {2, 3}]
                 + Phi[{4}, {1}, {2, 3, 5}, {6}]
                 + Phi[{4}, {1}, {2, 5}, {3, 6}]
                 + Phi[{4}, {1}, {2, 5}, {6}, {3}]
                 + Phi[{4}, {1, 5}, {2, 3, 6}] + Phi[{4}, {1, 5}, {2, 6}, {3}]
                 + Phi[{4}, {1, 5}, {6}, {2, 3}] + Phi[{4, 5}, {1}, {2, 3, 6}]
                 + Phi[{4, 5}, {1}, {2, 6}, {3}] + Phi[{4, 5}, {1, 6}, {2, 3}]
                 + Phi[{4, 5}, {6}, {1}, {2, 3}]
                sage: A.product_on_basis(y, y)
                Phi[{1, 2, 3, 4}] + Phi[{1, 3}, {2, 4}] + Phi[{1, 3, 4}, {2}]
                 + Phi[{3}, {1, 2, 4}] + Phi[{3}, {1, 4}, {2}]
                 + Phi[{3, 4}, {1, 2}]

            TESTS::

                sage: one = OrderedSetPartition([])
                sage: all(A.product_on_basis(one, z) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
                sage: all(A.product_on_basis(z, one) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
                sage: M = algebras.WQSym(QQ).M()
                sage: x = A[[2, 4], [1, 3]]
                sage: y = A[[1, 3], [2]]
                sage: A(M(x) * M(y)) == x * y  # long time
                True
                sage: A(M(x) ** 2) == x**2 # long time
                True
                sage: A(M(y) ** 2) == y**2
                True
            """
            K = self.basis().keys()
            if not x:
                return self.monomial(y)
            if not y:
                return self.monomial(x)
            xlist = [(j, (k == 0))
                     for part in x
                     for (k, j) in enumerate(sorted(part))]
            # xlist is a list of the form
            # [(e_1, s_1), (e_2, s_2), ..., (e_n, s_n)],
            # where e_1, e_2, ..., e_n are the entries of the parts of
            # x in the order in which they appear in x (reading each
            # part from bottom to top), and where s_i = True if e_i is
            # the smallest element of its part and False otherwise.
            m = max(max(part) for part in x)  # The degree of x
            ylist = [(m + j, (k == 0))
                     for part in y
                     for (k, j) in enumerate(sorted(part))]
            # ylist is like xlist, but for y instead of x, and with
            # a shift by m.

            def digest(s):
                # Turn a shuffle of xlist with ylist into the appropriate
                # ordered set partition.
                s0 = [p[0] for p in s]
                s1 = [p[1] for p in s]
                N = len(s)
                bars = [False] * N
                for i in range(N - 1):
                    s0i = s0[i]
                    s0i1 = s0[i + 1]
                    if s0i <= m and s0i1 <= m:
                        bars[i + 1] = s1[i + 1]
                    elif s0i > m and s0i1 > m:
                        bars[i + 1] = s1[i + 1]
                    elif s0i > m and s0i1 <= m:
                        bars[i + 1] = True
                blocks = []
                block = []
                for i in range(N):
                    if bars[i]:
                        blocks.append(block)
                        block = [s0[i]]
                    else:
                        block.append(s0[i])
                blocks.append(block)
                return K(blocks)
            return self.sum_of_monomials(digest(s) for s in ShuffleProduct(xlist, ylist))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of ``self`` on the basis element
            indexed by the ordered set partition ``x``.

            The coproduct of the basis element `\Phi_x` indexed by
            an ordered set partition `x` of `[n]` can be computed by the
            following formula ([NovThi06]_):

            .. MATH::

                \Delta \Phi_x
                = \sum \Phi_y \otimes \Phi_z ,

            where the sum ranges over all pairs `(y, z)` of ordered set
            partitions `y` and `z` such that:

            * `y` and `z` are ordered set partitions of two complementary
              subsets of `[n]`;

            * `x` is obtained either by concatenating `y` and `z`, or by
              first concatenating `y` and `z` and then merging the two
              "middle blocks" (i.e., the last block of `y` and the first
              block of `z`); in the latter case, the maximum of the last
              block of `y` has to be smaller than the minimum of the first
              block of `z` (so that when merging these blocks, their
              entries don't need to be sorted).

            EXAMPLES::

                sage: Phi = algebras.WQSym(QQ).Phi()

                sage: Phi.coproduct(Phi.one())  # indirect doctest
                Phi[] # Phi[]
                sage: Phi.coproduct( Phi([[1]]) )  # indirect doctest
                Phi[] # Phi[{1}] + Phi[{1}] # Phi[]
                sage: Phi.coproduct( Phi([[1,2]]) )
                Phi[] # Phi[{1, 2}] + Phi[{1}] # Phi[{1}] + Phi[{1, 2}] # Phi[]
                sage: Phi.coproduct( Phi([[1], [2]]) )
                Phi[] # Phi[{1}, {2}] + Phi[{1}] # Phi[{1}] + Phi[{1}, {2}] # Phi[]
                sage: Phi[[1,2],[3],[4]].coproduct()
                Phi[] # Phi[{1, 2}, {3}, {4}] + Phi[{1}] # Phi[{1}, {2}, {3}]
                 + Phi[{1, 2}] # Phi[{1}, {2}] + Phi[{1, 2}, {3}] # Phi[{1}]
                 + Phi[{1, 2}, {3}, {4}] # Phi[]

            TESTS::

                sage: M = algebras.WQSym(QQ).M()
                sage: x = Phi[[2, 4], [6], [1, 3], [5, 7]]
                sage: MM = M.tensor(M); AA = Phi.tensor(Phi)
                sage: AA(M(x).coproduct()) == x.coproduct()
                True
            """
            if not len(x):
                return self.one().tensor(self.one())
            K = self.indices()

            def standardize(P):  # standardize an ordered set partition
                base = sorted(sum((list(part) for part in P), []))
                # base is the ground set of P, as a sorted list.
                d = {val: i + 1 for i, val in enumerate(base)}
                # d is the unique order isomorphism from base to
                # {1, 2, ..., |base|} (encoded as dict).
                return K([[d[x] for x in part] for part in P])
            deconcatenates = [(x[:i], x[i:]) for i in range(len(x) + 1)]
            for i in range(len(x)):
                xi = sorted(x[i])
                for j in range(1, len(xi)):
                    left = K(list(x[:i]) + [xi[:j]])
                    right = K([xi[j:]] + list(x[i + 1:]))
                    deconcatenates.append((left, right))
            T = self.tensor_square()
            return T.sum_of_monomials((standardize(left), standardize(right))
                                      for (left, right) in deconcatenates)

        class Element(WQSymBasis_abstract.Element):
            def algebraic_complement(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the algebraic complement involution.

                See
                :meth:`WQSymBases.ElementMethods.algebraic_complement`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`coalgebraic_complement`, :meth:`star_involution`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: Phi = WQSym.Phi()
                    sage: Phi[[1],[2,4],[3]].algebraic_complement()
                    -Phi[{3}, {2}, {4}, {1}] + Phi[{3}, {2, 4}, {1}] + Phi[{3}, {4}, {2}, {1}]
                    sage: Phi[[1],[2,3],[4]].algebraic_complement()
                    -Phi[{4}, {2}, {3}, {1}] + Phi[{4}, {2, 3}, {1}] + Phi[{4}, {3}, {2}, {1}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(Phi[A]).algebraic_complement()
                    ....:     == M(Phi[A].algebraic_complement())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.algebraic_complement doc
                # for the formula we're using here.
                BR = self.base_ring()
                one = BR.one()
                mine = -one
                Phi = self.parent()
                OSPs = Phi.basis().keys()
                from sage.data_structures.blas_dict import linear_combination

                def img(A):
                    # The image of the basis element Phi[A], written as a
                    # dictionary (of its coordinates in the Phi-basis).
                    Rs = [Rr.reversed() for Rr in A.strongly_finer()]
                    return {OSPs(P): (one if (len(R) % 2 == len(P) % 2)
                                      else mine)
                            for R in Rs for P in R.strongly_finer()}
                return Phi._from_dict(linear_combination((img(A), c) for (A, c) in self))

            def coalgebraic_complement(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the coalgebraic complement involution.

                See
                :meth:`WQSymBases.ElementMethods.coalgebraic_complement`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`algebraic_complement`, :meth:`star_involution`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: Phi = WQSym.Phi()
                    sage: Phi[[1],[2],[3,4]].coalgebraic_complement()
                    -Phi[{4}, {3}, {1}, {2}] + Phi[{4}, {3}, {1, 2}] + Phi[{4}, {3}, {2}, {1}]
                    sage: Phi[[2],[1,4],[3]].coalgebraic_complement()
                    -Phi[{3}, {1}, {4}, {2}] + Phi[{3}, {1, 4}, {2}] + Phi[{3}, {4}, {1}, {2}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(Phi[A]).coalgebraic_complement()
                    ....:     == M(Phi[A].coalgebraic_complement())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.coalgebraic_complement doc
                # for the formula we're using here.
                BR = self.base_ring()
                one = BR.one()
                mine = -one
                Phi = self.parent()
                OSPs = Phi.basis().keys()
                from sage.data_structures.blas_dict import linear_combination

                def img(A):
                    # The image of the basis element Phi[A], written as a
                    # dictionary (of its coordinates in the Phi-basis).
                    Rs = [Rr.complement() for Rr in A.strongly_finer()]
                    return {OSPs(P): (one if (len(R) % 2 == len(P) % 2)
                                      else mine)
                            for R in Rs for P in R.strongly_finer()}
                return Phi._from_dict(linear_combination((img(A), c) for (A, c) in self))

            def star_involution(self):
                r"""
                Return the image of the element ``self`` of `WQSym`
                under the star involution.

                See
                :meth:`WQSymBases.ElementMethods.star_involution`
                for a definition of the involution and for examples.

                .. SEEALSO::

                    :meth:`algebraic_complement`, :meth:`coalgebraic_complement`

                EXAMPLES::

                    sage: WQSym = algebras.WQSym(ZZ)
                    sage: Phi = WQSym.Phi()
                    sage: Phi[[1,2],[5,6],[3,4]].star_involution()
                    Phi[{3, 4}, {1, 2}, {5, 6}]
                    sage: Phi[[3], [1, 2], [4]].star_involution()
                    Phi[{1}, {3, 4}, {2}]

                TESTS::

                    sage: M = WQSym.M()
                    sage: all(M(Phi[A]).star_involution() == M(Phi[A].star_involution())
                    ....:     for A in OrderedSetPartitions(4))
                    True
                """
                # See the WQSymBases.ElementMethods.star_involution doc
                # for the formula we're using here.
                Phi = self.parent()
                OSPs = Phi.basis().keys()
                return Phi._from_dict({OSPs(A.complement().reversed()): c for (A, c) in self},
                                      remove_zeros=False)

    Phi = StronglyFiner


WQSymBasis_abstract.options = WordQuasiSymmetricFunctions.options


class WQSymBases(Category_realization_of_parent):
    r"""
    The category of bases of `WQSym`.
    """
    def __init__(self, base, graded):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- an instance of `WQSym`
        - ``graded`` -- boolean; if the basis is graded or filtered

        TESTS::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: bases = WQSymBases(WQSym, True)
            sage: WQSym.M() in bases
            True
        """
        self._graded = graded
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: WQSymBases(WQSym, True)
            Category of graded bases of Word Quasi-symmetric functions over Integer Ring
            sage: WQSymBases(WQSym, False)
            Category of filtered bases of Word Quasi-symmetric functions over Integer Ring
        """
        if self._graded:
            type_str = "graded"
        else:
            type_str = "filtered"
        return "Category of {} bases of {}".format(type_str, self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: bases = WQSymBases(WQSym, True)
            sage: bases.super_categories()
            [Category of realizations of Word Quasi-symmetric functions over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring
                 and Category of graded coalgebras over Integer Ring,
             Category of graded connected hopf algebras with basis over Integer Ring]

            sage: bases = WQSymBases(WQSym, False)
            sage: bases.super_categories()
            [Category of realizations of Word Quasi-symmetric functions over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring
                 and Category of graded coalgebras over Integer Ring,
             Join of Category of filtered connected hopf algebras with basis over Integer Ring
                 and Category of graded algebras over Integer Ring
                 and Category of graded coalgebras over Integer Ring]
        """
        R = self.base().base_ring()
        cat = HopfAlgebras(R).Graded().WithBasis()
        if self._graded:
            cat = cat.Graded()
        else:
            cat = cat.Filtered()
        return [self.base().Realizations(),
                HopfAlgebras(R).Graded().Realizations(),
                cat.Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `WQSym`.

            EXAMPLES::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: WQSym.M()
                Word Quasi-symmetric functions over Integer Ring in the Monomial basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, p):
            """
            Return the basis element indexed by ``p``.

            INPUT:

            - ``p`` -- an ordered set partition

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()
                sage: M[1, 2, 1]  # pass a word
                M[{1, 3}, {2}]
                sage: _ == M[[1, 2, 1]] == M[Word([1,2,1])]
                True
                sage: M[[1, 2, 3]]
                M[{1}, {2}, {3}]

                sage: M[[[1, 2, 3]]]  # pass an ordered set partition
                M[{1, 2, 3}]
                sage: _ == M[[1,2,3],] == M[OrderedSetPartition([[1,2,3]])]
                True
                sage: M[[1,3],[2]]
                M[{1, 3}, {2}]

            TESTS::

                sage: M[[]]
                M[]
                sage: M[1, 2, 1] == M[Word([2,3,2])] == M[Word('aca')]
                True
                sage: M[[[1,2]]] == M[1,1] == M[1/1,2/2] == M[2/1,2/1] == M['aa']
                True
                sage: M[1] == M[1,] == M[Word([1])] == M[OrderedSetPartition([[1]])] == M[[1],]
                True
            """
            if p in ZZ:
                p = [ZZ(p)]
            if all(s in ZZ for s in p):
                return self.monomial(self._indices.from_finite_word([ZZ(s) for s in p]))

            if all(isinstance(s, str) for s in p):
                return self.monomial(self._indices.from_finite_word(p))
            try:
                return self.monomial(self._indices(p))
            except TypeError:
                raise ValueError("cannot convert %s into an element of %s" % (p, self._indices))

        def is_field(self, proof=True):
            """
            Return whether ``self`` is a field.

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()
                sage: M.is_field()
                False
            """
            return False

        def is_commutative(self):
            """
            Return whether ``self`` is commutative.

            EXAMPLES::

                sage: M = algebras.WQSym(ZZ).M()
                sage: M.is_commutative()
                False
            """
            return self.base_ring().is_zero()

        def one_basis(self):
            """
            Return the index of the unit.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
                sage: A.one_basis()
                []
            """
            OSP = self.basis().keys()
            return OSP([])

        def degree_on_basis(self, t):
            """
            Return the degree of an ordered set partition in
            the algebra of word quasi-symmetric functions.

            This is the sum of the sizes of the blocks of the
            ordered set partition.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
                sage: u = OrderedSetPartition([[2,1],])
                sage: A.degree_on_basis(u)
                2
                sage: u = OrderedSetPartition([[2], [1]])
                sage: A.degree_on_basis(u)
                2
            """
            return sum(len(part) for part in t)

    class ElementMethods:
        def algebraic_complement(self):
            r"""
            Return the image of the element ``self`` of `WQSym`
            under the algebraic complement involution.

            If `u = (u_1, u_2, \ldots, u_n)` is a packed word
            that contains the letters `1, 2, \ldots, k` and no
            others, then the *complement* of `u` is defined to
            be the packed word
            `\overline{u} := (k+1 - u_1, k+1 - u_2, \ldots, k+1 - u_n)`.

            The algebraic complement involution is defined as the
            linear map `WQSym \to WQSym` that sends each basis
            element `\mathbf{M}_u` of the monomial basis of `WQSym`
            to the basis element `\mathbf{M}_{\overline{u}}`.
            This is a graded algebra automorphism and a coalgebra
            anti-automorphism of `WQSym`.
            Denoting by `\overline{f}` the image of an element
            `f \in WQSym` under the algebraic complement involution,
            it can be shown that every packed word `u` satisfies

            .. MATH::

                \overline{\mathbf{M}_u} = \mathbf{M}_{\overline{u}},
                \qquad \overline{X_u} = X_{\overline{u}},

            where standard notations for classical bases of `WQSym`
            are being used (that is, `\mathbf{M}` for the monomial
            basis, and `X` for the characteristic basis).

            This can be restated in terms of ordered set partitions:
            For any ordered set partition `R = (R_1, R_2, \ldots, R_k)`,
            let `R^r` denote the ordered set partition
            `(R_k, R_{k-1}, \ldots, R_1)`; this is known as
            the *reversal* of `R`. Then,

            .. MATH::

                \overline{\mathbf{M}_A} = \mathbf{M}_{A^r}, \qquad
                \overline{X_A} = X_{A^r}

            for any ordered set partition `A`.

            The formula describing algebraic complements on the Q basis
            (:class:`WordQuasiSymmetricFunctions.StronglyCoarser`)
            is more complicated, and requires some definitions.
            We define a partial order `\leq` on the set of all ordered
            set partitions as follows: `A \leq B` if and only if
            `A` is strongly finer than `B` (see
            :meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.is_strongly_finer`
            for a definition of this).
            The *length* `\ell(R)` of an ordered set partition `R` shall
            be defined as the number of parts of `R`.
            Use the notation `Q` for the Q basis.
            For any ordered set partition `A` of `[n]`, we have

            .. MATH::

                \overline{Q_A} = \sum_P c_{A, P} Q_P,

            where the sum is over all ordered set partitions `P` of
            `[n]`, and where the coefficient `c_{A, P}` is defined
            as follows:

            * If there exists an ordered set partition `R` satisfying
              `R \leq P` and `A \leq R^r`, then this `R` is unique,
              and `c_{A, P} = \left(-1\right)^{\ell(R) - \ell(P)}`.

            * If there exists no such `R`, then `c_{A, P} = 0`.

            The formula describing algebraic complements on the `\Phi`
            basis (:class:`WordQuasiSymmetricFunctions.StronglyFiner`)
            is identical to the above formula for the Q basis, except
            that the `\leq` sign has to be replaced by `\geq` in the
            definition of the coefficients `c_{A, P}`. In fact, both
            formulas are particular cases of a general formula for
            involutions:
            Assume that `V` is an (additive) abelian group, and that
            `I` is a poset. For each `i \in I`, let `M_i` be an element
            of `V`. Also, let `\omega` be an involution of the set `I`
            (not necessarily order-preserving or order-reversing),
            and let `\omega'` be an involutive group endomorphism of
            `V` such that each `i \in I` satisfies
            `\omega'(M_i) = M_{\omega(i)}`.
            For each `i \in I`, let `F_i = \sum_{j \geq i} M_j`,
            where we assume that the sum is finite.
            Then, each `i \in I` satisfies

            .. MATH::

                \omega'(F_i)
                = \sum_j \sum_{\substack{k \leq j; \\ \omega(k) \geq i}}
                  \mu(k, j) F_j,

            where `\mu` denotes the Möbius function. This formula becomes
            particularly useful when the `k` satisfying `k \leq j`
            and `\omega(k) \geq i` is unique (if it exists).
            In our situation, `V` is `WQSym`, and `I` is the set of
            ordered set partitions equipped either with the `\leq` partial
            order defined above or with its opposite order.
            The `M_i` is the `\mathbf{M}_A`, whereas the `F_i` is either
            `Q_i` or `\Phi_i`.

            If we denote the star involution
            (:meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution`)
            of the quasisymmetric functions by `f \mapsto f^{\ast}`,
            and if we let `\pi` be the canonical projection
            `WQSym \to QSym`, then each `f \in WQSym` satisfies
            `\pi(\overline{f}) = (\pi(f))^{\ast}`.

            .. SEEALSO::

                :meth:`coalgebraic_complement`, :meth:`star_involution`

            EXAMPLES:

            Recall that the index set for the bases of `WQSym` is
            given by ordered set partitions, not packed words.
            Translated into the language of ordered set partitions,
            the algebraic complement involution acts on the
            Monomial basis by reversing the ordered set partition.
            In other words, we have

            .. MATH::

                \overline{\mathbf{M}_{(P_1, P_2, \ldots, P_k)}}
                = \mathbf{M}_{(P_k, P_{k-1}, \ldots, P_1)}

            for any standard ordered set partition
            `(P_1, P_2, \ldots, P_k)`. Let us check this in practice::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: M = WQSym.M()
                sage: M[[1,3],[2]].algebraic_complement()
                M[{2}, {1, 3}]
                sage: M[[1,4],[2,5],[3,6]].algebraic_complement()
                M[{3, 6}, {2, 5}, {1, 4}]
                sage: (3*M[[1]] - 4*M[[]] + 5*M[[1],[2]]).algebraic_complement()
                -4*M[] + 3*M[{1}] + 5*M[{2}, {1}]
                sage: X = WQSym.X()
                sage: X[[1,3],[2]].algebraic_complement()
                X[{2}, {1, 3}]
                sage: C = WQSym.C()
                sage: C[[1,3],[2]].algebraic_complement()
                -C[{1, 2, 3}] - C[{1, 3}, {2}] + C[{2}, {1, 3}]
                sage: Q = WQSym.Q()
                sage: Q[[1,2],[5,6],[3,4]].algebraic_complement()
                Q[{3, 4}, {1, 2, 5, 6}] + Q[{3, 4}, {5, 6}, {1, 2}] - Q[{3, 4, 5, 6}, {1, 2}]
                sage: Phi = WQSym.Phi()
                sage: Phi[[2], [1,3]].algebraic_complement()
                -Phi[{1}, {3}, {2}] + Phi[{1, 3}, {2}] + Phi[{3}, {1}, {2}]

            The algebraic complement involution intertwines the antipode
            and the inverse of the antipode::

                sage: all( M(I).antipode().algebraic_complement().antipode()  # long time
                ....:      == M(I).algebraic_complement()
                ....:      for I in OrderedSetPartitions(4) )
                True

            Testing the `\pi(\overline{f}) = (\pi(f))^{\ast}` relation::

                sage: all( M[I].algebraic_complement().to_quasisymmetric_function()
                ....:      == M[I].to_quasisymmetric_function().star_involution()
                ....:      for I in OrderedSetPartitions(4) )
                True

            .. TODO::

                Check further commutative squares.
            """
            # Convert to the Monomial basis, there apply the algebraic
            # complement componentwise, then convert back.
            parent = self.parent()
            M = parent.realization_of().M()
            dct = {I.reversed(): coeff for (I, coeff) in M(self)}
            return parent(M._from_dict(dct, remove_zeros=False))

        def coalgebraic_complement(self):
            r"""
            Return the image of the element ``self`` of `WQSym`
            under the coalgebraic complement involution.

            If `u = (u_1, u_2, \ldots, u_n)` is a packed word,
            then the *reversal* of `u` is defined to be the
            packed word `(u_n, u_{n-1}, \ldots, u_1)`.
            This reversal is denoted by `u^r`.

            The coalgebraic complement involution is defined as the
            linear map `WQSym \to WQSym` that sends each basis
            element `\mathbf{M}_u` of the monomial basis of `WQSym`
            to the basis element `\mathbf{M}_{u^r}`.
            This is a graded coalgebra automorphism and an algebra
            anti-automorphism of `WQSym`.
            Denoting by `f^r` the image of an element `f \in WQSym`
            under the coalgebraic complement involution,
            it can be shown that every packed word `u` satisfies

            .. MATH::

                (\mathbf{M}_u)^r = \mathbf{M}_{u^r}, \qquad
                (X_u)^r = X_{u^r},

            where standard notations for classical bases of `WQSym`
            are being used (that is, `\mathbf{M}` for the monomial
            basis, and `X` for the characteristic basis).

            This can be restated in terms of ordered set partitions:
            For any ordered set partition `R` of `[n]`, let
            `\overline{R}` denote the complement of `R` (defined in
            :meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.complement`).
            Then,

            .. MATH::

                (\mathbf{M}_A)^r = \mathbf{M}_{\overline{A}}, \qquad
                (X_A)^r = X_{\overline{A}}

            for any ordered set partition `A`.

            Recall that `WQSym` is a subring of the ring of all
            bounded-degree noncommutative power series in countably many
            indeterminates. The latter ring has an obvious continuous
            algebra anti-endomorphism which sends each letter `x_i` to
            `x_i` (and thus sends each monomial
            `x_{i_1} x_{i_2} \cdots x_{i_n}` to
            `x_{i_n} x_{i_{n-1}} \cdots x_{i_1}`).
            This anti-endomorphism is actually an involution.
            The coalgebraic complement involution is simply the
            restriction of this involution to the subring `WQSym`.

            The formula describing coalgebraic complements on the Q basis
            (:class:`WordQuasiSymmetricFunctions.StronglyCoarser`)
            is more complicated, and requires some definitions.
            We define a partial order `\leq` on the set of all ordered
            set partitions as follows: `A \leq B` if and only if
            `A` is strongly finer than `B` (see
            :meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.is_strongly_finer`
            for a definition of this).
            The *length* `\ell(R)` of an ordered set partition `R` shall
            be defined as the number of parts of `R`.
            Use the notation `Q` for the Q basis.
            For any ordered set partition `A` of `[n]`, we have

            .. MATH::

                (Q_A)^r = \sum_P c_{A, P} Q_P ,

            where the sum is over all ordered set partitions `P` of
            `[n]`, and where the coefficient `c_{A, P}` is defined
            as follows:

            * If there exists an ordered set partition `R` satisfying
              `R \leq P` and `A \leq \overline{R}`, then this `R` is
              unique,
              and `c_{A, P} = \left(-1\right)^{\ell(R) - \ell(P)}`.

            * If there exists no such `R`, then `c_{A, P} = 0`.

            The formula describing coalgebraic complements on the `\Phi`
            basis (:class:`WordQuasiSymmetricFunctions.StronglyFiner`)
            is identical to the above formula for the Q basis, except
            that the `\leq` sign has to be replaced by `\geq` in the
            definition of the coefficients `c_{A, P}`. In fact, both
            formulas are particular cases of the general formula for
            involutions described in the documentation of
            :meth:`algebraic_complement`.

            If we let `\pi` be the canonical projection
            `WQSym \to QSym`, then each `f \in WQSym` satisfies
            `\pi(f^r) = \pi(f)`.

            .. SEEALSO::

                :meth:`algebraic_complement`, :meth:`star_involution`

            EXAMPLES:

            Recall that the index set for the bases of `WQSym` is
            given by ordered set partitions, not packed words.
            Translated into the language of ordered set partitions,
            the coalgebraic complement involution acts on the
            Monomial basis by complementing the ordered set partition.
            In other words, we have

            .. MATH::

                (\mathbf{M}_A)^r = \mathbf{M}_{\overline{A}}

            for any standard ordered set partition `P`.
            Let us check this in practice::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: M = WQSym.M()
                sage: M[[1,3],[2]].coalgebraic_complement()
                M[{1, 3}, {2}]
                sage: M[[1,2],[3]].coalgebraic_complement()
                M[{2, 3}, {1}]
                sage: M[[1], [4], [2,3]].coalgebraic_complement()
                M[{4}, {1}, {2, 3}]
                sage: M[[1,4],[2,5],[3,6]].coalgebraic_complement()
                M[{3, 6}, {2, 5}, {1, 4}]
                sage: (3*M[[1]] - 4*M[[]] + 5*M[[1],[2]]).coalgebraic_complement()
                -4*M[] + 3*M[{1}] + 5*M[{2}, {1}]
                sage: X = WQSym.X()
                sage: X[[1,3],[2]].coalgebraic_complement()
                X[{1, 3}, {2}]
                sage: C = WQSym.C()
                sage: C[[1,3],[2]].coalgebraic_complement()
                C[{1, 3}, {2}]
                sage: Q = WQSym.Q()
                sage: Q[[1,2],[5,6],[3,4]].coalgebraic_complement()
                Q[{1, 2, 5, 6}, {3, 4}] + Q[{5, 6}, {1, 2}, {3, 4}] - Q[{5, 6}, {1, 2, 3, 4}]
                sage: Phi = WQSym.Phi()
                sage: Phi[[2], [1,3]].coalgebraic_complement()
                -Phi[{2}, {1}, {3}] + Phi[{2}, {1, 3}] + Phi[{2}, {3}, {1}]

            The coalgebraic complement involution intertwines the antipode
            and the inverse of the antipode::

                sage: all( M(I).antipode().coalgebraic_complement().antipode()  # long time
                ....:      == M(I).coalgebraic_complement()
                ....:      for I in OrderedSetPartitions(4) )
                True

            Testing the `\pi(f^r) = \pi(f)` relation above::

                sage: all( M[I].coalgebraic_complement().to_quasisymmetric_function()
                ....:      == M[I].to_quasisymmetric_function()
                ....:      for I in OrderedSetPartitions(4) )
                True

            .. TODO::

                Check further commutative squares.
            """
            # Convert to the Monomial basis, there apply the coalgebraic
            # complement componentwise, then convert back.
            parent = self.parent()
            M = parent.realization_of().M()
            dct = {I.complement(): coeff for (I, coeff) in M(self)}
            return parent(M._from_dict(dct, remove_zeros=False))

        def star_involution(self):
            r"""
            Return the image of the element ``self`` of `WQSym`
            under the star involution.

            The star involution is the composition of the
            algebraic complement involution
            (:meth:`algebraic_complement`) with the coalgebraic
            complement involution (:meth:`coalgebraic_complement`).
            The composition can be performed in either order, as the
            involutions commute.

            The star involution is a graded Hopf algebra
            anti-automorphism of `WQSym`.
            Let `f^{\ast}` denote the image of an element
            `f \in WQSym` under the star involution.
            Let `\mathbf{M}`, `X`, `Q` and `\Phi` stand for the
            monomial, characteristic, Q and Phi bases of `WQSym`.
            For any ordered set partition `A` of `[n]`, we let
            `A^{\ast}` denote the complement
            (:meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.complement`)
            of the reversal
            (:meth:`~sage.combinat.set_partition_ordered.OrderedSetPartition.reversed`)
            of `A`. Then, for any ordered set partition `A` of `[n]`,
            we have

            .. MATH::

                (\mathbf{M}_A)^{\ast} = \mathbf{M}_{A^{\ast}}, \qquad
                (X_A)^{\ast} = X_{A^{\ast}}, \qquad
                (Q_A)^{\ast} = Q_{A^{\ast}}, \qquad
                (\Phi_A)^{\ast} = \Phi_{A^{\ast}} .

            The star involution
            (:meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution`)
            on the ring of noncommutative symmetric functions is a
            restriction of the star involution on `WQSym`.

            If we denote the star involution
            (:meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution`)
            of the quasisymmetric functions by `f \mapsto f^{\ast}`,
            and if we let `\pi` be the canonical projection
            `WQSym \to QSym`, then each `f \in WQSym` satisfies
            `\pi(f^{\ast}) = (\pi(f))^{\ast}`.

            .. TODO::

                More commutative diagrams?
                FQSym and FSym need their own star_involution
                methods defined first.

            .. SEEALSO::

                :meth:`algebraic_complement`, :meth:`coalgebraic_complement`

            EXAMPLES:

            Keep in mind that the default input method for basis keys
            of `WQSym` is by entering an ordered set partition, not a
            packed word. Let us check the basis formulas for the
            star involution::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: M = WQSym.M()
                sage: M[[1,3], [2,4,5]].star_involution()
                M[{1, 2, 4}, {3, 5}]
                sage: M[[1,3],[2]].star_involution()
                M[{2}, {1, 3}]
                sage: M[[1,4],[2,5],[3,6]].star_involution()
                M[{1, 4}, {2, 5}, {3, 6}]
                sage: (3*M[[1]] - 4*M[[]] + 5*M[[1],[2]]).star_involution()
                -4*M[] + 3*M[{1}] + 5*M[{1}, {2}]
                sage: X = WQSym.X()
                sage: X[[1,3],[2]].star_involution()
                X[{2}, {1, 3}]
                sage: C = WQSym.C()
                sage: C[[1,3],[2]].star_involution()
                -C[{1, 2, 3}] - C[{1, 3}, {2}] + C[{2}, {1, 3}]
                sage: Q = WQSym.Q()
                sage: Q[[1,3], [2,4,5]].star_involution()
                Q[{1, 2, 4}, {3, 5}]
                sage: Phi = WQSym.Phi()
                sage: Phi[[1,3], [2,4,5]].star_involution()
                Phi[{1, 2, 4}, {3, 5}]

            Testing the formulas for `(Q_A)^{\ast}` and `(\Phi_A)^{\ast}`::

                sage: all(Q[A].star_involution() == Q[A.complement().reversed()] for A in OrderedSetPartitions(4))
                True
                sage: all(Phi[A].star_involution() == Phi[A.complement().reversed()] for A in OrderedSetPartitions(4))
                True

            The star involution commutes with the antipode::

                sage: all( M(I).antipode().star_involution()  # long time
                ....:      == M(I).star_involution().antipode()
                ....:      for I in OrderedSetPartitions(4) )
                True

            Testing the `\pi(f^{\ast}) = (\pi(f))^{\ast}` relation::

                sage: all( M[I].star_involution().to_quasisymmetric_function()
                ....:      == M[I].to_quasisymmetric_function().star_involution()
                ....:      for I in OrderedSetPartitions(4) )
                True

            Testing the fact that the star involution on the
            noncommutative symmetric functions is a restriction of
            the star involution on `WQSym`::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: R = NCSF.R()
                sage: all(R[I].star_involution().to_fqsym().to_wqsym()
                ....:     == R[I].to_fqsym().to_wqsym().star_involution()
                ....:     for I in Compositions(4))
                True

            .. TODO::

                Check further commutative squares.
            """
            # Convert to the Monomial basis, there apply the algebraic
            # complement componentwise, then convert back.
            parent = self.parent()
            M = parent.realization_of().M()
            dct = {I.reversed().complement(): coeff for (I, coeff) in M(self)}
            return parent(M._from_dict(dct, remove_zeros=False))

        def to_quasisymmetric_function(self):
            r"""
            The projection of ``self`` to the ring `QSym` of
            quasisymmetric functions.

            There is a canonical projection `\pi : WQSym \to QSym`
            that sends every element `\mathbf{M}_P` of the monomial
            basis of `WQSym` to the monomial quasisymmetric function
            `M_c`, where `c` is the composition whose parts are the
            sizes of the blocks of `P`.
            This `\pi` is a ring homomorphism.

            OUTPUT:

            - an element of the quasisymmetric functions in the monomial basis

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()
                sage: M[[1,3],[2]].to_quasisymmetric_function()
                M[2, 1]
                sage: (M[[1,3],[2]] + 3*M[[2,3],[1]] - M[[1,2,3],]).to_quasisymmetric_function()
                4*M[2, 1] - M[3]
                sage: X, Y = M[[1,3],[2]], M[[1,2,3],]
                sage: X.to_quasisymmetric_function() * Y.to_quasisymmetric_function() == (X*Y).to_quasisymmetric_function()
                True

                sage: C = algebras.WQSym(QQ).C()
                sage: C[[2,3],[1,4]].to_quasisymmetric_function() == M(C[[2,3],[1,4]]).to_quasisymmetric_function()
                True

                sage: C2 = algebras.WQSym(GF(2)).C()
                sage: C2[[1,2],[3,4]].to_quasisymmetric_function()
                M[2, 2]
                sage: C2[[2,3],[1,4]].to_quasisymmetric_function()
                M[4]
            """
            from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
            M = QuasiSymmetricFunctions(self.parent().base_ring()).Monomial()
            MW = self.parent().realization_of().M()
            return M.sum_of_terms((i.to_composition(), coeff)
                                  for (i, coeff) in MW(self))
