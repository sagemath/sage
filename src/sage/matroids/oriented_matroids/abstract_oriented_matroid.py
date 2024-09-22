r"""
Abstract class for oriented matroids

AUTHORS:

- Aram Dermenjian (2019-07-12): Initial version
- Elizabeth Flight (2023-08-01): Beta version
- Tudor Tanasa (2023-08-01): Beta version
"""

# *****************************************************************************
#  Copyright (C) 2019 Aram Dermenjian <aram.dermenjian.math at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
from sage.structure.global_options import GlobalOptions
from sage.matroids.oriented_matroids.signed_subset_element import SignedSubsetElement


class AbstractOrientedMatroid(UniqueRepresentation, Parent):
    r"""
    Abstract class for oriented matroids.

    .. SEEALSO::

        :class:`oriented_matroids.oriented_matroid.OrientedMatroid`
    """

    # List of all possible keys
    keys = ['circuit', 'covector', 'vector', 'real_hyperplane_arrangement']

    Element = SignedSubsetElement

    class options(GlobalOptions):
        r"""
        Options for oriented matroids.

        @OPTIONS@

        .. NOTE::

            Changing the ``convention`` for tableaux also changes the
            ``convention`` for partitions.
        """
        NAME = 'OrientedMatroids'
        display = dict(default="set",
                       description='Changes how signed subsets are displayed.',
                       values=dict(set='display as sets',
                                   vector='display as vectors',
                                   ),
                       )

    def __init__(self, category=None):
        if category is None:
            category = Sets()
        Parent.__init__(self, category=category)

    @abstract_method
    def is_valid(self) -> bool:
        r"""
        Return whether ``self`` satisfies the oriented matroid axioms.

        Given a set of objects, this method tests against
        a provided set of axioms for a given representation
        to ensure that we actually do have an oriented matroid.
        """
        pass

    def groundset(self):
        """
        Return the groundset of ``self``.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(2)
            sage: M = OrientedMatroid(A); M.groundset()
            (Hyperplane t0 - t1 + 0,)
        """
        return self._groundset

    def elements(self):
        """
        Return elements.

        The elements of an oriented matroid are the "defining" elements of
        the oriented matroid. For example, covectors are the elements of
        an oriented matroid defined using covectors.
        """
        return self._elements

    def circuits(self):
        """
        Return all circuits.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1]], key='circuit')
            sage: M.circuits()
            [+: 0
             -:
             0: ,
             +:
             -: 0
             0: ]
        """
        if hasattr(self, "_circuits"):
            return self._circuits
        raise NotImplementedError("circuits not implemented")

    def cocircuits(self):
        """
        Return all cocircuits.
        """
        if hasattr(self, "_cocircuits"):
            return self._cocircuits
        raise NotImplementedError("cocircuits not implemented")

    def vectors(self):
        """
        Return all vectors.
        """
        if hasattr(self, "_vectors"):
            return self._vectors
        raise NotImplementedError("vectors not implemented")
        pass

    def covectors(self):
        """
        Return all covectors.
        """
        if hasattr(self, "_covectors"):
            return self._covectors
        raise NotImplementedError("covectors not implemented")

    def convert_to(self, new_type=None):
        '''
        Return an oriented matroid of type specified.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1], [0]], key='vector')
            sage: M.convert_to('circuit')
            Circuit oriented matroid of rank 0
            sage: M.convert_to()
            Traceback (most recent call last):
            ...
            TypeError: must be given a type to convert to
        '''
        from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        if new_type is None:
            raise TypeError("must be given a type to convert to")
        elif new_type in AbstractOrientedMatroid.keys:
            if hasattr(self, new_type + 's'):
                els = getattr(self, new_type + 's')()
            else:
                raise NotImplementedError("no %ss() method found in oriented matroid" % (new_type,))
            return OrientedMatroid(els,
                                   key=new_type,
                                   groundset=self.groundset())
        else:
            raise NotImplementedError("type %s not implemented" % (new_type,))

    def dual(self):
        """
        Return the dual oriented matroid.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A)
            sage: m = M.matroid()
            sage: m.dual()
            Dual of 'Matroid of rank 2 on 3 elements with 5 flats'
            sage: m.dual() is m
            False
        """
        pass

    @cached_method
    def matroid(self):
        r"""
        Return the underlying matroid.
        """
        pass

    # @cached_method
    def rank(self):
        r"""
        Return the rank.

        The *rank* of an oriented matroid is the rank of its underlying
        matroid.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A); M.rank()
            2
            sage: A = hyperplane_arrangements.braid(4)
            sage: M = OrientedMatroid(A); M.rank()
            3
        """
        return self.matroid().rank()

    def an_element(self):
        """
        Return an arbitrary element.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: D = DiGraph({'v1': {'v2': 1, 'v3': 2, 'v4': 3},
            ....:              'v2': {'v3': 4, 'v4': 5},
            ....:              'v3': {'v4': 6}})
            sage: M = OrientedMatroid(D, key="circuit")
            sage: M.an_element() in M.circuits()
            True
        """
        from sage.misc.prandom import randint
        els = self.elements()
        i = randint(1, len(els))
        return els[i - 1]

    def face_poset(self, facade=False):
        r"""
        Return the (big) face poset.

        The *(big) face poset* is the poset on covectors such that `X \leq Y```self`` `
        the if and only if `S(X,Y) = \emptyset` and
        `\underline{Y} \subseteq \underline{X}`.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: C = [[1,1,1], [1,1,0], [1,1,-1], [1,0,-1], [1,-1,-1],
            ....:      [0,-1,-1], [-1,-1,-1], [0,1,1], [-1,1,1],
            ....:      [-1,0,1], [-1,-1,1], [-1,-1,0], [0,0,0]]
            sage: M = OrientedMatroid(C, key='covector')
            sage: M.face_poset()
            Finite meet-semilattice containing 13 elements
        """
        from sage.combinat.posets.lattices import MeetSemilattice
        els = self.covectors()
        rels = [
            (Y, X)
            for X in els
            for Y in els
            if Y.is_conformal_with(X) and Y.support().issubset(X.support())
        ]
        return MeetSemilattice((els, rels), cover_relations=False, facade=facade)

    def face_lattice(self, facade=False):
        r"""
        Return the (big) face lattice.

        The *(big) face lattice* is the (big) face poset with a top element
        added.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: C = [[1,1,1], [1,1,0], [1,1,-1], [1,0,-1], [1,-1,-1],
            ....:      [0,-1,-1], [-1,-1,-1], [0,1,1], [-1,1,1],
            ....:      [-1,0,1], [-1,-1,1], [-1,-1,0], [0,0,0]]
            sage: M = OrientedMatroid(C, key='covector')
            sage: M.face_lattice()
            Finite lattice containing 14 elements
        """
        from sage.combinat.posets.lattices import LatticePoset
        els = self.covectors()
        rels = [
            (Y, X)
            for X in els
            for Y in els
            if Y.is_conformal_with(X) and Y.support().issubset(X.support())
        ]

        # Add top element
        for i in els:
            rels.append((i, 1))
        els.append(1)
        return LatticePoset((els, rels), cover_relations=False, facade=facade)

    def topes(self):
        r"""
        Return the topes.

        A *tope* is the maximal covector in the face poset.
        """
        return self.face_poset(facade=True).maximal_elements()

    def tope_poset(self, base_tope, facade=False):
        r"""
        Return the tope poset.

        The tope poset is the poset `(\mathcal{T}, B)` where `\mathcal{T}`
        is the set of topes and `B` is a distinguished tope called the
        *base tope*. The order is given by inclusion of separation sets
        from the base tope: `X \leq Y` if and only if
        `S(B, X) \subseteq S(B, Y)`.
        """
        from sage.combinat.posets.posets import Poset
        els = self.topes()
        rels = [
            (X, Y)
            for X in els
            for Y in els
            if base_tope.separation_set(X).issubset(base_tope.separation_set(Y))
        ]

        return Poset((els, rels), cover_relations=False, facade=facade)

    def is_simplicial(self):
        r"""
        Return if the oriented matroid is simplicial.

        An oriented matroid is *simplicial* if every tope is simplicial.

        .. SEEALSO::

            :meth:`~sage.matroids.oriented_matroids.signed_subset_element.SignedSubsetElement.is_simplicial`
        """
        for t in self.topes():
            if not t.is_simplicial():
                return False
        return True

    def is_acyclic(self):
        r"""
        Return if oriented matroid is acyclic.

        A covector oriented matroid is *acyclic* if there exists a positive
        tope where a *positive tope* is defined as a tope with no
        negative part.
        """
        for t in self.topes():
            if len(t.negatives()) == 0:
                return True
        return False

    def deletion(self, change_set):
        r"""
        Return a covector oriented matroid of a deletion.

        Let `M = (E, \mathcal{L})` be an oriented matroid over a set `E`
        and a set of covectors `\mathcal{L}`. Given `A \subseteq E`, the
        *deletion* is the (covector) oriented matroid
        `M\backslash A = (E \backslash A, \mathcal{L} \backslash A)` where

        .. MATH::

            \mathcal{C} \backslash A = \left\{ X\mid_{E \backslash A} : X \in \mathcal{C}\right\}

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        """
        if change_set in self.groundset():
            change_set = set([change_set])
        else:
            change_set = set(change_set)

        from sage.matroids.oriented_matroids.oriented_matroid import deep_tupler
        groundset = set(self.groundset()).difference(change_set)
        groundset = deep_tupler(groundset)
        data = []
        for c in self.covectors():
            p = tuple(c.positives().difference(change_set))
            n = tuple(c.negatives().difference(change_set))
            z = tuple(c.zeros().difference(change_set))
            data.append((p, n, z))
        data = deep_tupler(data)

        from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        return OrientedMatroid(data, key='covector', groundset=groundset)

    def restriction(self, change_set):
        r"""
        Return a covector oriented matroid of a restriction.

        Given an oriented matroid `M = (E, \mathcal{L})` where `E` is a
        set and `\mathcal{L}` is the set of covectors. Given
        `A \subseteq E`, the *restriction* is the (covector) oriented
        matroid `M / A = (E \backslash A, \mathcal{C} / A)` where

        .. MATH::

            \mathcal{C} / A = \left\{ X\mid_{E \backslash A} : X \in \mathcal{C} \text{ and} A \subseteq X^0 \right\}

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A); M
            Hyperplane arrangement oriented matroid of rank 2
            sage: R = M.restriction(M.groundset()[1]); R
            Covector oriented matroid of rank 1
            sage: R.elements()
            [+:
             -:
             0: Hyperplane 0*t0 + t1 - t2 + 0,Hyperplane t0 + 0*t1 - t2 + 0,
             +: Hyperplane 0*t0 + t1 - t2 + 0,Hyperplane t0 + 0*t1 - t2 + 0
             -:
             0: ,
             +:
             -: Hyperplane 0*t0 + t1 - t2 + 0,Hyperplane t0 + 0*t1 - t2 + 0
             0: ]
        """
        if change_set in self.groundset():
            change_set = set([change_set])
        else:
            change_set = set(change_set)

        from sage.matroids.oriented_matroids.oriented_matroid import deep_tupler
        groundset = set(self.groundset()).difference(change_set)
        groundset = deep_tupler(groundset)
        data = []
        for c in self.covectors():
            p = tuple(c.positives().difference(change_set))
            n = tuple(c.negatives().difference(change_set))
            z = tuple(c.zeros().difference(change_set))
            if change_set.issubset(c.zeros()):
                data.append((p, n, z))
        data = deep_tupler(data)

        from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        return OrientedMatroid(data, key='covector', groundset=groundset)

    def loops(self):
        r"""
        Return the loops of an oriented matroid.

        A *loop* is an element `e \in E` such that there is a
        tope `T \in \mathcal{T}` with `T(e) = 0`. In particular
        if `T(e) = 0` for some `T`, then it is true for all
        `T \in \mathcal{T}`.
        """
        T = self.topes()[0]
        loops = []
        gs = self.groundset()
        for i, j in enumerate(T):
            if T(j) == 0:
                loops.append(gs[i])
        return loops

    def are_parallel(self, e, f):
        r"""
        Return whether two elements in groundset are parallel.

        Two elements in the groundset `e, f \in E` are parallel if they
        are not loops and for all `X \in \mathcal{C}`, `X(e) = 0`
        implies `X(f) = 0`. See Lemma 4.1.10 [BLSWZ1999]_ .
        """
        gs = set(self.groundset()).difference(set(self.loops()))
        if e not in gs or f not in gs:
            raise ValueError(
                "Elements must be in groundset and must not be loops")
        for i in self.elements():
            if i(e) == 0 and i(f) != 0:
                return False
        return True

    def is_simple(self):
        r"""
        Return if the oriented matroid is simple.

        An oriented matroid is *simple* if there are no loops
        and no parallel elements.
        """
        from sage.combinat.subset import Subsets
        if len(self.loops()) > 0:
            return False
        for i in Subsets(self.groundset(), 2):
            if self.are_parallel(i[0], i[1]):
                return False
        return True

    def _element_constructor_(self, x):
        r"""
        Determine if ``x`` may be viewed as belonging to ``self``.
        """
        try:
            if x in self.elements():
                return x
            return False
        except ValueError:
            return False
