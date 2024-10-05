r"""
Oriented matroids construction

Theory
======

Oriented matroids are a generalization of directed graphs, central hyperplane
arrangements, vector arrangements, and other mathematical objects. There are
many cryptomorphic definitions of oriented matroids. Precise definitions
for each type can be found in that class's directory.

See :wikipedia:`Oriented_matroid` for more details.


Built-in oriented matroids
==========================

The types of oriented matroids currently implemented into sage are:

- :class:`Circuit Oriented Matroids<sage.matroids.oriented_matroids.circuit_oriented_matroid.CircuitOrientedMatroid>`
- :class:`Covector Oriented Matroids<sage.matroids.oriented_matroids.covector_oriented_matroid.CovectorOrientedMatroid>`
- :class:`Vector Oriented Matroids<sage.matroids.oriented_matroids.vector_oriented_matroid.VectorOrientedMatroid>`
- :class:`(Real) Hyperplane Arrangement Oriented Matroids<sage.matroids.oriented_matroids.real_hyperplane_arrangement_oriented_matroid.RealHyperplaneArrangementOrientedMatroid>`

Constructing oriented matroids
==============================

To define your own oriented matroid, you can call the function
`OrientedMatroid(data, key)` where `data` is the data of the oriented matroid
and the `key` is the type of oriented matroid you are constructing. In case you
pass an object as the data (such as a hyperplane arrangement, digraph,
etc.) the code will try and create an oriented matroid for you.


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

from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
from sage.structure.global_options import GlobalOptions
from sage.matroids.oriented_matroids.signed_subset_element import SignedSubsetElement

from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangementElement
from sage.geometry.triangulation.point_configuration import PointConfiguration
from sage.graphs.digraph import DiGraph
from sage.structure.element import Matrix
import copy


class OrientedMatroid(Parent, metaclass=ClasscallMetaclass):
    r"""
    Construct an oriented matroid.

    The implementation of the oriented matroid determines which axiom set will
    be used.

    INPUT:

    - ``data`` -- (default: ``None``) the data that defines the oriented
      matroid. It can be one of the following:

      + Objects

          + Hyperplane Arrangement
          + Point Configuration
          + Digraph
          + Matrix (not yet implemented)

      + A list or tuple of

          + :class:`SignedSubsetElement`
          + A tuple with positive, negative, and zero sets.

    - ``groundset`` -- (default: ``None``) the groundset of the oriented
      matroid

    - ``key`` -- (default: ``None``) the representation of the oriented
      matroid. It can be one of the following:

      + ``'covector'`` - uses covector axioms with covectors
      + ``'vector'`` - uses vector axioms with signed subsets
      + ``'circuit'`` - uses circuit axioms with signed subsets
      + ``None`` - try and guess key.

    EXAMPLES::

        sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        sage: A = hyperplane_arrangements.braid(3)
        sage: M = OrientedMatroid(A); M
        Hyperplane arrangement oriented matroid of rank 2
        sage: M.groundset()
        (Hyperplane 0*t0 + t1 - t2 + 0,
         Hyperplane t0 - t1 + 0*t2 + 0,
         Hyperplane t0 + 0*t1 - t2 + 0)
        sage: OrientedMatroid.options.display='vector'
        sage: M.elements()
        [(0,0,0),
         (0,1,1),
         (0,-1,-1),
         (1,0,1),
         (1,1,1),
         (1,-1,0),
         (1,-1,1),
         (1,-1,-1),
         (-1,0,-1),
         (-1,1,0),
         (-1,1,1),
         (-1,1,-1),
         (-1,-1,-1)]

        sage: D = DiGraph({'v1': {'v2': 1, 'v3': 2,'v4': 3},
        ....:              'v2': {'v3': 4, 'v4': 5},
        ....:              'v3': {'v4': 6}})
        sage: M = OrientedMatroid(D, key="circuit"); M
        Circuit oriented matroid of rank 3
        sage: len(M.circuits())
        14

        sage: PC = PointConfiguration([[1, 0, 0], [0, 1, 0], [0, 0, 1],
        ....:                          [1/2, 1/2, 0], [0, 1/2, 1/2], [1/3, 1/3, 1/3]])
        sage: M = OrientedMatroid(PC); M
        Circuit oriented matroid of rank 3
        sage: M.matroid()
        Matroid of rank 3 on 6 elements with 7 circuits

        sage: OrientedMatroid([[0]], key='covector')
        Covector oriented matroid of rank 0
        sage: M = OrientedMatroid([[0]], key='circuit')
        sage: M.is_valid(certificate=True)
        (False, 'empty set not allowed')
        sage: OrientedMatroid.options.display='set'

    OUTPUT:

    An oriented matroid whose axioms are determined by the type.

    .. TODO::

        - Currently, chirotopes are not implemented.

    REFERENCES:

    For more information see [BLSWZ1999]_.
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

    def __classcall_private__(self, data=None, groundset=None, key=None, **kwds):
        # Instantiate oriented matroid
        OM = None

        # If we have a hyperplane arrangement we need to force the key to be an
        # arrangement,
        if isinstance(data, HyperplaneArrangementElement):
            if key != 'covector' and key is not None:
                raise ValueError(
                    'hyperplane arrangements are currently only implemented using covector axioms')
            key = 'real_hyperplane_arrangement'
        elif isinstance(data, PointConfiguration):
            if key != 'circuit' and key is not None:
                raise ValueError(
                    'point configurations are currently only implemented using circuit axioms')
            key = 'circuit'
            # PC circuits are given as (+, 0, -); and only half are given
            ci = [(C[0], C[2], C[1]) for C in data.circuits()]
            ci += [(C[2], C[0], C[1]) for C in data.circuits()]
            data = ci
        elif isinstance(data, DiGraph):
            if key != 'circuit' and key is not None:
                raise ValueError(
                    'digraphs are currently only implemented using circuit axioms')
            key = 'circuit'

            # we need to add negative edges in order to do all simple cycles
            digraph = copy.copy(data)
            edges = copy.copy(list(digraph.edges(sort=True)))
            groundset = []
            if len(edges) != len(set(edges)):
                raise ValueError('edge labels need to be unique')
            if None in digraph.edge_labels():
                raise ValueError('edge labels must be set for all edges')

            # Add minus edges to properly get cycles
            for e in edges:
                digraph.add_edge(e[1], e[0], "NEG_" + str(e[2]))
                groundset.append(str(e[2]))
            # Each cycle defines a circuit
            data = []
            for c in digraph.all_cycles_iterator(simple=True):
                p = set([])
                n = set([])
                for e in range(len(c) - 1):
                    e = str(digraph.edge_label(c[e], c[e + 1]))
                    if e.startswith('NEG_'):
                        n.add(e.strip('NEG_'))
                    else:
                        p.add(e)
                # If an edge exists in both sets, then this is a false cycle.
                # This implies we have ee^-1 which is why it's false.
                # So we only add the true ones.
                if len(p.intersection(n)) == 0:
                    data.append([p, n])
        elif isinstance(data, Matrix):
            if key != 'chirotope' and key is not None:
                raise ValueError(
                    'matrices are currently only implemented using chirotope axioms')
            key = 'chirotope'

        if key not in OrientedMatroid.keys:
            raise ValueError("invalid type key")

        # In the following cases, deep_tupler is used since we are using
        # UniqueRepresentation Which doesn't allow us to have non-hashable things.
        if key == "covector":
            from sage.matroids.oriented_matroids.covector_oriented_matroid \
                import CovectorOrientedMatroid
            data = deep_tupler(data)
            if groundset is not None:
                groundset = deep_tupler(groundset)
            OM = CovectorOrientedMatroid(data, groundset=groundset)
        elif key == "circuit":
            from sage.matroids.oriented_matroids.circuit_oriented_matroid \
                import CircuitOrientedMatroid
            data = deep_tupler(data)
            if groundset is not None:
                groundset = deep_tupler(groundset)
            OM = CircuitOrientedMatroid(data, groundset=groundset)
        elif key == "vector":
            from sage.matroids.oriented_matroids.vector_oriented_matroid \
                import VectorOrientedMatroid
            data = deep_tupler(data)
            if groundset is not None:
                groundset = deep_tupler(groundset)
            OM = VectorOrientedMatroid(data, groundset=groundset)
        elif key == "real_hyperplane_arrangement":
            from sage.matroids.oriented_matroids.real_hyperplane_arrangement_oriented_matroid \
                import RealHyperplaneArrangementOrientedMatroid
            A = copy.copy(data)
            if groundset is None:
                groundset = deep_tupler(A.hyperplanes())
            else:
                groundset = deep_tupler(groundset)
            OM = RealHyperplaneArrangementOrientedMatroid(A, groundset=groundset)

        if OM is None:
            raise NotImplementedError(
                f"oriented matroid of type {key} is not implemented")

        return OM


    def __init__(self, category=None):
        if category is None:
            category = Sets()
        Parent.__init__(self, category=category)

    @abstract_method
    def is_valid(self, certificate=False) -> bool | tuple[bool, str]:
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
        Return all elements.

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
        elif new_type in OrientedMatroid.keys:
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
            raise ValueError("elements must be in groundset and must not be loops")
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




def deep_tupler(obj):
    r"""
    Change a (nested) list or set into a (nested) tuple to be hashable.
    """
    if isinstance(obj, list) or isinstance(obj, set):
        return tuple([deep_tupler(i) for i in obj])
    return obj
