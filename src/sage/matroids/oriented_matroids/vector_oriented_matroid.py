r"""
Oriented matroid with vector axioms

This implements an oriented matroid using the vector axioms.

AUTHORS:

- Aram Dermenjian (2019-07-12): Initial version
"""

# ****************************************************************************
#      Copyright (C) 2019   Aram Dermenjian <aram.dermenjian.math at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid


class VectorOrientedMatroid(OrientedMatroid):
    r"""
    An oriented matroid implemented using vector axioms.

    According to  Definition 3.7.1 in [BLSWZ1999]_ a *vector* of an oriented
    matroid is any composition of circuits.

    This implements an oriented matroid using the vectors axioms. For this
    let `\mathcal{V}` be a set of signed subsets and `E` a groundset. Then
    a pair `M = (E,\mathcal{V})` is an oriented matroid using the vector
    axioms if (see Theorem 3.7.5 and Corollary 3.7.9 in [BLSWZ1999]_):

    - `\emptyset \in \mathcal{V}`
    - `\mathcal{V} = -\mathcal{V}`
    - For all `X,Y \in \mathcal{V}`, `X \circ Y \in \mathcal{V}`
    - For all `X,Y \in \mathcal{V}` and `e \in X^+ \cap Y^-` there exists
      a `Z \in \mathcal{V}` such that
      `Z^+ \subseteq (X^+ \cup Y^+) \backslash \left\{e\right\}` and
      `Z^- \subseteq (X^- \cup Y^-) \backslash \left\{e\right\}` and
      `(X \backslash Y) \cup (Y \backslash X) \cup (X^+ \cap Y^+) \cup (X^- \cap Y^-) \subseteq Z`.

    INPUT:

    - ``data`` -- a tuple containing :class:`SignedSubsetElement` elements or
      data that can be used to construct :class:`SignedSubsetElement` elements
    - ``groundset`` -- (default: ``None``) the groundset for the data; if not
      provided, we grab the data from the signed subsets

    EXAMPLES::

        sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        sage: M = OrientedMatroid([[1],[-1],[0]], key='vector'); M
        Vector oriented matroid of rank 0
        sage: M.groundset()
        (0,)
        sage: M = OrientedMatroid([[1],[-1],[0]], key='vector', groundset=['e'])
        sage: M
        Vector oriented matroid of rank 0
        sage: M.groundset()
        ('e',)

    .. SEEALSO::

        - :class:`~sage.oriented_matroids.oriented_matroid.OrientedMatroid`
        - :class:`~sage.oriented_matroids.signed_subset_element.SignedSubsetElement`
    """
    def __init__(self, data, groundset=None, category=None):
        """
        Return a ``VectorOrientedMatroid`` object.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1],[-1],[0]], key='vector'); M
            Vector oriented matroid of rank 0
            sage: TestSuite(M).run()
        """
        OrientedMatroid.__init__(self, category=category)

        # Set up our vectors
        vectors = []
        for d in data:
            # Use the appropriate element
            vectors.append(self.element_class(
                self, data=d, groundset=groundset))

        # If our groundset is none, make sure the groundsets are the same for
        # all elements
        if groundset is None and bool(vectors):
            groundset = vectors[0].groundset()
            for X in vectors:
                if X.groundset() != groundset:
                    raise ValueError("groundsets must be the same")

        self._vectors = vectors
        self._elements = vectors

        if groundset is None:
            self._groundset = groundset
        else:
            self._groundset = tuple(groundset)

    def is_valid(self, certificate=False) -> bool | tuple[bool, str]:
        """
        Return whether our vectors satisfy the vector axioms.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: V2 = [[1, 1]]
            sage: M2 = OrientedMatroid(V2, key='vector')
            sage: M2.is_valid(certificate=True)
            (False, 'every element needs an opposite')

            sage: V3 = [[1,1], [-1,-1], [0,-1], [0,1], [-1,0], [1,0]]
            sage: M3 = OrientedMatroid(V3, key='vector')
            sage: M3.is_valid(certificate=True)
            (False, 'composition must be in vectors')

            sage: V4 = [[1,1], [-1,-1]]
            sage: M4 = OrientedMatroid(V4, key='vector')
            sage: M4.is_valid(certificate=True)
            (False, 'vector elimination failed')
        """
        vectors = self.vectors()

        zero_found = False
        for X in vectors:
            # Axiom 1: Make sure empty is not present
            if X.is_zero():
                zero_found = True
            # Axiom 2: Make sure negative exists
            if -X not in vectors:
                if certificate:
                    return (False, "every element needs an opposite")
                return False
            for Y in vectors:
                # Axiom 3: Closed under composition
                if X.composition(Y) not in vectors:
                    if certificate:
                        return (False, "composition must be in vectors")
                    return False
                # Axiom 4: Vector elimination
                E = X.positives().intersection(Y.negatives())

                ze1 = X.support().difference(Y.support())
                ze2 = Y.support().difference(X.support())
                ze3 = X.positives().intersection(Y.positives())
                ze4 = X.negatives().intersection(Y.negatives())
                ze = ze1.union(ze2).union(ze3).union(ze4)
                for e in E:
                    p = X.positives().union(Y.positives())
                    p.discard(e)
                    n = X.negatives().union(Y.negatives())
                    n.discard(e)
                    found = False
                    for Z in vectors:
                        if found:
                            break
                        if Z.positives().issubset(p) \
                                and Z.negatives().issubset(n) \
                                and ze.issubset(Z.support()):
                            found = True
                    if not found:
                        if certificate:
                            return (False, "vector elimination failed")
                        return False

        if not zero_found:
            if certificate:
                return (False, "empty set is required")
            return False

        if certificate:
            return (True, "")
        return True

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: V = [[1,1], [-1,-1], [0,0]]
            sage: M = OrientedMatroid(V, key='vector'); M
            Vector oriented matroid of rank 1
        """
        try:
            rep = f"Vector oriented matroid of rank {self.rank()}"
        except ValueError:
            rep = "Vector oriented matroid"
        return rep

    def matroid(self):
        r"""
        Return the underlying matroid.

        Given an oriented matroid defined using vector, the *underlying
        matroid* is the (circuit) matroid whose groundset is the groundset of
        the oriented matroid and the circuits are the set of supports of all
        the signed subsets.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: V = [[1,1],[-1,-1],[0,0]]
            sage: M = OrientedMatroid(V, key='vector')
            sage: M.matroid()
            Matroid of rank 1 on 2 elements with 1 circuits
        """
        circOM = self.convert_to('circuit')
        return circOM.matroid()

    def circuits(self):
        r"""
        Return the circuits.

        Given a vector oriented matroid, the set of circuits is the set
        `\min(V)` which denotes the set of inclusion-minimal (nonempty) signed
        subsets.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1], [0]], key='vector')
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
        from sage.combinat.posets.posets import Poset
        # remove 0
        vecs = [v for v in self.vectors() if not v.is_zero()]
        P = Poset([vecs, lambda x, y: x.is_restriction_of(y)])
        self._circuits = P.minimal_elements()
        return self._circuits