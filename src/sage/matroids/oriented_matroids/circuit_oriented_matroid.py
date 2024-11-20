r"""
Oriented matroid with circuit axioms

This implements an oriented matroid using the circuit axioms.

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
from sage.matroids.oriented_matroids.signed_subset_element import SignedSubsetElement

class CircuitOrientedMatroid(OrientedMatroid):
    r"""
    An oriented matroid implemented using circuit axioms.

    This implements an oriented matroid using the circuit axioms. For this
    let `\mathcal{C}` be a set of signed subsets and `E` a groundset. Then
    a pair `M = (E, \mathcal{C})` is an oriented matroid using the circuit
    axioms if (see Definition 3.2.1 in [BLSWZ1999]_):

    - `\emptyset \notin \mathcal{C}`
    - `\mathcal{C} = -\mathcal{C}`
    - For all `X, Y \in \mathcal{C}`, if the support of `X` is contained
      in the support of `Y` then `X = Y` or `X = -Y`
    - For all `X, Y \in \mathcal{C}`, `X \neq -Y`, and
      `e \in X^+ \cap Y^-` there exists a `Z \in \mathcal{C}` such that
      `Z^+ \subseteq (X^+ \cup Y^+) \backslash \left\{e\right\}` and
      `Z^- \subseteq (X^- \cup Y^-) \backslash \left\{e\right\}`.

    INPUT:

    - ``data`` -- a tuple containing :class:`SignedSubsetElement` elements or
      data that can be used to construct :class:`SignedSubsetElement` elements
    - ``groundset`` -- (default: ``None``) the groundset for the data; if not
      provided, we grab the data from the signed subsets

    EXAMPLES::

        sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        sage: M = OrientedMatroid([[1], [-1]], key='circuit'); M
        Circuit oriented matroid of rank 0
        sage: M.groundset()
        (0,)

        sage: C = [((1,4),(2,3)), ((2,3),(1,4))]
        sage: M = OrientedMatroid(C, key='circuit'); M
        Circuit oriented matroid of rank 3
        sage: M.groundset()
        (1, 2, 3, 4)
        sage: M.an_element() in M.elements()
        True
        sage: M.elements() == M.circuits()
        True

        sage: C5 = [((1,),(3,),(2,)), ((1,2),(3,),(4,)),
        ....:       ((3,),(1,),(2,)), ((3,),(1,2),(4,))]
        sage: OrientedMatroid(C5, key='circuit')
        Traceback (most recent call last):
        ...
        ValueError: groundsets must be the same

    .. SEEALSO::

        - :class:`~sage.oriented_matroids.oriented_matroid.OrientedMatroid`
    """
    def __init__(self, data, groundset=None):
        """
        Return a ``CircuitOrientedMatroid`` object.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1]], key='circuit'); M
            Circuit oriented matroid of rank 0
            sage: TestSuite(M).run()
        """
        OrientedMatroid.__init__(self)

        # Set up our circuits
        circuits = []
        if data:
            for d in data:
                # Convert to the appropriate element class
                circuits.append(SignedSubsetElement(self, data=d, groundset=groundset))

        # If our groundset is none, make sure the groundsets are the same for
        # all elements
        if groundset is None and len(circuits) > 0:
            if len(data[0]) < 3:
                groundset = []
                for X in circuits:
                    groundset = list(set(groundset + X.groundset()))
            else:
                groundset = circuits[0].groundset()
                for X in circuits:
                    if X.groundset() != groundset:
                        raise ValueError("groundsets must be the same")

        self._circuits = circuits
        self._elements = circuits

        if groundset is None:
            self._groundset = None
        else:
            self._groundset = tuple(groundset)

    def is_valid(self, certificate=False) -> bool | tuple[bool, dict]:
        """
        Return whether our circuits satisfy the circuit axioms.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: C = [((1,4),(2,3)), ((2,3),(1,4))]
            sage: M = OrientedMatroid(C, key='circuit')
            sage: M.is_valid(certificate=True)
            (True, {})
            sage: M.is_valid()
            True

            sage: C2 = [((1,4),(2,3)), ((1,3),(2,4)), ((2,3),(1,4))]
            sage: M2 = OrientedMatroid(C2, key='circuit')
            sage: M2.is_valid(certificate=True)
            (False,
             {'elt': (+: 1,4
               -: 2,3
               0: ,
               +: 1,3
               -: 2,4
               0: ),
              'msg': 'only same/opposites can have same support'})
            sage: M2.is_valid()
            False

            sage: C3 = [((),()), ((1,4),(2,3)), ((2,3),(1,4))]
            sage: M3 = OrientedMatroid(C3, key='circuit', groundset=[1,2,3,4])
            sage: M3.is_valid(certificate=True)
            (False,
             {'elt': +: 
               -: 
               0: 1,2,3,4,
              'msg': 'empty set not allowed'})

            sage: C4= [((1,),()), ((1,4),(2,3)), ((2,3),(1,4))]
            sage: M4 = OrientedMatroid(C4, key='circuit', groundset=[1,2,3,4])
            sage: M4.is_valid(certificate=True)
            (False,
             {'elt': +: 1
               -: 
               0: 2,3,4,
              'msg': 'every element needs an opposite'})
        """
        circuits = self.circuits()

        for X in circuits:
            # Axiom 1: Make sure empty is not present
            if X.is_zero():
                if certificate:
                    error_info = {
                        'msg': "empty set not allowed",
                        'elt': X
                        }
                    return (False, error_info)
                return False
            # Axiom 2: (symmetry) Make sure negative exists
            if -X not in circuits:
                if certificate:
                    error_info = {
                        'msg': "every element needs an opposite",
                        'elt': X
                        }
                    return (False, error_info)
                return False
            for Y in circuits:
                # Axiom 3: (incomparability) supports must not be contained
                if X.support().issubset(Y.support()):
                    if X != Y and X != -Y:
                        if certificate:
                            error_info = {
                                'msg': "only same/opposites can have same support",
                                'elt': (X, Y)
                                }
                            return (
                                False,
                                error_info
                            )
                        return False
                # Axiom 4: Weak elimination
                if X != -Y:
                    E = X.positives().intersection(Y.negatives())
                    for e in E:
                        p = X.positives().union(Y.positives())
                        p.discard(e)
                        n = X.negatives().union(Y.negatives())
                        n.discard(e)
                        found = False
                        for Z in circuits:
                            if found:
                                break
                            if Z.positives().issubset(p) and Z.negatives().issubset(n):
                                found = True
                        if not found:
                            if certificate:
                                error_info = {
                                    'msg': "weak elimination failed",
                                    'elt': (X, Y)
                                    }
                                return (False, error_info)
                            return False
        if certificate:
            return (True, {})
        return True

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: C = [((1,),(2,)), ((2,),(1,)), ((3,),(4,)), ((4,),(3,))]
            sage: OrientedMatroid(C, key='circuit', groundset=[1,2,3,4])
            Circuit oriented matroid of rank 2
        """
        try:
            rep = f"Circuit oriented matroid of rank {self.rank()}"
        except ValueError:
            rep = "Circuit oriented matroid"
        return rep

    def matroid(self):
        r"""
        Return the underlying matroid.

        Given an oriented matroid defined using circuits, the *underlying
        matroid* is the (circuit) matroid whose groundset is the groundset of
        the oriented matroid and the circuits are the set of supports of all
        the signed subsets.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: C = [((1,),(2,)), ((2,),(1,)), ((3,),(4,)), ((4,),(3,))]
            sage: M = OrientedMatroid(C, key='circuit', groundset=[1,2,3,4])
            sage: M.matroid()
            Matroid of rank 2 on 4 elements with 2 circuits
        """
        from sage.matroids.constructor import Matroid
        circs = list(set([frozenset(X.support()) for X in self.elements()]))
        return Matroid(groundset=self.groundset(), circuits=circs)