r"""
Oriented matroid with covector axioms

This implements an oriented matroid using the covector axioms.

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
from __future__ import annotations


class CovectorOrientedMatroid(OrientedMatroid):
    r"""
    An oriented matroid implemented using covector axioms.

    According to  Definition 3.7.1 in [BLSWZ1999]_ a *covector* of an oriented
    matroid is any composition of cocircuits.

    This implements an oriented matroid using the covector axioms. For this
    let `\mathcal{L}` be a set of covectors and `E` a groundset. Then
    a pair `M = (E,\mathcal{L})` is an oriented matroid using the covectors
    axioms if (see Theorem 4.1.1 in [BLSWZ1999]_):

    - `0 \in \mathcal{L}`
    - `X \in \mathcal{L}` implies `-X \in \mathcal{L}`
    - For all `X,Y \in \mathcal{L}`, `X \circ Y \in \mathcal{L}`
    - For all `X,Y \in \mathcal{L}` and `e \in S(X,Y)` there exists a
      `Z \in \mathcal{L}` such that `Z(e) = 0` and
      `Z(f) = (X \circ Y)(f) = (Y \circ X)(f)` for all `f \notin S(X,Y)`.

    INPUT:

    - ``data`` -- a tuple containing :class:`SignedSubsetElement` elements or
      data that can be used to construct :class:`SignedSubsetElement` elements
    - ``groundset`` -- (default: ``None``) the groundset for the data; if not
      provided, we grab the data from the signed subsets

    EXAMPLES::

        sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        sage: M = OrientedMatroid([[1], [-1], [0]], groundset=['e'], key='covector')
        sage: M.groundset()
        ('e',)

        sage: C = [[1,1,1], [1,1,0], [1,1,-1], [1,0,-1], [1,-1,-1], [0,-1,-1],
        ....:      [-1,-1,-1], [0,1,1], [-1,1,1], [-1,0,1], [-1,-1,1], [-1,-1,0], [0,0,0]]
        sage: M = OrientedMatroid(C, key='covector'); M
        Covector oriented matroid of rank 2
        sage: M.groundset()
        (0, 1, 2)
        sage: gs = ['h1', 'h2', 'h3']
        sage: M = OrientedMatroid(C, key='covector', groundset=gs)
        sage: M.groundset()
        ('h1', 'h2', 'h3')

    .. SEEALSO::

        - :class:`~sage.oriented_matroids.oriented_matroid.OrientedMatroid`
    """
    def __init__(self, data, groundset=None):
        """
        Return a ``CovectorOrientedMatroid`` object.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1], [0]], groundset=['e'], key='covector'); M
            Covector oriented matroid of rank 1
            sage: TestSuite(M).run()
        """
        OrientedMatroid.__init__(self)

        # Set up our covectors
        covectors = []
        for d in data:
            # Ensure we're using the right type.
            covectors.append(SignedSubsetElement(self, data=d, groundset=groundset))
        # If our groundset is None, make sure the groundsets are the same for
        # all elements
        if groundset is None and len(covectors) > 0:
            groundset = covectors[0].groundset()
            for X in covectors:
                if X.groundset() != groundset:
                    raise ValueError("groundsets must be the same")

        self._covectors = covectors
        self._elements = covectors

        if groundset is None:
            self._groundset = groundset
        else:
            self._groundset = tuple(groundset)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1], [0]], groundset=['e'], key='covector')
            sage: M
            Covector oriented matroid of rank 1
        """
        try:
            rep = f"Covector oriented matroid of rank {self.rank()}"
        except ValueError:
            rep = "Covector oriented matroid"
        return rep

    def is_valid(self, certificate=False) -> bool | tuple[bool, dict]:
        """
        Return whether our covectors satisfy the covector axioms.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: M = OrientedMatroid([[1], [-1], [0]], groundset=['e'], key='covector')
            sage: M.is_valid()
            True

            sage: C2 = [[0,0], [1,1]]
            sage: M2 = OrientedMatroid(C2, key='covector')
            sage: M2.is_valid(certificate=True)
            (False,
             {'elt': +: 0,1
               -:
               0: ,
              'msg': 'every element needs an opposite'})

            sage: C3 = [[1,1], [-1,-1], [0,1], [1,0], [-1,0], [0,-1]]
            sage: M3 = OrientedMatroid(C3, key='covector')
            sage: M3.is_valid(certificate=True)
            (False,
             {'elt': (+: 1
               -:
               0: 0,
               +:
               -: 0,1
               0: ),
              'msg': 'composition must be in vectors'})

            sage: C4 = [[0,0], [1,1], [-1,-1], [1,-1], [-1,1]]
            sage: M4 = OrientedMatroid(C4, key='covector')
            sage: M4.is_valid(certificate=True)
            (False,
             {'elt': (+: 0,1
               -:
               0: ,
               +: 0
               -: 1
               0: ),
              'msg': 'weak elimination failed'})
        """
        covectors = self.covectors()

        zero_found = False
        for X in covectors:
            # Axiom 1: Make sure empty is not present
            if X.is_zero():
                zero_found = True
            # Axiom 2: Make sure negative exists
            if -X not in covectors:
                if certificate:
                    error_info = {
                        'msg': "every element needs an opposite",
                        'elt': X
                        }
                    return (False, error_info)
                return False
            for Y in covectors:
                # Axiom 3: Closed under composition
                if X.composition(Y) not in covectors:
                    if certificate:
                        error_info = {
                            'msg': "composition must be in vectors",
                            'elt': (X, Y)
                            }
                        return (False, error_info)
                    return False
                # Axiom 4: Weak elimination axiom
                E = X.separation_set(Y)
                ze = set(self.groundset()).difference(E)
                xy = X.composition(Y)
                for e in E:
                    found = False
                    for Z in covectors:
                        if found:
                            break
                        if Z(e) == 0:
                            found = True
                            for f in ze:
                                if Z(f) != xy(f):
                                    found = False
                    if not found:
                        if certificate:
                            error_info = {
                                'msg': "weak elimination failed",
                                'elt': (X, Y)
                                }
                            return (False, error_info)
                        return False

        if not zero_found:
            if certificate:
                error_info = {
                    'msg': "all zero covector is required",
                    'elt': None
                    }
                return (False, error_info)
            return False

        if certificate:
            return (True, {})
        return True

    def matroid(self):
        r"""
        Return the underlying matroid.

        Given a covector oriented matroid, the *underlying matroid* is the
        matroid whose collection of flats is given by the set of zeros of all
        signed vectors.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: C = [[1,1,1], [1,1,0], [1,1,-1], [1,0,-1], [1,-1,-1],
            ....:      [0,-1,-1], [-1,-1,-1], [0,1,1], [-1,1,1], [-1,0,1],
            ....:      [-1,-1,1], [-1,-1,0], [0,0,0]]
            sage: M = OrientedMatroid(C, key='covector')
            sage: M.matroid()
            Matroid of rank 2 on 3 elements with 5 flats
        """
        from sage.matroids.constructor import Matroid
        from sage.matroids.flats_matroid import FlatsMatroid
        flats = list(set([frozenset(X.zeros()) for X in self.elements()]))
        return FlatsMatroid(groundset=self.groundset(), flats=flats)