r"""
Minimum Weighted Feedback Vertex Set on Interval Graphs
=======================================================

This module implements a dynamic programming algorithm for computing a
Minimum Weighted Feedback Vertex Set (MWFVS) on interval graphs.

The algorithm is based on:

    Lu & Tang (1997),
    "A linear-time algorithm for the weighted feedback vertex problem
    on interval graphs", Information Processing Letters 61.

NOTE:
    This implementation uses an explicit dynamic programming formulation
    and runs in O(n^2) time and space. The original paper describes an
    O(n + m) algorithm using a more optimized structure.

DEFINITIONS:
    - ``CVS`` -- Cycle-Free Vertex Set: a subset `S` of `V` such that `G[S]` contains no cycle
    - ``FVS`` -- Feedback Vertex Set: a subset `F` of `V` such that `V \\ F` is a ``CVS``
    - ``MWCVS`` -- maximum-weight ``CVS``
    - ``MWFVS`` -- minimum-weight FVS = `V \\ MWCVS` 

EXAMPLES::

    sage: from sage.graphs.interval_mwfvs import mwfvs_interval
    sage: intervals = [(1,4),(2,6),(3,5),(7,9)]
    sage: weights = [3,1,4,2]
    sage: mwfvs_interval(intervals, weights)
    {1}

    sage: mwfvs_interval([], [])
    set()
"""

from typing import List, Tuple, Set, Dict, FrozenSet


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _set_weight(s: FrozenSet[int], weights: List[float]) -> float:
    """Return total weight of a set of vertices."""
    return sum(weights[v] for v in s)


def _max_weight_set(
    s1: FrozenSet[int],
    s2: FrozenSet[int],
    weights: List[float],
) -> FrozenSet[int]:
    """Return the set with maximum weight."""
    return s1 if _set_weight(s1, weights) >= _set_weight(s2, weights) else s2


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def mwcvs_interval(
    intervals: List[Tuple[float, float]],
    weights: List[float],
) -> Set[int]:
    """
    Compute a Maximum Weighted Cycle-Free Vertex Set (MWCVS)
    of a weighted interval graph.

    INPUT:

    - ``intervals`` -- list of endpoint pairs (left, right)
    - ``weights`` -- list of vertex weights

    OUTPUT:

    - set of vertex indices (0-based) forming MWCVS

    EXAMPLES::

        sage: intervals = [(1,4),(2,6),(3,5),(7,9)]
        sage: weights = [3,1,4,2]
        sage: mwcvs_interval(intervals, weights)
        {0, 2, 3}

        
        sage: intervals = [(1,5),(2,6),(3,7)]
        sage: weights = [5,1,4]
        sage: mwcvs_interval(intervals, weights)
        {0, 2}

        
        sage: intervals = [(1,4),(5,8),(9,12)]
        sage: weights = [2,3,4]
        sage: mwcvs_interval(intervals, weights)
        {0, 1, 2}
    """
    n = len(intervals)

    if n == 0:
        return set()

    if len(weights) != n:
        raise ValueError("Length of weights must match intervals")

    # ------------------------------------------------------------------
    # Enforce distinct endpoints
    # ------------------------------------------------------------------
    endpoints = [ep for iv in intervals for ep in iv]
    if len(endpoints) != len(set(endpoints)):
        raise ValueError("All interval endpoints must be distinct")

    # ------------------------------------------------------------------
    # Sort intervals by right endpoint
    # ------------------------------------------------------------------
    order = sorted(range(n), key=lambda i: intervals[i][1])

    max_ep = max(max(iv) for iv in intervals)

    # Augmented arrays index 0..n+1  (0 = left sentinel, n+1 = right sentinel)
    A = [0.0] * (n + 2) # left endpoints
    B = [0.0] * (n + 2) # right endpoints
    W = [0.0] * (n + 2) # weights

    # Sentinel 0
    A[0], B[0], W[0] = -1.0, 0.0, 1.0

    for k, idx in enumerate(order):
        A[k + 1] = intervals[idx][0]
        B[k + 1] = intervals[idx][1]
        W[k + 1] = weights[idx]

    # Sentinel n+1
    A[n + 1] = 2 * max_ep + 1
    B[n + 1] = 2 * max_ep + 2
    W[n + 1] = 1.0

    # ------------------------------------------------------------------
    # Compute PRED which is defined: 
    # PRED(i) = largest k with B[k] < A[i]
    # ------------------------------------------------------------------
    events = []
    for i in range(n + 2):
        events.append((A[i], 0, i))
        events.append((B[i], 1, i))
    events.sort()

    PRED = [0] * (n + 2)
    cur = 0

    for _, typ, i in events:
        if typ == 0:
            PRED[i] = cur
        else:
            cur = i

    # ------------------------------------------------------------------
    # Main DP 
    # ------------------------------------------------------------------
    # Each table is indexed by (i, j) representing subproblem 
    # I_ij: subproblem defined by intervals between i and j
    #
    # MA[i, j]: MWCVS in subproblem I_ij
    #
    # MB[i, j]: Same as MA[i,j], but MUST include both i and j
    #
    # MC[i, j]: Same as MA[i,j], but MUST include j

    MA: Dict[Tuple[int, int], FrozenSet[int]] = {}
    MB: Dict[Tuple[int, int], FrozenSet[int]] = {}
    MC: Dict[Tuple[int, int], FrozenSet[int]] = {}

    for i in range(1, n + 2):
        MA[0, i] = frozenset([0, i]) if W[i] > 0 else frozenset([0])
        MB[0, i] = frozenset([0, i])
        MC[0, i] = frozenset([0, i])


    for j in range(2, n + 2):
        for i in range(max(1, PRED[j]), j):

            # Case 1: No overlap
            # 
            # Independent (non-overlapping) intervals cannot form a cycle
            if B[i] < A[j]: 
                MB[i, j] = MC[i - 1, i] | frozenset([j])

                MA[i, j] = (
                    MA[i - 1, i] | frozenset([j])
                    if W[j] > 0 else MA[i - 1, i]
                )

                MC[i, j] = MA[i - 1, i] | frozenset([j])

            elif A[j] < A[i]:
                # Case 2: Overlapping intervals where j contains i
                
                k = PRED[i]

                MB[i, j] = MC[k, j] | frozenset([i])

                MA[i, j] = _max_weight_set(
                    _max_weight_set(MA[i - 1, i], MA[i - 1, j], W),
                    MB[i, j],
                    W,
                )

                MC[i, j] = _max_weight_set(MC[i - 1, j], MB[i, j], W)

            else:
                # Case 3: Partially Overlapping intervals 
                k = PRED[j]

                MB[i, j] = MC[k, i] | frozenset([j])

                # Only MB differs because it enforces inclusion of both i and j

                # Common recurrence for MA and MC because:
                # MA(i, j) considers all possibilities:
                # - exclude j
                # - exclude i
                # - include both i and j
                
                MA[i, j] = _max_weight_set(
                    _max_weight_set(MA[i - 1, i], MA[i - 1, j], W),
                    MB[i, j],
                    W,
                )

                MC[i, j] = _max_weight_set(MC[i - 1, j], MB[i, j], W)

    # ------------------------------------------------------------------
    # Extract solution
    # ------------------------------------------------------------------
    result = MA[n, n + 1] - frozenset([0, n + 1])
    return {order[i - 1] for i in result}


def mwfvs_interval(
    intervals: List[Tuple[float, float]],
    weights: List[float],
) -> Set[int]:
    """
    Compute a Minimum Weighted Feedback Vertex Set (MWFVS)
    of an interval graph.

    INPUT:

    - ``intervals`` -- list of endpoint pairs (left, right)
    - ``weights`` -- list of vertex weights

    OUTPUT:

    - set of vertex indices (0-based)

    EXAMPLES::

        sage: intervals = [(1,4),(2,6),(3,5),(7,9)]
        sage: weights = [3,1,4,2]
        sage: mwfvs_interval(intervals, weights)
        {1}

        
        sage: intervals = [(1,5),(2,6),(3,7)]
        sage: weights = [5,1,4]
        sage: mwfvs_interval(intervals, weights)
        {1}

        
        sage: intervals = [(1,4),(5,8),(9,12)]
        sage: weights = [2,3,4]
        sage: mwfvs_interval(intervals, weights)
        set()
    """
    n = len(intervals)
    cvs = mwcvs_interval(intervals, weights)
    return set(range(n)) - cvs
