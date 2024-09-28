r"""
Functions for computing the profiles of sub-Hopf algebras of the
mod p Steenrod algebra.

In particular, the main functions are :func:`profile_elt` and
:func:`enveloping_profile_elements`, which compute the profile
function for the smallest sub-Hopf algebra of the Steenrod algebra
containing the given elements.

See :func:`SteenrodAlgebra
<sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
:mod:`sage.algebras.steenrod.steenrod_algebra` for information about
profile functions.

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): initial version
- John Palmieri (2022): cleanup, modifications
"""

# ****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import zip_longest
from sage.rings.integer_ring import ZZ


def profile_elt(elt, char=2):
    """
    Return the smallest sub-Hopf algebra containing ``elt``.

    INPUT:

    - ``elt`` -- element of the Steenrod algebra (or a sub-Hopf algebra
      of it) or list(s) representing it
    - ``char`` -- (default: 2) the characteristic

    ``elt`` could also be a list (when ``char=2``) or a pair of lists
    (otherwise), in which case it is treated as corresponding to an
    element of the Steenrod algebra: ``(a, b, c) <-> Sq(a, b, c)`` or
    ``((a, b, c), (x, y, z)) <-> Q_a Q_b Q_c P(x, y, z)``.

    OUTPUT: the profile function corresponding to the smallest
    sub-Hopf algebra containing the element passed

    EXAMPLES::

        sage: from sage.modules.fp_graded.steenrod.profile import profile_elt
        sage: A2 = SteenrodAlgebra(2)
        sage: profile_elt(A2.Sq(2))
        (2, 1)
        sage: profile_elt(A2.Sq(4,8))
        (3, 4, 3, 2, 1)

        sage: B = SteenrodAlgebra(3)
        sage: b = B.an_element(); b
        2 Q_1 Q_3 P(2,1)
        sage: profile_elt(b, char=3)
        ((1, 1), (1, 2, 1, 2))
        sage: profile_elt(B.P(2,1), char=3)
        ((1, 1), ())
        sage: profile_elt(B.Q(2), char=3)
        ((0,), (1, 1, 2))
    """
    try:
        # Convert element of the Steenrod algebra to iterable(s).
        elt = elt.leading_support()
    except AttributeError:
        # Not in the Steenrod algebra so presumably already in the
        # right form.
        pass

    if char == 2:
        minprofile = [max(0, ZZ(n).exact_log(char) + 1) for n in elt]
        return find_min_profile(minprofile, char)

    # odd primes:
    alistQ, alistP = elt
    minprofileP = [max(0, ZZ(n).exact_log(char) + 1) for n in alistP]
    if not alistQ:
        minpQ = []
    else:
        minpQ = [1] * (max(alistQ) + 1)
        for j in alistQ:
            minpQ[j] = 2
    return find_min_profile((minprofileP, minpQ), char=char)


def enveloping_profile_elements(alist, char=2):
    r"""
    Return the profile function for the smallest sub-Hopf algebra
    containing the list of elements passed.

    INPUT:

    - ``alist`` -- list of Steenrod algebra elements
    - ``char`` -- (default: 2) the characteristic

    As with :func:`profile_elt`, the entries of ``alist`` could also
    be iterables or pairs of iterables.

    OUTPUT: the profile function for the minimum sub-algebra
    containing all the elements of ``alist``

    EXAMPLES::

        sage: from sage.modules.fp_graded.steenrod.profile import enveloping_profile_elements
        sage: enveloping_profile_elements([Sq(2),Sq(4)])
        (3, 2, 1)
        sage: enveloping_profile_elements([Sq(4)])
        (3, 2, 1)
        sage: enveloping_profile_elements([Sq(2,1,2),Sq(7)])
        (3, 2, 2, 1)

        sage: B = SteenrodAlgebra(3)
        sage: enveloping_profile_elements([B.P(2,1), B.P(0,0,3)], char=3)
        ((1, 1, 2, 1), ())
        sage: enveloping_profile_elements([B.P(3,1)], char=3)
        ((2, 1), ())
        sage: enveloping_profile_elements([B.P(2,1), B.P(0,0,3), B.Q(2)], char=3)
        ((1, 1, 2, 1), (1, 1, 2))
    """
    if char == 2:
        profiles = [profile_elt(x) for x in alist if x != 0]
        if not profiles:
            return (0,)
        # zip_longest doesn't do the right thing with lists of length
        # 1, but this case is simple:
        if len(profiles) == 1:
            return profiles[0]
        return find_min_profile(max(*a) for a in zip_longest(*profiles, fillvalue=0))

    # odd primes:
    profiles = [profile_elt(x, char) for x in alist if x != 0]
    if len(profiles) == 1:
        return profiles[0]
    profiles_P = [x[0] for x in profiles]
    profiles_Q = [x[1] for x in profiles]
    if not profiles_P and not profiles_Q:
        return ((0,), (0,))
    else:
        maxP = [max(*a) for a in zip_longest(*profiles_P, fillvalue=0)]
        maxQ = [max(*a) for a in zip_longest(*profiles_Q, fillvalue=0)]
    return find_min_profile([maxP, maxQ], char=char)


def find_min_profile(prof, char=2):
    r"""
    Return the smallest valid profile function containing a tuple of
    nonnegative integers.

    INPUT:

    - ``prof`` -- list or tuple of nonnegative integers
    - ``char`` -- (default: 2) the characteristic

    OUTPUT: a valid profile containing ``prof``

    A profile function `e` must satisfy `e(r) \geq \min( e(r-i) - i,
    e(i))` for all `0 < i < r`, and at odd primes, if `k(i+j) = 1`,
    then either `e(i) \leq j` or `k(j) = 1` for all `i \geq 1`, `j
    \geq 0`. We use these inequalities to generate the smallest
    profile function `e` satisfying `e(r) \geq prof(r)` for each `r`
    when `char=2`, and similarly at odd primes.

    EXAMPLES::

        sage: from sage.modules.fp_graded.steenrod.profile import find_min_profile
        sage: find_min_profile([1,2])
        (1, 2, 1)
        sage: find_min_profile([2,1])
        (2, 1)
        sage: find_min_profile([1,2,3])
        (1, 2, 3, 1, 1)
        sage: find_min_profile([4])
        (4, 3, 2, 1)

        sage: find_min_profile([[4], []], char=3)
        ((4, 3, 2, 1), ())
        sage: find_min_profile([[1], [2]], char=3)
        ((1,), (2, 2))
        sage: find_min_profile([[], [2,1,1,2]], char=3)
        ((0,), (2, 1, 1, 2))
    """
    if char == 2:
        if not prof:
            return (0,)
        # Add a zero to the front so that the relevant part of the
        # pseudo-profile new is indexed starting with 1.
        new = [0] + list(prof)
        # new is probably too short. Increase its length by using the
        # defining inequality for a profile function:
        #    e(s+t) \geq min(e(s)-t, e(t)).
        # So if min(e(s)-t, e(t)) > 0, then the length must be at
        # least s+t.
        pad = 0
        new += [0] * max(new)
        for t in range(len(new)):
            for s in range(len(new)):
                if min(new[s] - t, new[t]) > 0:
                    pad = max(pad, s + t)
        e = [0] * len(new)
        # Now compute the new profile e.
        for r in range(len(e)):
            m = max((min(e[r - i] - i, e[i]) for i in range(1, r)), default=0)
            e[r] = max(m, new[r])
        # Strip trailing zeroes.
        while e[-1] == 0:
            e = e[:-1]
        return tuple(e[1:])

    # odd primes:
    pP = list(prof[0])
    pQ = list(prof[1])
    P = find_min_profile(pP, char=2)
    if not pQ:
        return (P, tuple(pQ))
    # newQ: dictionary of the form {index: value} where value is
    # either 1 or 2.
    maxP = max(P)
    newQ = list(pQ) + [None] * maxP
    for j in range(len(pQ)):
        if newQ[j] == 2:
            for i in range(maxP):
                if P[i] > j:
                    newQ[i + 1 + j] = 2
    # Strip all of the None values
    # Do it from the back to minimize reshuffles and keep the index matching
    for i in range(len(newQ) - 1, len(pQ) - 1, -1):
        if newQ[i] is None:
            del newQ[i]
    return (P, tuple(newQ))
