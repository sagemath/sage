r"""
Parameters for hypergeometric functions

AUTHORS:

- Xavier Caruso, Florian Fürnsinn (2025-10): initial version
"""

# ***************************************************************************
#    Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#                       Florian Fürnsinn <florian.fuernsinn@univie.ac.at>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.misc.cachefunc import cached_method

from sage.functions.other import ceil
from sage.arith.functions import lcm

from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.semirings.tropical_semiring import TropicalSemiring

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix


# HypergeometricParameters of hypergeometric functions
########################################

class HypergeometricParameters():
    r"""
    Class for parameters of hypergeometric functions.
    """
    def __init__(self, top, bottom, add_one=True):
        r"""
        Initialize this set of parameters.

        INPUT:

        - ``top`` -- list of top parameters

        - ``bottom`` -- list of bottom parameters

        - ``add_one`` -- boolean (default: ``True``),
          if ``True``, add an additional one to the bottom
          parameters.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: type(pa)
            <class 'sage.functions.hypergeometric_parameters.HypergeometricParameters'>

        By default, parameters are sorted, duplicates are removed and
        a trailing `1` is added to the bottom parameters::

            sage: HypergeometricParameters([1/2, 1/3, 2/3], [2/3])
            ((1/3, 1/2), (1,))

        We can avoid adding the trailing `1` by passing ``add_one=False``::

            sage: HypergeometricParameters([1/2, 1/3, 2/3], [2/3], add_one=False)
            ((1/3, 1/2, 1), (1,))
        """
        try:
            top = sorted([QQ(a) for a in top if a is not None])
            bottom = sorted([QQ(b) for b in bottom if b is not None])
        except TypeError:
            raise NotImplementedError("parameters must be rational numbers")
        i = j = 0
        while i < len(top) and j < len(bottom):
            if top[i] == bottom[j]:
                del top[i]
                del bottom[j]
            elif top[i] > bottom[j]:
                j += 1
            else:
                i += 1
        if add_one:
            bottom.append(QQ(1))
        else:
            try:
                i = bottom.index(QQ(1))
                bottom.append(QQ(1))
                del bottom[i]
            except ValueError:
                bottom.append(QQ(1))
                top.append(QQ(1))
                top.sort()
        self.top = tuple(top)
        self.bottom = tuple(bottom)
        if len(top) == 0 and len(bottom) == 0:
            self.d = 1
            self.bound = 1
        else:
            self.d = lcm([ a.denominator() for a in top ]
                       + [ b.denominator() for b in bottom ])
            self.bound = 2 * self.d * max(top + bottom) + 1

    def __repr__(self):
        r"""
        Return a string representation of these parameters.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa  # indirect doctest
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
        """
        return "(%s, %s)" % (self.top, self.bottom)

    def __hash__(self):
        return hash((self.top, self.bottom))

    def __eq__(self, other):
        return (isinstance(other, HypergeometricParameters)
            and self.top == other.top and self.bottom == other.bottom)

    def is_balanced(self):
        r"""
        Return ``True`` if there are as many top parameters as bottom
        parameters; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.is_balanced()
            True

        ::

            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5], add_one=False)
            sage: pa
            ((1/4, 1/3, 1/2, 1), (2/5, 3/5, 1))
            sage: pa.is_balanced()
            False
        """
        return len(self.top) == len(self.bottom)

    @cached_method
    def christol_sorting(self, c=1):
        r"""
        Return a sorted list of triples, where each triple is associated to one
        of the parameters a, and consists of the decimal part of d*c*a (where
        integers are assigned 1 instead of 0), the negative value of a, and a
        sign (plus or minus 1), where top parameters are assigned -1 and bottom
        parameters +1. Sorting the list lexecographically according to the
        first two entries of the tuples sorts the corresponing parameters
        according to the total ordering (defined on p.6 in [Chr1986]_).

        INPUT:

        - ``c`` -- an integer (default: ``1``)

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.christol_sorting(7)
            [(12, -3/5, 1),
             (20, -1/3, -1),
             (30, -1/2, -1),
             (45, -1/4, -1),
             (48, -2/5, 1),
             (60, -1, 1)]
        """
        d = self.d
        A = [(d - (-d*c*a) % d, -a, -1) for a in self.top]
        B = [(d - (-d*c*b) % d, -b, 1) for b in self.bottom]
        return sorted(A + B)

    def parenthesis_criterion(self, c):
        r"""
        Return ``True`` if in each prefix of the list
        ``self.christol_sorting(c)`` there are at least as many triples with
        third entry -1 as triples with third entry +1. Return ``False``
        otherwise.

        INPUT:

        - ``c`` -- an integer

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.christol_sorting(7)
            [(12, -3/5, 1),
             (20, -1/3, -1),
             (30, -1/2, -1),
             (45, -1/4, -1),
             (48, -2/5, 1),
             (60, -1, 1)]
            sage: pa.parenthesis_criterion(7)
            False
            sage: pa.christol_sorting(1)
            [(15, -1/4, -1),
             (20, -1/3, -1),
             (24, -2/5, 1),
             (30, -1/2, -1),
             (36, -3/5, 1),
             (60, -1, 1)]
            sage: pa.parenthesis_criterion(1)
            True
        """
        parenthesis = 0
        for _, _, paren in self.christol_sorting(c):
            parenthesis += paren
            if parenthesis > 0:
                return False
        return parenthesis <= 0

    def interlacing_criterion(self, c):
        r"""
        Return ``True`` if the sorted lists of the decimal parts (where integers
        are assigned 1 instead of 0) of c*a and c*b for a in the top parameters
        and b in the bottom parameters interlace, i.e., the entries in the sorted
        union of the two lists alternate between entries from the first and from
        the second list. Used to determine algebraicity of the hypergeometric
        function with these parameters with the Beukers-Heckman criterion.

        INPUT:

        - ``c`` -- an integer

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/3, 2/3], [1/2])
            sage: pa
            ((1/3, 2/3), (1/2, 1))
            sage: pa.interlacing_criterion(1)
            True
            sage: pa.interlacing_criterion(5)
            True

        ::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/8, 3/8, 5/8], [1/4, 1/2])
            sage: pa
            ((1/8, 3/8, 5/8), (1/4, 1/2, 1))
            sage: pa.interlacing_criterion(1)
            True
            sage: pa.interlacing_criterion(3)
            False
        """
        previous_paren = 1
        for _, _, paren in self.christol_sorting(c):
            if paren == previous_paren:
                return False
            previous_paren = paren
        return True

    def q_christol_sorting(self, q):
        r"""
        Return a sorted list of pairs, one associated to each top parameter a,
        and one associated to each bottom parameter b where the pair is either
        (1/2 + (-a) % q, -1) or (1 + (-b) % q, 1).

        INPUT:

        - ``q`` -- an integer

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.q_christol_sorting(7)
            [(2, 1), (2.5, -1), (3.5, -1), (5.5, -1), (6, 1), (7, 1)]
        """
        A = [(1/2 + (-a) % q, -1) for a in self.top]
        B = [(1 + (-b) % q, 1) for b in self.bottom]
        return sorted(A + B)

    def q_parenthesis(self, q):
        r"""
        Return maximal value of the sum of all the second entries of the pairs
        in a prefix of ``self.q_christol_sorting(q)`` and the first entry of
        the last pair in the prefix of smallest length where this value is
        attained.

        INPUT:

        - ``q`` -- an integer

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.q_christol_sorting(7)
            [(2, 1), (2.5, -1), (3.5, -1), (5.5, -1), (6, 1), (7, 1)]
            sage: pa.q_parenthesis(7)
            (2, 1)
        """
        parenthesis = maximum = shift = 0
        for s, paren in self.q_christol_sorting(q):
            parenthesis += paren
            if parenthesis > maximum:
                maximum = parenthesis
                shift = s
        return shift, maximum

    def q_parenthesis_criterion(self, q):
        r"""
        Return ``True`` if in each prefix of the list
        ``self.q_christol_sorting(q)`` there are at least as many pairs with
        second entry -1 as pairs with second entry +1. Return ``False``
        otherwise.

        INPUT:

        - ``q`` -- an integer

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.q_christol_sorting(7)
            [(2, 1), (2.5, -1), (3.5, -1), (5.5, -1), (6, 1), (7, 1)]
            sage: pa.q_parenthesis_criterion(7)
            False
            sage: pa.q_christol_sorting(61)
            [(15.5, -1), (20.5, -1), (25, 1), (30.5, -1), (37, 1), (61, 1)]
            sage: pa.q_parenthesis_criterion(61)
            True
        """
        parenthesis = 0
        for _, paren in self.q_christol_sorting(q):
            parenthesis += paren
            if parenthesis > 0:
                return False
        return parenthesis <= 0

    def q_interlacing_number(self, q):
        r"""
        Return the number of pairs in the list ``self.q_christol_sorting(q)``
        with second entry 1, that were preceded by a pair with second entry
        -1.

        INPUT:

        - ``q`` -- an integer.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.q_christol_sorting(7)
            [(2, 1), (2.5, -1), (3.5, -1), (5.5, -1), (6, 1), (7, 1)]
            sage: pa.q_interlacing_number(7)
            1
        """
        interlacing = 0
        previous_paren = 1
        for _, paren in self.q_christol_sorting(q):
            if paren == 1 and previous_paren == -1:
                interlacing += 1
            previous_paren = paren
        return interlacing

    def remove_positive_integer_differences(self):
        r"""
        Return parameters, where pairs consisting of a top parameter
        and a bottom parameter with positive integer differences are
        removed, starting with pairs of minimal positive integer
        difference.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([5/2, -1/2, 5/3], [3/2, 1/3])
            sage: pa
            ((-1/2, 5/3, 5/2), (1/3, 3/2, 1))
            sage: pa.remove_positive_integer_differences()
            ((-1/2, 5/3), (1/3, 1))

        The choice of which pair with integer differences to remove first
        is important::

            sage: pa = HypergeometricParameters([4, 2, 1/2], [1, 3])
            sage: pa
            ((1/2, 2, 4), (1, 3, 1))
            sage: pa.remove_positive_integer_differences()
            ((1/2,), (1,))
        """
        differences = []
        top = list(self.top)
        bottom = list(self.bottom)
        for i in range(len(top)):
            for j in range(len(bottom)):
                diff = top[i] - bottom[j]
                if diff in ZZ and diff > 0:
                    differences.append((diff, i, j))
        for _, i, j in sorted(differences):
            if top[i] is not None and bottom[j] is not None:
                top[i] = None
                bottom[j] = None
        return HypergeometricParameters(top, bottom, add_one=False)

    def has_negative_integer_differences(self):
        r"""
        Return ``True`` if there exists a pair of a top parameter and a bottom
        parameter, such that the top one minus the bottom one is a negative integer;
        return ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.has_negative_integer_differences()
            False

        ::

            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/2])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/2, 1))
            sage: pa.has_negative_integer_differences()
            True
        """
        return any(a - b in ZZ and a < b for a in self.top for b in self.bottom)

    def shift(self, s):
        r"""
        Return the parameters obtained by adding s to each of them.

        INPUT:

        - ``s`` -- a rational number

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.shift(2)
            ((1, 9/4, 7/3, 5/2), (12/5, 13/5, 3, 1))
        """
        top = [a+s for a in self.top]
        bottom = [b+s for b in self.bottom]
        return HypergeometricParameters(top, bottom, add_one=False)

    def decimal_part(self):
        r"""
        Return the parameters obtained by taking the decimal part of each of
        the parameters, where integers are assigned 1 instead of 0.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([5/4, 1/3, 2], [2/5, -2/5])
            sage: pa
            ((1/3, 5/4, 2), (-2/5, 2/5, 1))
            sage: pa.decimal_part()
            ((1/4, 1/3, 1), (2/5, 3/5, 1))
        """
        top = [1 + a - ceil(a) for a in self.top]
        bottom = [1 + b - ceil(b) for b in self.bottom]
        return HypergeometricParameters(top, bottom, add_one=False)

    def valuation_position(self, p, drift=0):
        r"""
        If the `h_k`s are the coefficients of the hypergeometric
        series corresponding to these parameters and `\delta` is
        the drift, return the smallest value of

        .. MATH::

            \text{val}_p(h_k) + \delta k

        and the first index `k` where this minimum is reached.

        INPUT:

        - ``p`` -- a prime number

        - ``drift`` -- a rational number (default: ``0``)

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/5, 1/5, 1/5], [1/3, 3^10/5])
            sage: pa.valuation_position(3)
            (-9, 1)

        When the relevant sequence is not bounded from below, the
        tuple ``(-Infinity, None)`` is returned::

            sage: pa.valuation_position(5)
            (-Infinity, None)

        An example with a drift::

            sage: pa.valuation_position(3, drift=-7/5)
            (-54/5, 7)
        """
        # We treat the case of parameters having p in the denominator
        # When it happens, we remove the parameter and update the drift accordingly
        top = []
        for a in self.top:
            if a in ZZ and a <= 0:
                raise NotImplementedError  # TODO
            v = a.valuation(p)
            if v < 0:
                drift += v
            else:
                top.append(a)
        bottom = []
        for b in self.bottom:
            v = b.valuation(p)
            if v < 0:
                drift -= v
            else:
                bottom.append(b)

        # We check that we are inside the disk of convergence
        diff = len(top) - len(bottom)
        growth = (p-1)*drift + diff
        if (growth, drift) < (0, 0):
            return -infinity, None

        # Main part: computation of the valuation
        # We use Christol's formula (see Lemma 3.1.10 of [CFVM2025])
        # with modification in order to take the drift into account
        TSR = TropicalSemiring(QQ)
        n = len(top) + len(bottom) + 1
        parameters = HypergeometricParameters(top, bottom)
        order = IntegerModRing(parameters.d)(p).multiplicative_order()
        valuation = position = ZZ(0)
        valfinal = breaks = None
        indices = {}
        count = 0
        q = 1
        TM = identity_matrix(TSR, n)
        while True:
            # We take into account the contribution of V({k/p^r}, p^r).
            # We represent the partial sum until r by the table new_breaks.
            # Its entries are triples (valuation, position, parameter):
            # - parameter is the parameter corresponding to a point of
            #   discontinuity of the last summand V({k/p^r}, p^r)
            # - valuation is the minimum of the partial sum on the
            #   range starting at this discontinuity point (included)
            #   and ending at the next one (excluded)
            # - position is the first position the minimum is reached
            # (The dictionary new_indices allows for finding rapidly
            # an entry in new_breaks with a given parameter.)
            # The table breaks and the dictionary indices correspond
            # to the same data for r-1.

            pq = p * q

            # We compute the point of discontinuity of V({k/p^r}, p^r)
            # and store them in AB
            # Each entry of AB has the form (x, dw, parameter) where:
            # - x is the position of the discontinuity point
            # - dw is the *opposite* of the jump of V({k/p^r}, p^r)
            #   at this point
            #   (taking the opposite is useful for sorting reasons)
            # - parameter is the corresponding parameter
            A = [(1 + (-a) % pq, -1, a) for a in top]
            B = [(1 + (-b) % pq, 1, b) for b in bottom]
            AB = [(0, 0, 0)] + sorted(A + B) + [(pq, None, 0)]

            # We compute new_breaks
            new_breaks = []
            new_indices = {}
            w = 0
            TMstep = matrix(TSR, n)
            for i in range(n):
                x, dw, param = AB[i]    # discontinuity point
                y, _, right = AB[i+1]   # next discontinuity point
                w -= dw   # the value of V({k/p^r}, p^r) on this interval
                new_indices[param] = len(new_breaks)
                if x == y:
                    # Case of empty interval
                    val = infinity
                    pos = None
                elif breaks is None:
                    # Case r = 1
                    if drift < 0:
                        pos = y - 1
                    else:
                        pos = x
                    val = drift * pos
                else:
                    # Case r > 1
                    # The variable complete stores whether the interval
                    # [x,y] covers all [0, p^(r-1)) modulo p^(r-1)
                    complete = (y - x >= q)
                    if complete and drift < 0:
                        interval = ((y-1) // q) - 1
                        j = j0 = indices[right]
                    else:
                        interval = max(0, (x-1) // q)
                        j = j0 = indices[param]
                    val = infinity
                    while True:
                        valj, posj, paramj = breaks[j]
                        valj += drift * interval
                        TMstep[i,j] = TSR(drift*interval + w)
                        if valj < val:
                            val = valj
                            pos = posj + q * interval
                        j += 1
                        if j >= len(breaks):
                            j = 0
                            interval += 1
                        if j == j0 or (not complete and breaks[j][2] == right):
                            break
                new_breaks.append((val + w, pos, param))

            # The halting criterion
            minimum = min(new_breaks)
            valuation, position, _ = minimum
            if drift >= 0 and q > parameters.bound:
                if growth == 0:
                    if count < order:
                        TM = TMstep * TM
                    count += 1
                    if count == order:
                        try:
                            TM = TM.weak_transitive_closure()
                        except ValueError:
                            return -infinity, None
                        valfinal = min(TM[i,j].lift() + new_breaks[j][0]
                                       for i in range(n) for j in range(n))
                    if valuation == valfinal:
                        return valuation, position
                elif drift >= len(top) and minimum == new_breaks[0]:
                    return valuation, position

            # We update the values for the next r
            q = pq
            drift = p*drift + diff
            breaks = new_breaks
            indices = new_indices

    def dwork_image(self, p):
        r"""
        Return the parameters obtained by applying the Dwork map to each of
        the parameters. The Dwork map `D_p(x)` of a p-adic integer x is defined
        as the unique p-adic integer such that `p D_p(x) - x` is a nonnegative
        integer smaller than p.

        INPUT:

        - ``p`` -- a prime number

        EXAMPLE::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.dwork_image(7)
            ((1/3, 1/2, 3/4), (1/5, 4/5, 1))

        If `p` is not coprime to the common denominators of the parameters,
        a ``ValueError`` is raised::

            sage: pa.dwork_image(3)
            Traceback (most recent call last):
            ...
            ValueError: denominators of parameters are not coprime to p
        """
        try:
            top = [(a + (-a) % p) / p for a in self.top]
            bottom = [(b + (-b) % p) / p for b in self.bottom]
        except ZeroDivisionError:
            raise ValueError("denominators of parameters are not coprime to p")
        return HypergeometricParameters(top, bottom, add_one=False)

    def frobenius_order(self, p):
        r"""
        Return the Frobenius order of the hypergeometric function with this set
        of parameters, that is the order of the Dwork map acting on the decimal
        parts of the parameters.

        INPUT:

        - ``p`` -- a prime number

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import HypergeometricParameters
            sage: pa = HypergeometricParameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.frobenius_order(7)
            2
        """
        param = self.decimal_part()
        iter = param.dwork_image(p)
        i = 1
        while param != iter:
            iter = iter.dwork_image(p)
            i += 1
        return i
