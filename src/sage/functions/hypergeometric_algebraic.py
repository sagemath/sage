r"""
Algebraic properties of hypergeometric functions.

[Tutorial]

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

import operator

from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name
from sage.misc.cachefunc import cached_method

from sage.misc.misc_c import prod
from sage.misc.functional import log
from sage.functions.other import ceil
from sage.functions.hypergeometric import hypergeometric
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.matrix.constructor import matrix

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element import coerce_binop
from sage.structure.category_object import normalize_names

from sage.categories.action import Action
from sage.categories.pushout import pushout
from sage.categories.map import Map
from sage.categories.finite_fields import FiniteFields

from sage.matrix.special import companion_matrix
from sage.matrix.special import identity_matrix

from sage.symbolic.ring import SR
from sage.combinat.subset import Subsets
from sage.rings.infinity import infinity
from sage.sets.primes import Primes
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.padics.factory import Qp
from sage.rings.number_field.number_field import CyclotomicField

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing


# Helper functions
##################

def insert_zeroes(P, n):
    cs = P.list()
    coeffs = n * len(cs) * [0]
    for i in range(len(cs)):
        coeffs[n*i] = cs[i]
    return P.parent()(coeffs)


def kernel(M, repeat=2):
    n = M.nrows()
    m = M.ncols()
    if n > m + 1:
        raise RuntimeError
    if n <= m:
        K = M.base_ring().base_ring()
        for _ in range(repeat):
            a = K.random_element()
            Me = matrix(n, m, [f(a) for f in M.list()])
            if Me.rank() == n:
                return
    for J in Subsets(range(m), n-1):
        MJ = M.matrix_from_columns(J)
        minor = MJ.delete_rows([0]).determinant()
        if minor.is_zero():
            continue
        ker = [minor]
        for i in range(1, n):
            minor = MJ.delete_rows([i]).determinant()
            ker.append((-1)**i * minor)
        Z = matrix(ker) * M
        if not Z.is_zero():
            return
        g = ker[0].leading_coefficient() * gcd(ker)
        ker = [c//g for c in ker]
        return ker


# Parameters of hypergeometric functions
########################################

class Parameters():
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: type(pa)
            <class 'sage.functions.hypergeometric_algebraic.Parameters'>

        By default, parameters are sorted, duplicates are removed and
        a trailing `1` is added to the bottom parameters::

            sage: Parameters([1/2, 1/3, 2/3], [2/3])
            ((1/3, 1/2), (1,))

        We can avoid adding the trailing `1` by passing ``add_one=False``::

            sage: Parameters([1/2, 1/3, 2/3], [2/3], add_one=False)
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa  # indirect doctest
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
        """
        return "(%s, %s)" % (self.top, self.bottom)

    def __hash__(self):
        return hash((self.top, self.bottom))

    def __eq__(self, other):
        return (isinstance(other, Parameters)
            and self.top == other.top and self.bottom == other.bottom)

    def is_balanced(self):
        r"""
        Return ``True`` if there are as many top parameters as bottom
        parameters; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.is_balanced()
            True

        ::

            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5], add_one=False)
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/3, 2/3], [1/2])
            sage: pa
            ((1/3, 2/3), (1/2, 1))
            sage: pa.interlacing_criterion(1)
            True
            sage: pa.interlacing_criterion(5)
            True

        ::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/8, 3/8, 5/8], [1/4, 1/2])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([5/2, -1/2, 5/3], [3/2, 1/3])
            sage: pa
            ((-1/2, 5/3, 5/2), (1/3, 3/2, 1))
            sage: pa.remove_positive_integer_differences()
            ((-1/2, 5/3), (1/3, 1))

        The choice of which pair with integer differences to remove first
        is important::

            sage: pa = Parameters([4, 2, 1/2], [1, 3])
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
        return Parameters(top, bottom, add_one=False)

    def has_negative_integer_differences(self):
        r"""
        Return ``True`` if there exists a pair of a top parameter and a bottom
        parameter, such that the top one minus the bottom one is a negative integer;
        return ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.has_negative_integer_differences()
            False

        ::

            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/2])
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

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: pa
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: pa.shift(2)
            ((1, 9/4, 7/3, 5/2), (12/5, 13/5, 3, 1))
        """
        top = [a+s for a in self.top]
        bottom = [b+s for b in self.bottom]
        return Parameters(top, bottom, add_one=False)

    def decimal_part(self):
        r"""
        Return the parameters obtained by taking the decimal part of each of
        the parameters, where integers are assigned 1 instead of 0.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([5/4, 1/3, 2], [2/5, -2/5])
            sage: pa
            ((1/3, 5/4, 2), (-2/5, 2/5, 1))
            sage: pa.decimal_part()
            ((1/4, 1/3, 1), (2/5, 3/5, 1))
        """
        top = [1 + a - ceil(a) for a in self.top]
        bottom = [1 + b - ceil(b) for b in self.bottom]
        return Parameters(top, bottom, add_one=False)

    def valuation_position(self, p, drift=0):
        top = []
        for a in self.top:
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
        diff = len(top) - len(bottom)
        if ((p-1)*drift + diff, drift) < (0, 0):
            return -infinity, None

        parameters = Parameters(top, bottom)
        order = IntegerModRing(parameters.d)(p).multiplicative_order()
        q = 1
        valuation = position = 0
        breaks = [(0, 0, 0)]
        indices = None
        count = 0
        while True:
            pq = p * q
            A = [(1 + (-a) % pq, -1, a) for a in top]
            B = [(1 + (-b) % pq, 1, b) for b in bottom]
            AB = [(0, 0, 0)] + sorted(A + B) + [(pq, None, None)]
            new_breaks = []
            new_indices = {}
            w = 0
            for i in range(len(AB) - 1):
                x, dw, param = AB[i]
                y, _, right = AB[i+1]
                w -= dw
                new_indices[param] = len(new_breaks)
                complete = (y - x >= q)
                if complete and drift < 0:
                    interval = (y // q) - 1
                else:
                    interval = x // q
                if x == y:
                    val = infinity
                    pos = x
                elif indices is None:
                    val = drift * interval
                    pos = q * interval
                else:
                    val = infinity
                    if complete and drift < 0:
                        j = j0 = indices[right]
                    else:
                        j = j0 = indices[param]
                    while True:
                        valj, posj, paramj = breaks[j]
                        valj += drift * interval
                        if valj < val:
                            val = valj
                            pos = posj + q * interval
                        j += 1
                        if j >= len(breaks):
                            if right is None:
                                break
                            j = 0
                            interval += 1
                        if (not complete and paramj == right) or (complete and j == j0):
                            break
                new_breaks.append((val + w, pos, param))
            breaks = new_breaks
            indices = new_indices
            minimum = min(breaks)
            if drift >= 0 and q > parameters.bound:
                # Not sure at all about this criterion
                if minimum == breaks[0] and minimum[0] == valuation and minimum[1] == position:
                    count += 1
                    if count >= order:
                        return valuation, position
                else:
                    return -infinity, None
            q = pq
            drift = p*drift + diff
            valuation, position, _ = minimum

    def dwork_image(self, p):
        r"""
        Return the parameters obtained by applying the Dwork map to each of
        the parameters. The Dwork map `D_p(x)` of a p-adic integer x is defined
        as the unique p-adic integer such that `p D_p(x) - x` is a nonnegative
        integer smaller than p.

        INPUT:

        - ``p`` -- a prime number

        EXAMPLE::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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
        return Parameters(top, bottom, add_one=False)

    def frobenius_order(self, p):
        r"""
        Return the Frobenius order of the hypergeometric function with this set
        of parameters, that is the order of the Dwork map acting on the decimal
        parts of the parameters.

        INPUT:

        - ``p`` -- a prime number

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: pa = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
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


# Hypergeometric functions
##########################

# Do we want to implement polynomial linear combinaison
# of hypergeometric functions?
# Advantages:
#  . reductions mod p of hypergeometric functions have this form in general
#  . many methods can be extended to this context
# Difficulty:
#  . not sure we can handle easily simplifications!

class HypergeometricAlgebraic(Element):
    r"""
    Class for hypergeometric functions over arbitrary base rings.
    """
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        r"""
        Initialize this hypergeometric function.

        INPUT:

        - ``parent`` -- the parent of this function

        - ``arg1``, ``arg2`` -- arguments defining this hypergeometric
          function, they can be:
          - the top and bottom paramters
          - a hypergeometric function and ``None``
          - an instance of the class :class:`Parameters` and ``None``

        - ``scalar`` -- an element in the base ring, the scalar by
          which the hypergeometric function is multiplied

        TESTS::

            sage: S.<x> = QQ[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: type(h)
            <class 'sage.functions.hypergeometric_algebraic.HypergeometricFunctions.element_class'>
            sage: TestSuite(h).run()
        """
        Element.__init__(self, parent)
        base = parent.base_ring()
        if scalar is None:
            scalar = base.one()
        else:
            scalar = base(scalar)
        if scalar == 0:
            parameters = None
        elif isinstance(arg1, HypergeometricAlgebraic):
            parameters = arg1._parameters
            scalar *= base(arg1._scalar)
        elif isinstance(arg1, Parameters):
            parameters = arg1
        else:
            parameters = Parameters(arg1, arg2)
        char = self.parent()._char
        if scalar:
            if any(b in ZZ and b < 0 for b in parameters.bottom):
                raise ValueError("the parameters %s do not define a hypergeometric function" % parameters)
            if char > 0:
                val, _ = parameters.valuation_position(char)
                if val < 0:
                    raise ValueError("the parameters %s do not define a hypergeometric function in characteristic %s" % (parameters, char))
        self._scalar = scalar
        self._parameters = parameters
        self._coeffs = [scalar]
        self._char = char

    def __hash__(self):
        return hash((self.base_ring(), self._parameters, self._scalar))

    def __eq__(self, other):
        return (isinstance(other, HypergeometricAlgebraic)
            and self.base_ring() is other.base_ring()
            and self._parameters == other._parameters
            and self._scalar == other._scalar)

    def _repr_(self):
        if self._parameters is None:
            return "0"
        scalar = self._scalar
        if scalar == 1:
            s = ""
        elif scalar._is_atomic():
            scalar = str(scalar)
            if scalar == "-1":
                s = "-"
            else:
                s = scalar + "*"
        else:
            s = "(%s)*" % scalar
        s += "hypergeometric(%s, %s, %s)" % (self.top(), self.bottom(), self.parent().variable_name())
        return s

    def _latex_(self):
        if self._parameters is None:
            return "0"
        scalar = self._scalar
        if scalar == 1:
            s = ""
        elif scalar._is_atomic():
            scalar = latex(scalar)
            if scalar == "-1":
                s = "-"
            else:
                s = scalar
        else:
            s = r"\left(%s\right)" % scalar
        top = self.top()
        bottom = self.bottom()
        s += r"\,_{%s} F_{%s} " % (len(top), len(bottom))
        s += r"\left(\begin{matrix} "
        s += ",".join(latex(a) for a in top)
        s += r"\\"
        s += ",".join(latex(b) for b in bottom)
        s += r"\end{matrix}; %s \right)" % self.parent().latex_variable_name()
        return s

    def base_ring(self):
        r"""
        Return the ring over which this hypergeometric function is defined.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.base_ring()
            Rational Field

        ::

            sage: T.<y> = Qp(5)[]
            sage: g = hypergeometric([1/3, 2/3], [1/2], y)
            sage: g.base_ring()
            5-adic Field with capped relative precision 20

        ::

            sage: U.<z> = GF(5)[]
            sage: h = hypergeometric([1/3, 2/3], [1/2], z)
            sage: h.base_ring()
            Finite Field of size 5

        ::

            sage: V.<w> = CC[]
            sage: k = hypergeometric([1/3, 2/3], [1/2], w)
            sage: k.base_ring()
            Complex Field with 53 bits of precision
        """
        return self.parent().base_ring()

    def top(self):
        r"""
        Return the top parameters of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.top()
            (1/3, 2/3)
        """
        return self._parameters.top

    def bottom(self):
        r"""
        Return the bottom parameters of this hypergeometric function (excluding
        the extra ``1``).

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.bottom()
            (1/2,)
        """
        return self._parameters.bottom[:-1]

    def scalar(self):
        r"""
        Return the scalar of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.scalar()
            1
            sage: g = 4*f
            sage: g.scalar()
            4
        """
        return self._scalar

    def change_ring(self, R):
        r"""
        Return this hypergeometric function with changed base ring.

        INPUT:

        - ``R`` -- a commutative ring

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.base_ring()
            Rational Field
            sage: g = f.change_ring(Qp(5))
            sage: g.base_ring()
            5-adic Field with capped relative precision 20
        """
        H = self.parent().change_ring(R)
        return H(self._parameters, None, self._scalar)

    def change_variable_name(self, name):
        r"""
        Return this hypergeometric function with changed variable name

        INPUT:

        - ``name`` -- a string, the new variable name

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: T.<y> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f
            hypergeometric((1/3, 2/3), (1/2,), x)
            sage: g = f.change_variable_name('y')
            sage: g
            hypergeometric((1/3, 2/3), (1/2,), y)
        """
        H = self.parent().change_variable_name(name)
        return H(self._parameters, None, self._scalar)

    def _add_(self, other):
        r"""
        Return the (formal) sum of the hypergeometric function
        and ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f + g
            3/2*hypergeometric((1/3, 2/3), (1/2,), x)
            sage: f + h
            hypergeometric((1/3, 2/3), (1/2,), x) + hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: f + cos(x)
            cos(x) + hypergeometric((1/3, 2/3), (1/2,), x)
        """
        if self._parameters is None:
            return other
        if isinstance(other, HypergeometricAlgebraic):
            if other._parameters is None:
                return self
            if self._parameters == other._parameters:
                scalar = self._scalar + other._scalar
                return self.parent()(self._parameters, scalar=scalar)
        return SR(self) + SR(other)

    def _neg_(self):
        r"""
        Return the negative of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = 2*hypergeometric([1/3, 2/3], [1/2], x)
            sage: -f
            -2*hypergeometric((1/3, 2/3), (1/2,), x)
        """
        if self._parameters is None:
            return self
        return self.parent()(self._parameters, scalar=-self._scalar)

    def _sub_(self, other):
        r"""
        Return the (formal) difference of the hypergeometric function
        with ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function or a formal expression

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f - g
            1/2*hypergeometric((1/3, 2/3), (1/2,), x)
            sage: f - h
            hypergeometric((1/3, 2/3), (1/2,), x) - hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: f - sin(x)
            hypergeometric((1/3, 2/3), (1/2,), x) - sin(x)
        """
        if self._parameters is None:
            return other
        if isinstance(other, HypergeometricAlgebraic):
            if other._parameters is None:
                return self
            if self._parameters == other._parameters:
                scalar = self._scalar - other._scalar
                return self.parent()(self._parameters, scalar=scalar)
        return SR(self) - SR(other)

    def _mul_(self, other):
        r"""
        Return the (formal) product of the hypergeometric function
        and ``other``

        INPUT:

        - ``other`` -- a hypergeometric function or a formal expression

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f*g
            1/2*hypergeometric((1/3, 2/3), (1/2,), x)^2
            sage: f*h
            hypergeometric((1/3, 2/3), (1/2,), x)*hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: sin(x)*f + x
            hypergeometric((1/3, 2/3), (1/2,), x)*sin(x) + x
        """
        return SR(self) * SR(other)

    def __call__(self, x):
        r"""
        Return the value of this hypergeometric function at ``x``.

        INPUT:

        - ``x`` -- an element

        EXAMPLES::

            sage: S.<x> = RR[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f(0.5)
            1.36602540378444

        ::

            sage: g = 2*f
            sage: g(0.2)
            2.20941633798502
        """
        scalar = self._scalar
        if scalar == 0:
            return self.base_ring().zero()
        X = SR('X')
        h = hypergeometric(self.top(), self.bottom(), X)
        if scalar != 1:
            h *= scalar
        return h(X=x)

    def _compute_coeffs(self, prec):
        r"""
        Compute the coefficients of the series representation of this
        hypergeometric function up to a given precision, and store
        them in ``self._coeffs``.

        INPUT:

        - ``prec`` -- a positive integer

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f._coeffs
            [1]
            sage: f._compute_coeffs(3)
            sage: f._coeffs
            [1, 4/9, 80/243]
        """
        coeffs = self._coeffs
        start = len(coeffs) - 1
        c = coeffs[-1]
        for i in range(start, prec - 1):
            for a in self._parameters.top:
                c *= a + i
            for b in self._parameters.bottom:
                c /= b + i
            coeffs.append(c)

    def power_series(self, prec):
        r"""
        Return the power series representation of this hypergeometric
        function up to a given precision.

        INPUT:

        - ``prec`` -- a positive integer

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.power_series(3)
            1 + 4/9*x + 80/243*x^2 + O(x^3)
        """
        S = self.parent().power_series_ring()
        self._compute_coeffs(prec)
        return S(self._coeffs, prec=prec)

    def shift(self, s):
        r"""
        Return this hypergeometric function, where each parameter
        (including the additional ``1`` as a bottom parameter) is
        increased by ``s``.

        INPUT:

        - ``s`` -- a rational number

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = f.shift(3/2)
            sage: g
            hypergeometric((1, 11/6, 13/6), (2, 5/2), x)
        """
        return self.parent()(self._parameters.shift(s), scalar=self._scalar)

    @coerce_binop
    def hadamard_product(self, other):
        r"""
        Return the hadamard product of the hypergeometric function
        and ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = 1/2*hypergeometric([1/5, 2/5], [3/5], x)
            sage: f.hadamard_product(h)
            1/2*hypergeometric((1/5, 1/3, 2/5, 2/3), (1/2, 3/5, 1), x)
        """
        if self._scalar == 0:
            return self
        if other._scalar == 0:
            return other
        top = self.top() + other.top()
        bottom = self._parameters.bottom + other.bottom()
        scalar = self._scalar * other._scalar
        return self.parent()(top, bottom, scalar=scalar)

    def _div_(self, other):
        r"""
        Return the (formal) quotient of the hypergeometric function
        and ``other``

        INPUT:

        - ``other`` -- a hypergeometric function or a formal expression

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f/g
            2
            sage: f/h
            hypergeometric((1/3, 2/3), (1/2,), x)/hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: f/sin(x) + x
            x + hypergeometric((1/3, 2/3), (1/2,), x)/sin(x)
        """
        return SR(self) / SR(other)

    def denominator(self):
        r"""
        Return the smallest common denominator of the parameters.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.denominator()
            6
        """
        return self._parameters.d

    def differential_operator(self, var='d'):
        r"""
        Return the hypergeometric differential operator that annihilates
        this hypergeometric function as an Ore polynomial in the variable
        ``var``.

        INPUT:

        - ``var`` -- a string (default: ``d``), the variable name of
          the derivation

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.differential_operator(var='D')
            (-x^2 + x)*D^2 + (-2*x + 1/2)*D - 2/9

        Note that this does not necessarily give the minimal differential
        operator annihilating this hypergeometric function: in the example
        below, this method returns an operator of order `3` where `g` is
        solution of a differential equation of order `2`::

            sage: g = hypergeometric([1/3, 2/3, 6/5], [1/5, 1/2], x)
            sage: L = g.differential_operator()
            sage: L.degree()
            3
            sage: gs = g.power_series(100)
            sage: (72*x^3 - 234*x^2 + 162*x)*gs.derivative(2) + (144*x^2 - 450*x + 81)*gs.derivative() + (16*x - 216)*gs
            O(x^99)
        """
        S = self.parent().polynomial_ring()
        x = S.gen()
        D = OrePolynomialRing(S, S.derivation(), names=var)
        if self._scalar == 0:
            return D.one()
        t = x * D.gen()
        A = D.one()
        for a in self._parameters.top:
            A *= t + S(a)
        B = D.one()
        for b in self._parameters.bottom:
            B *= t + S(b-1)
        L = B - x*A
        return D([c//x for c in L.list()])

    def derivative(self):
        r"""
        Return the derivative of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.derivative()
            4/9*hypergeometric((4/3, 5/3), (3/2,), x)
        """
        top = [a+1 for a in self.top()]
        bottom = [b+1 for b in self.bottom()]
        scalar = prod(self._parameters.top) / prod(self._parameters.bottom)
        scalar = self.base_ring()(scalar) * self._scalar
        return self.parent()(top, bottom, scalar)


# Over the rationals

class HypergeometricAlgebraic_QQ(HypergeometricAlgebraic):
    def __mod__(self, p):
        r"""
        Return the reduction of the hypergeometric function modulo ``p``.

        INPUT:

        - ``p`` -- a prime number.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = f % 5
            sage: g
            hypergeometric((1/3, 2/3), (1/2,), x)
            sage: g.base_ring()
            Finite Field of size 5
        """
        k = FiniteField(p)
        val = self._scalar.valuation(p)
        if val == 0:
            return self.change_ring(k)
        h = self.change_ring(Qp(p, 1))
        return h.residue()

    def valuation(self, p):
        r"""
        Return the p-adic valuation of this hypergeometric function, i.e., the
        maximal s, such that p^(-s) times this hypergeometric function has
        p-integral coefficients.

        INPUT:

        - ``p`` -- a prime number

        EXAMPLES::

            sage: S.<x> = QQ[x]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.valuation(5)
            0
            sage: g = 5*f
            sage: g.valuation(5)
            1
        """
        val, _ = self._parameters.valuation_position(p)
        return val + self._scalar.valuation(p)

    def has_good_reduction(self, p):
        r"""
        Return ``True`` if the p-adic valuation of this hypergeometric function
        is non-negative, i.e., if its reduction modulo ``p`` is well-defined.pAdicGeneric

        INPUT:

        - ``p`` -- a prime number

        EXAMPLES::

            sage: S.<x> = QQ[x]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.valuation(5)
            0
            sage: f.has_good_reduction(5)
            True
            sage: g = 1/5*f
            sage: g.has_good_reduction(5)
            False
        """
        return self.valuation(p) >= 0

    def good_reduction_primes(self):
        r"""
        Return the set of prime numbers modulo which this hypergeometric
        function can be reduced, i.e., the p-adic valuation is positive.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.good_reduction_primes()
            Set of all prime numbers with 3 excluded: 2, 5, 7, 11, ...

        ALGORITHM:

        We rely on Christol's criterion ([Chr1986]_, Prop. 1) for globally
        bounded hypergeometric function, from which a criterion can be deduced
        modulo which primes a hypergeometric function can be reduced
        ([CFV2025]_, Thm. 3.1.3). For small primes `p`, we compute the `p`-adic
        valuation of the hypergeometric function individually.
        """
        params = self._parameters
        d = params.d

        # We check the parenthesis criterion for c=1
        if not params.parenthesis_criterion(1):
            return Primes(modulus=0)

        # We check the parenthesis criterion for other c
        # and derive congruence classes with good reduction
        goods = {c: None for c in range(d) if d.gcd(c) == 1}
        goods[1] = True
        for c in goods.keys():
            if goods[c] is not None:
                continue
            cc = c
            goods[c] = True
            while cc != 1:
                if goods[cc] is False or not params.parenthesis_criterion(cc):
                    goods[c] = False
                    break
                cc = (cc * c) % d
            if goods[c]:
                cc = c
                while cc != 1:
                    goods[cc] = True
                    cc = (cc * c) % d

        # We treat exceptional primes
        bound = params.bound
        exceptions = {}
        for p in Primes():
            if p > bound:
                break
            if d % p == 0 and self.valuation(p) >= 0:
                exceptions[p] = True
            if d % p == 0 or not goods[p % d]:
                continue
            if self.valuation(p) < 0:
                exceptions[p] = False

        goods = [c for c, v in goods.items() if v]
        return Primes(modulus=d, classes=goods, exceptions=exceptions)

    def is_algebraic(self):
        if any(a in ZZ and a <= 0 for a in self.top()):
            return True
        if not self._parameters.is_balanced():
            return False
        simplified_parameters = self._parameters.remove_positive_integer_differences()
        if simplified_parameters.has_negative_integer_differences():
            return False
        d = simplified_parameters.d
        return all(simplified_parameters.interlacing_criterion(c)
                   for c in range(d) if d.gcd(c) == 1)

    def is_globally_bounded(self, include_infinity=True):
        if include_infinity and len(self.top()) > len(self.bottom()) + 1:
            return False
        d = self.denominator()
        for c in range(d):
            if d.gcd(c) == 1:
                if not self._parameters.parenthesis_criterion(c):
                    return False
        return True

    def p_curvature_ranks(self):
        raise NotImplementedError

    def monodromy(self, x=0, var='z'):
        params = self._parameters
        if not params.is_balanced():
            raise ValueError("hypergeometric equation is not Fuchsian")
        d = params.d
        K = CyclotomicField(d, names=var)
        z = K.gen()
        S = PolynomialRing(K, names='X')
        X = S.gen()
        if x == 0:
            B = prod(X - z**(b*d) for b in params.bottom)
            return companion_matrix(B, format='right').inverse()
        elif x == 1:
            A = prod(X - z**(a*d) for a in params.top)
            B = prod(X - z**(b*d) for b in params.bottom)
            return companion_matrix(A, format='right').inverse() * companion_matrix(B, format='right')
        elif x is infinity:
            A = prod(X - z**(a*d) for a in params.top)
            return companion_matrix(A, format='right')
        else:
            n = len(params.top)
            return identity_matrix(QQ, n)

    def is_maximum_unipotent_monodromy(self):
        return all(b in ZZ for b in self.bottom())

    is_mum = is_maximum_unipotent_monodromy


# Over the p-adics

class HypergeometricAlgebraic_padic(HypergeometricAlgebraic):
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar)
        K = self.base_ring()
        self._p = K.prime()
        self._e = K.e()

    def residue(self):
        k = self.base_ring().residue_field()
        if self._scalar.valuation() == 0:
            return self.change_base(k)
        val, pos = self._parameters.valuation_position(self._p)
        if val < 0:
            raise ValueError("bad reduction")
        if val > 0:
            H = self.parent().change_ring(k)
            return H(self._parameters, scalar=0)
        raise NotImplementedError("the reduction is not a hypergeometric function")
        # In fact, it is x^s * h[s] * h, with
        # . s = pos
        # . h = self.shift(s)

    def _val_pos(self):
        p = self._p
        d = self.denominator()
        parameters = self._parameters
        if d.gcd(p) > 1:
            T = []
            B = []
            difference = 0
            for j in self.top():
                d = j.denominator().valuation(p)
                difference += d
                if not d:
                    T += [j]
            for j in self.bottom():
                d = j.denominator().valuation(p)
                difference -= d
                if not d:
                    B += [j]
            if difference > 0:
                return -infinity, None
            if difference < 0:
                _, prec = self.parent()(T, B, self.scalar())._val_pos()
                #This is the case when the p-adic valuation of the coefficients
                #goes to +infinity, but if you remove all the coefficients with p
                #in the denominator it goes to -infinity.
                #The following claim is almost surely wrong.
                if prec is None:
                    prec = self._parameters.bound
                #Here I just check the valuation of the first few coefficients.
                #There should be something better using q_g:parenthesis.
                L = self.change_ring(Qp(p, 1)).power_series(prec+1).coefficients()
                val = + infinity
                pos = 0
                for j in range(len(L)):
                    v = L[j].valuation()
                    if v < val:
                        val = v
                        pos = j
                return val, pos
            if difference == 0:
                return self.parent()(T, B, self.scalar())._val_pos()
        u = 1
        if not parameters.parenthesis_criterion(u):
            return -infinity, None
        u = p % d
        while u != 1:
            if not parameters.parenthesis_criterion(u):
                return -infinity, None
            u = p*u % d
        # From here, it is absolutely conjectural!
        # ... and probably not quite correct
        val = self._scalar.valuation()
        pos = 0
        q = 1
        while True:
            s, v = parameters.q_parenthesis(p)
            if v == 0:
                break
            val -= self._e * v
            pos += q*s
            q *= p
            parameters = parameters.shift(s).dwork_image(p)
        return val, pos

    def log_radius_of_convergence(self):
        p = self._p
        step = self._e / (p - 1)
        log_radius = 0
        for a in self._parameters.top:
            v = a.valuation(p)
            if v < 0:
                log_radius += v
            else:
                log_radius += step
        for b in self._parameters.bottom:
            v = b.valuation(p)
            if v < 0:
                log_radius -= v
            else:
                log_radius -= step
        return log_radius

    def valuation(self, log_radius=0):
        drift = -log_radius / self._e
        val, _ = self._parameters.valuation_position(self._p, drift)
        return val

    def newton_polygon(self, log_radius):
        raise NotImplementedError

    def tate_series(self):
        raise NotImplementedError

    def __call__(self, x):
        raise NotImplementedError


# Over prime finite fields

class HypergeometricAlgebraic_GFp(HypergeometricAlgebraic):
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        # TODO: do we want to simplify automatically if the
        # hypergeometric series is a polynomial?
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar)
        self._p = p = self.base_ring().cardinality()
        self._coeffs = [Qp(p, 1)(self._scalar)]

    def power_series(self, prec):
        S = self.parent().power_series_ring()
        self._compute_coeffs(prec)
        try:
            f = S(self._coeffs, prec=prec)
        except ValueError:
            raise ValueError("denominator appears in the series at the required precision")
        return f

    def is_almost_defined(self):
        p = self._char
        d = self.denominator()
        if d.gcd(p) > 1:
            return False
        u = 1
        if not self._parameters.parenthesis_criterion(u):
            return False
        u = p % d
        while u != 1:
            if not self._parameters.parenthesis_criterion(u):
                return False
            u = p*u % d
        return True

    def is_defined(self):
        p = self._char
        if not self.is_almost_defined():
            return False
        bound = self._parameters.bound
        if bound < p:
            return True
        prec = 1 + p ** ceil(log(self._parameters.bound, p))
        try:
            self.series(prec)
        except ValueError:
            return False
        return True

    def is_defined_conjectural(self):
        p = self._char
        if not self.is_almost_defined():
            return False
        bound = self._parameters.bound
        q = p
        while q <= bound:
            if not self._parameters.q_parenthesis_criterion(q):
                return False
            q *= p
        return True

    def __call__(self, x):
        return self.polynomial()(x)

    def is_polynomial(self):
        raise NotImplementedError

    def degree(self):
        raise NotImplementedError

    def polynomial(self):
        raise NotImplementedError

    def is_algebraic(self):
        return True

    def p_curvature(self):
        L = self.differential_operator()
        K = L.base_ring().fraction_field()
        S = OrePolynomialRing(K, L.parent().twisting_derivation().extend_to_fraction_field(), names='d')
        L = S(L.list())
        d = S.gen()
        p = self._char
        rows = [ ]
        n = L.degree()
        for i in range(p, p + n):
            Li = d**i % L
            rows.append([Li[j] for j in range(n)])
        return matrix(rows)

    def p_curvature_corank(self):  # maybe p_curvature_rank is preferable?
        # TODO: check if it is also correct when the parameters are not balanced
        return self._parameters.q_interlacing_number(self._char)

    def dwork_relation(self):
        r"""
        Return (P1, h1), ..., (Ps, hs) such that

            self = P1*h1^p + ... + Ps*hs^p
        """
        parameters = self._parameters
        if not parameters.is_balanced():
            raise ValueError("the hypergeometric function must be a pFq with q = p-1")
        p = self._char
        H = self.parent()
        F = H.base_ring()
        Hp = H.change_ring(Qp(p, 1))
        x = H.polynomial_ring().gen()
        coeffs = self._coeffs
        Ps = {}
        for r in range(p):
            params = parameters.shift(r).dwork_image(p)
            _, s = Hp(params)._val_pos()
            h = H(params.shift(s))
            e = s*p + r
            if e >= len(coeffs):
                self._compute_coeffs(e + 1)
            c = F(coeffs[e])
            if c:
                if h in Ps:
                    Ps[h] += c * x**e
                else:
                    Ps[h] = c * x**e
        return Ps

    def annihilating_ore_polynomial(self, var='Frob'):
        # QUESTION: does this method actually return the
        # minimal Ore polynomial annihilating self?
        # Probably not :-(
        parameters = self._parameters
        if not parameters.is_balanced():
            raise NotImplementedError("the hypergeometric function is not a pFq with q = p-1")

        p = self._char
        S = self.parent().polynomial_ring()
        zero = S.zero()
        Frob = S.frobenius_endomorphism()
        Ore = OrePolynomialRing(S, Frob, names=var)

        # We remove the scalar
        if self._scalar == 0:
            return Ore.one()
        self = self.parent()(parameters)

        order = parameters.frobenius_order(p)
        bound = self.p_curvature_corank()

        rows = [{self: S.one()}]
        # If row is the i-th item of rows, we have:
        #   self = sum_g row[g] * g**(p**i)
        q = 1
        while True:
            row = {}
            previous_row = rows[-1]
            for _ in range(order):
                row = {}
                for g, P in previous_row.items():
                    for h, Q in g.dwork_relation().items():
                        # here g = sum(Q * h^p)
                        if h in row:
                            row[h] += P * insert_zeroes(Q, q)
                        else:
                            row[h] = P * insert_zeroes(Q, q)
                previous_row = row
                q *= p  # q = p**i
            rows.append(row)

            i = len(rows)
            Mrows = []
            Mqo = 1
            columns = {}
            for j in range(i-1, max(-1, i-2-bound), -1):
                for col in rows[j]:
                    columns[col] = None
            for j in range(i-1, max(-1, i-2-bound), -1):
                Mrow = []
                for col in columns:
                    Mrow.append(insert_zeroes(rows[j].get(col, zero), Mqo))
                Mrows.append(Mrow)
                Mqo *= p ** order
            M = matrix(S, Mrows)

            ker = kernel(M)
            if ker is not None:
                return insert_zeroes(Ore(ker), order)

    def is_lucas(self):
        p = self._char
        if self._parameters.frobenius_order(p) > 1:
            # TODO: check this
            return False
        S = self.parent().polynomial_ring()
        K = S.fraction_field()
        Ore = OrePolynomialRing(K, K.frobenius_endomorphism(), names='F')
        Z = Ore(self.annihilating_ore_polynomial())
        Ap = self.series(p).polynomial()
        F = Ap * Ore.gen() - 1
        return (Z % F).is_zero()


# Parent
########

class HypergeometricToSR(Map):
    def _call_(self, h):
        return h.scalar() * hypergeometric(h.top(), h.bottom(), SR.var(h.parent().variable_name()))


class ScalarMultiplication(Action):
    def _act_(self, scalar, h):
        return h.parent()(h, scalar=scalar)


class HypergeometricFunctions(Parent, UniqueRepresentation):
    def __init__(self, base, name, category=None):
        self._name = normalize_names(1, name)[0]
        self._latex_name = latex_variable_name(self._name)
        self._char = char = base.characteristic()
        if char == 0:
            base = pushout(base, QQ)
        if base in FiniteFields() and base.is_prime_field():
            self.Element = HypergeometricAlgebraic_GFp
        elif base is QQ:
            self.Element = HypergeometricAlgebraic_QQ
        elif isinstance(base, pAdicGeneric):
            self.Element = HypergeometricAlgebraic_padic
        else:
            self.Element = HypergeometricAlgebraic
        Parent.__init__(self, base, category=category)
        self.register_action(ScalarMultiplication(base, self, False, operator.mul))
        self.register_action(ScalarMultiplication(base, self, True, operator.mul))
        if char == 0:
            SR.register_coercion(HypergeometricToSR(self.Hom(SR)))

    def _repr_(self):
        return "Hypergeometric functions in %s over %s" % (self._name, self._base)

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def _coerce_map_from_(self, other):
        if (isinstance(other, HypergeometricFunctions)
        and other.has_coerce_map_from(self)):
            return True

    def _pushout_(self, other):
        if isinstance(other, HypergeometricFunctions) and self._name == other._name:
            base = pushout(self.base_ring(), other.base_ring())
            if base is not None:
                return HypergeometricFunctions(base, self._name)
        if SR.has_coerce_map_from(other):
            return SR

    def base_ring(self):
        return self._base

    def variable_name(self):
        return self._name

    def latex_variable_name(self):
        return self._latex_name

    def change_ring(self, R):
        return HypergeometricFunctions(R, self._name)

    def change_variable_name(self, name):
        return HypergeometricFunctions(self._base, name)

    def polynomial_ring(self):
        return PolynomialRing(self.base_ring(), self._name)

    def power_series_ring(self, default_prec=None):
        return PowerSeriesRing(self.base_ring(), self._name, default_prec=default_prec)
