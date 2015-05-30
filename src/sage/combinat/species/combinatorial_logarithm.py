"""
Combinatorial Logarithm

This file provides the cycle index series for the virtual species `\Omega`,
the 'combinatorial logarithm', defined to be the compositional inverse of
the species `E^{+}` of nonempty sets:

.. MATH::
    \Omega \circ E^{+} = E^{+} \circ \Omega = X.

AUTHORS:

- Andrew Gainer-Dewar (2013): initial version

TESTS::

    sage: from sage.combinat.species.combinatorial_logarithm import CombinatorialLogarithmSeries
    sage: CombinatorialLogarithmSeries().coefficients(5)
    [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3], -1/4*p[1, 1, 1, 1] + 1/4*p[2, 2]]

    sage: Eplus = sage.combinat.species.set_species.SetSpecies(min=1).cycle_index_series()
    sage: CombinatorialLogarithmSeries().compose(Eplus).coefficients(4)
    [0, p[1], 0, 0]

"""
#*****************************************************************************
#       Copyright (C) 2013 Andrew Gainer-Dewar <andrew.gainer.dewar@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_function
from sage.combinat.species.stream import _integers_from
from sage.combinat.sf.all import SymmetricFunctions
from sage.combinat.species.generating_series import CycleIndexSeriesRing
from sage.rings.all import RationalField, Integer, divisors

@cached_function
def _cl_term(n, R = RationalField()):
    """
    Compute the order-n term of the cycle index series of the virtual species `\Omega`,
    the compositional inverse of the species `E^{+}` of nonempty sets.

    EXAMPLES::

        sage: from sage.combinat.species.combinatorial_logarithm import _cl_term
        sage: [_cl_term(i) for i in range(4)]
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]
    """

    n = Integer(n) #check that n is an integer

    p = SymmetricFunctions(R).power()

    res = p.zero()
    if n == 1:
        res = p([1])
    elif n > 1:
        res = 1/n * ((-1)**(n-1) * p([1])**n - sum(d * p([Integer(n/d)]).plethysm(_cl_term(d, R)) for d in divisors(n)[:-1]))

    return res

def _cl_gen (R = RationalField()):
    """
    Produce a generator which yields the terms of the cycle index series of the virtual species
    `\Omega`, the compositional inverse of the species `E^{+}` of nonempty sets.

    EXAMPLES::

        sage: from sage.combinat.species.combinatorial_logarithm import _cl_gen
        sage: g = _cl_gen()
        sage: [next(g) for i in range(4)]
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]
    """
    return (_cl_term(i, R) for i in _integers_from(0))

@cached_function
def CombinatorialLogarithmSeries(R = RationalField()):
    """
    Return the cycle index series of the virtual species `\Omega`, the compositional inverse
    of the species `E^{+}` of nonempty sets.

    The notion of virtual species is treated thoroughly in [BLL]_. The specific algorithm used
    here to compute the cycle index of `\Omega` is found in [Labelle]_.

    EXAMPLES:

    The virtual species `\Omega` is 'properly virtual', in the sense that its cycle index
    has negative coefficients::

        sage: from sage.combinat.species.combinatorial_logarithm import CombinatorialLogarithmSeries
        sage: CombinatorialLogarithmSeries().coefficients(4)
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]

    Its defining property is that `\Omega \circ E^{+} = E^{+} \circ \Omega = X` (that is, that
    composition with `E^{+}` in both directions yields the multiplicative identity `X`)::

        sage: Eplus = sage.combinat.species.set_species.SetSpecies(min=1).cycle_index_series()
        sage: CombinatorialLogarithmSeries().compose(Eplus).coefficients(4)
        [0, p[1], 0, 0]

    REFERENCES:

    .. [BLL] F. Bergeron, G. Labelle, and P. Leroux. "Combinatorial species and tree-like structures". Encyclopedia of Mathematics and its Applications, vol. 67, Cambridge Univ. Press. 1998.

    .. [Labelle] G. Labelle. "New combinatorial computational methods arising from pseudo-singletons." DMTCS Proceedings 1, 2008.
    """
    CIS = CycleIndexSeriesRing(R)
    return CIS(_cl_gen(R))
