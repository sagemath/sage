"""
Combinatorial Logarithm

This file provides the cycle index series for the virtual species `\Omega`,
the 'combinatorial logarithm', defined to be the compositional inverse of
the species `E^{+}` of nonempty sets:

.. MATH::

    \Omega \circ E^{+} = E^{+} \circ \Omega = X.

.. warning::

    This module is now deprecated.  Please use
    :meth:`sage.combinat.species.generating_series.CycleIndexSeriesRing.exponential`
    instead of :func:`CombinatorialLogarithmSeries`.

AUTHORS:

- Andrew Gainer-Dewar (2013): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Andrew Gainer-Dewar <andrew.gainer.dewar@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.species.generating_series import CycleIndexSeriesRing, LogarithmCycleIndexSeries
from sage.rings.all import QQ
from sage.misc.cachefunc import cached_function
from sage.misc.superseded import deprecation

@cached_function
def CombinatorialLogarithmSeries(R=QQ):
    r"""
    Return the cycle index series of the virtual species `\Omega`, the compositional inverse
    of the species `E^{+}` of nonempty sets.

    The notion of virtual species is treated thoroughly in [BLL]_. The specific algorithm used
    here to compute the cycle index of `\Omega` is found in [Labelle]_.

    EXAMPLES:

    The virtual species `\Omega` is 'properly virtual', in the sense that its cycle index
    has negative coefficients::

        sage: from sage.combinat.species.combinatorial_logarithm import CombinatorialLogarithmSeries
        sage: CombinatorialLogarithmSeries().coefficients(4)
        doctest:...: DeprecationWarning: CombinatorialLogarithmSeries is deprecated, use CycleIndexSeriesRing(R).logarithm_series() or CycleIndexSeries().logarithm() instead
        See http://trac.sagemath.org/14846 for details.
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]

    Its defining property is that `\Omega \circ E^{+} = E^{+} \circ \Omega = X` (that is, that
    composition with `E^{+}` in both directions yields the multiplicative identity `X`)::

        sage: Eplus = sage.combinat.species.set_species.SetSpecies(min=1).cycle_index_series()
        sage: CombinatorialLogarithmSeries().compose(Eplus).coefficients(4)
        [0, p[1], 0, 0]
    """
    deprecation(14846, "CombinatorialLogarithmSeries is deprecated, use CycleIndexSeriesRing(R).logarithm_series() or CycleIndexSeries().logarithm() instead")
    return LogarithmCycleIndexSeries(R)
