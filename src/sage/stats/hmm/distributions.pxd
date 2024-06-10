#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

from sage.stats.time_series cimport TimeSeries
from sage.stats.intlist cimport IntList
from sage.misc.randstate cimport randstate

cdef class Distribution:
    pass

cdef class DiscreteDistribution(Distribution):
    cdef object v

cdef class GaussianDistribution(Distribution):
    pass

cdef class GaussianMixtureDistribution(Distribution):
    cdef TimeSeries c0, c1, param
    cdef IntList fixed

    cdef double _sample(self, randstate rstate) noexcept
    cpdef double prob(self, double x) noexcept
    cpdef double prob_m(self, double x, int m) noexcept
    cpdef is_fixed(self, i=?)
