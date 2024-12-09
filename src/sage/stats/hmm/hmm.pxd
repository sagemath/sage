#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL) v2+.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


from sage.stats.time_series cimport TimeSeries
from sage.stats.intlist cimport IntList

cdef class HiddenMarkovModel:
    cdef int N
    cdef TimeSeries A, pi

    cdef TimeSeries _baum_welch_gamma(self, TimeSeries alpha, TimeSeries beta)
