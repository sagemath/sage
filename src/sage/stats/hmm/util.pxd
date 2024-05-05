from sage.stats.time_series cimport TimeSeries

cdef class HMM_Util:
    cpdef normalize_probability_TimeSeries(self, TimeSeries T, Py_ssize_t i, Py_ssize_t j)
    cpdef TimeSeries initial_probs_to_TimeSeries(self, pi, bint normalize)
    cpdef TimeSeries state_matrix_to_TimeSeries(self, A, int N, bint normalize)

