# sage_setup: distribution = sagemath-graphs
# Cached informations about OA(k,n)
#
# - max_true: max k such that OA(k,n,existence=True) previously returned True
#
# - min_unknown: min k such that OA(k,n,existence=True) previously returned
#                Unknown
# - max_unknown: max k such that OA(k,n,existence=True) previously returned
#                Unknown
# - min_false: min k such that OA(k,n,existence=True) previously returned False

cdef struct cache_entry:
    unsigned short max_true
    unsigned short min_unknown
    unsigned short max_unknown
    unsigned short min_false

cdef cache_entry * _OA_cache
cdef int _OA_cache_size

cpdef _OA_cache_get(int k, int n)
