from sage.rings.integer cimport Integer
from sage.rings.real_mpfr cimport RealNumber
from sage.stats.distributions.dgs cimport dgs_disc_gauss_mp_t, dgs_disc_gauss_dp_t
from sage.structure.sage_object cimport SageObject


cdef class DiscreteGaussianDistributionIntegerSampler(SageObject):
    cdef readonly RealNumber sigma
    cdef readonly RealNumber c
    cdef readonly Integer tau
    cdef readonly object algorithm
    cdef dgs_disc_gauss_mp_t *_gen_mp
    cdef dgs_disc_gauss_dp_t *_gen_dp
