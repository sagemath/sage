# distutils: extra_compile_args = -D_XPG6
"""
Flint imports

TESTS:

Import this module::

    sage: import sage.libs.flint.flint_sage

We verify that :trac:`6919` is correctly fixed::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: A = 2^(2^17+2^15)
    sage: a = A * x^31
    sage: b = (A * x) * x^30
    sage: a == b
    True
"""

# cimport all .pxd files to make sure they compile
from .fmpz_poly_mat cimport *
from .fmpq_mpoly cimport *
from .nmod_poly_mat cimport *
from .perm cimport *
from .nmod_vec cimport *
from .fmpzi cimport *
from .mpoly cimport *
from .mpfr_mat cimport *
from .nmod_poly cimport *
from .long_extras cimport *
from .partitions cimport *
from .fq_nmod_poly_factor cimport *
from .mpf_mat cimport *
from .fmpz_poly_q cimport *
from .ca_vec cimport *
from .double_extras cimport *
from .fq_nmod_mpoly_factor cimport *
from .fmpz_mpoly cimport *
from .fmpq_vec cimport *
from .fq_poly cimport *
from .fq_nmod_mat cimport *
from .calcium cimport *
from .fmpz_vec cimport *
from .fexpr_builtin cimport *
from .ca_poly cimport *
from .fq_embed cimport *
from .nmod_poly_factor cimport *
from .fmpz_extras cimport *
from .acb_mat cimport *
from .gr_generic cimport *
from .gr_poly cimport *
from .fq_zech_poly_factor cimport *
from .fmpz_mpoly_q cimport *
from .nmod cimport *
from .ca_mat cimport *
from .fq_zech_vec cimport *
from .arb cimport *
from .nf_elem cimport *
from .fmpz_poly cimport *
from .fq_nmod_embed cimport *
from .arb_poly cimport *
from .fmpq cimport *
from .fmpq_mpoly_factor cimport *
from .acb_calc cimport *
from .fq_vec cimport *
from .padic_mat cimport *
from .fmpz cimport *
from .fmpz_mat cimport *
from .fmpz_mod_mat cimport *
from .fq_zech cimport *
from .double_interval cimport *
from .fq_default_mat cimport *
from .fmpz_mod_mpoly cimport *
from .qsieve cimport *
from .qfb cimport *
from .thread_pool cimport *
from .fmpz_mod_mpoly_factor cimport *
from .acb_modular cimport *
from .fq_nmod_poly cimport *
from .profiler cimport *
from .acb_hypgeom cimport *
from .d_mat cimport *
from .fq_zech_poly cimport *
from .fmpz_mod cimport *
from .qqbar cimport *
from .hypgeom cimport *
from .acb_elliptic cimport *
from .acb_dft cimport *
from .d_vec cimport *
from .ulong_extras cimport *
from .dlog cimport *
from .bool_mat cimport *
from .fmpq_poly cimport *
from .fexpr cimport *
from .mag cimport *
from .fmpz_mpoly_factor cimport *
from .mpfr_vec cimport *
from .ca_ext cimport *
from .gr cimport *
from .acb_poly cimport *
from .fft cimport *
from .padic cimport *
from .dirichlet cimport *
from .fmpz_mod_poly_factor cimport *
from .fq_default cimport *
from .fq cimport *
from .gr_vec cimport *
from .fq_nmod cimport *
from .fq_zech_embed cimport *
from .nmod_mpoly_factor cimport *
from .flint cimport *
from .fmpz_factor cimport *
from .qadic cimport *
from .nf cimport *
from .gr_mpoly cimport *
from .arb_fpwrap cimport *
from .fq_default_poly cimport *
from .aprcl cimport *
from .acf cimport *
from .gr_special cimport *
from .arf cimport *
from .fq_zech_mat cimport *
from .gr_mat cimport *
from .arb_fmpz_poly cimport *
from .fq_nmod_vec cimport *
from .fmpz_mod_poly cimport *
from .bernoulli cimport *
from .fq_default_poly_factor cimport *
from .arb_mat cimport *
from .arb_hypgeom cimport *
from .fq_poly_factor cimport *
from .nmod_mpoly cimport *
from .acb cimport *
from .fmpz_lll cimport *
from .fmpz_poly_factor cimport *
from .ca_field cimport *
from .acb_dirichlet cimport *
from .arith cimport *
from .fmpz_mod_vec cimport *
from .fmpq_mat cimport *
from .fq_mat cimport *
from .mpn_extras cimport *
from .mpf_vec cimport *
from .ca cimport *
from .padic_poly cimport *
from .nmod_mat cimport *
from .fq_nmod_mpoly cimport *
from .arb_calc cimport *

# Try to clean up after ourselves before sage terminates. This
# probably doesn't do anything if your copy of flint is re-entrant
# (and most are). Moreover it isn't strictly necessary, because the OS
# will reclaim these resources anyway after sage terminates. However
# this might reveal other bugs, and can help tools like valgrind do
# their jobs.
import atexit
atexit.register(_fmpz_cleanup_mpz_content)
