
from sage.misc.lazy_import import lazy_import

lazy_import('sage.functions.piecewise', 'piecewise')
lazy_import('sage.functions.error', ['erf', 'erfc', 'erfi', 'erfinv',
                                     'fresnel_sin', 'fresnel_cos'])

from sage.functions.trig import (sin, cos, sec, csc, cot, tan,
                                 asin, acos, atan,
                                 acot, acsc, asec,
                                 arcsin, arccos, arctan,
                                 arccot, arccsc, arcsec,
                                 arctan2, atan2)

from sage.functions.hyperbolic import (tanh, sinh, cosh, coth, sech, csch,
                                       asinh, acosh, atanh, acoth, asech, acsch,
                                       arcsinh, arccosh, arctanh, arccoth, arcsech, arccsch)

reciprocal_trig_functions = {'sec': cos, 'csc': sin,
                             'cot': tan, 'sech': cosh, 'csch': sinh, 'coth': tanh}


from sage.functions.other import (ceil, floor, abs_symbolic, sqrt, real_nth_root,
                                  arg, real_part, real, frac,
                                  factorial, binomial,
                                  imag_part, imag, imaginary, conjugate, cases,
                                  complex_root_of)

from sage.functions.log import (exp, exp_polar, log, ln,
                                polylog, dilog, lambert_w, harmonic_number)

from sage.functions.transcendental import (zeta, zetaderiv, zeta_symmetric, hurwitz_zeta,
                                           dickman_rho, stieltjes)

from sage.functions.bessel import (bessel_I, bessel_J, bessel_K, bessel_Y,
                                   Bessel, struve_H, struve_L, hankel1, hankel2,
                                   spherical_bessel_J, spherical_bessel_Y,
                                   spherical_hankel1, spherical_hankel2)

from sage.functions.special import (spherical_harmonic, elliptic_e,
                                    elliptic_f, elliptic_ec, elliptic_eu,
                                    elliptic_kc, elliptic_pi, elliptic_j)

from sage.functions.jacobi import (jacobi, inverse_jacobi, jacobi_nd, jacobi_ns, jacobi_nc,
                                   jacobi_dn, jacobi_ds, jacobi_dc, jacobi_sn, jacobi_sd,
                                   jacobi_sc, jacobi_cn, jacobi_cd, jacobi_cs, jacobi_am,
                                   inverse_jacobi_nd, inverse_jacobi_ns, inverse_jacobi_nc,
                                   inverse_jacobi_dn, inverse_jacobi_ds, inverse_jacobi_dc,
                                   inverse_jacobi_sn, inverse_jacobi_sd, inverse_jacobi_sc,
                                   inverse_jacobi_cn, inverse_jacobi_cd, inverse_jacobi_cs)

from sage.functions.orthogonal_polys import (chebyshev_T,
                                             chebyshev_U,
                                             gen_laguerre,
                                             gen_legendre_P,
                                             gen_legendre_Q,
                                             hermite,
                                             jacobi_P,
                                             laguerre,
                                             legendre_P,
                                             legendre_Q,
                                             ultraspherical,
                                             gegenbauer,
                                             krawtchouk,
                                             meixner,
                                             hahn)

from sage.functions.spike_function import spike_function

from sage.functions.prime_pi import legendre_phi, partial_sieve_function, prime_pi

from sage.functions.wigner import (wigner_3j, clebsch_gordan, racah, wigner_6j,
                                   wigner_9j, gaunt)

from sage.functions.generalized import (dirac_delta, heaviside, unit_step, sgn, sign,
                                        kronecker_delta)

from sage.functions.min_max import max_symbolic, min_symbolic

from sage.functions.airy import airy_ai, airy_ai_prime, airy_bi, airy_bi_prime

from sage.functions.exp_integral import (exp_integral_e, exp_integral_e1, log_integral, li, Li,
                                         log_integral_offset,
                                         sin_integral, cos_integral, Si, Ci,
                                         sinh_integral, cosh_integral, Shi, Chi,
                                         exponential_integral_1, Ei, exp_integral_ei)

from sage.functions.hypergeometric import hypergeometric, hypergeometric_M, hypergeometric_U

from sage.functions.gamma import (gamma, psi, beta, log_gamma,
                                  gamma_inc, gamma_inc_lower)

Γ = gamma
ψ = psi
ζ = zeta
del lazy_import
