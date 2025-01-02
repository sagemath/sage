# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Eisenstein series
"""
# ****************************************************************************
#       Copyright (C) 2004-2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.functions import lcm
from sage.arith.misc import bernoulli, divisors, is_squarefree
from sage.misc.lazy_import import lazy_import
from sage.misc.timing import cputime
from sage.modular.arithgroup.congroup_gammaH import GammaH_class
from sage.modular.dirichlet import DirichletGroup
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ

lazy_import('sage.modular.modform.eis_series_cython', ['eisenstein_series_poly', 'Ek_ZZ'])
lazy_import('sage.rings.number_field.number_field', 'CyclotomicField')


def eisenstein_series_qexp(k, prec=10, K=QQ, var='q', normalization='linear'):
    r"""
    Return the `q`-expansion of the normalized weight `k` Eisenstein series on
    `\SL_2(\ZZ)` to precision ``prec`` in the ring `K`. Three normalizations
    are available, depending on the parameter ``normalization``; the default
    normalization is the one for which the linear coefficient is 1.

    INPUT:

    - ``k`` -- an even positive integer

    - ``prec`` -- (default: 10) a nonnegative integer

    - ``K`` -- (default: `\QQ`) a ring

    - ``var`` -- (default: ``'q'``) variable name to use for `q`-expansion

    - ``normalization`` -- (default: ``'linear'``) normalization to use. If this
      is ``'linear'``, then the series will be normalized so that the linear
      term is 1. If it is ``'constant'``, the series will be normalized to have
      constant term 1. If it is ``'integral'``, then the series will be
      normalized to have integer coefficients and no common factor, and linear
      term that is positive. Note that ``'integral'`` will work over arbitrary
      base rings, while ``'linear'`` or ``'constant'`` will fail if the
      denominator (resp. numerator) of `B_k / 2k` is invertible.

    ALGORITHM:

        We know `E_k = \text{constant} + \sum_n \sigma_{k-1}(n) q^n`. So we
        compute all the `\sigma_{k-1}(n)` simultaneously, using the fact that
        `\sigma` is multiplicative.

    EXAMPLES::

        sage: eisenstein_series_qexp(2,5)
        -1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + O(q^5)
        sage: eisenstein_series_qexp(2,0)
        O(q^0)
        sage: eisenstein_series_qexp(2,5,GF(7))
        2 + q + 3*q^2 + 4*q^3 + O(q^5)
        sage: eisenstein_series_qexp(2,5,GF(7),var='T')
        2 + T + 3*T^2 + 4*T^3 + O(T^5)

    We illustrate the use of the ``normalization`` parameter::

        sage: eisenstein_series_qexp(12, 5, normalization='integral')
        691 + 65520*q + 134250480*q^2 + 11606736960*q^3 + 274945048560*q^4 + O(q^5)
        sage: eisenstein_series_qexp(12, 5, normalization='constant')
        1 + 65520/691*q + 134250480/691*q^2 + 11606736960/691*q^3 + 274945048560/691*q^4 + O(q^5)
        sage: eisenstein_series_qexp(12, 5, normalization='linear')
        691/65520 + q + 2049*q^2 + 177148*q^3 + 4196353*q^4 + O(q^5)
        sage: eisenstein_series_qexp(12, 50, K=GF(13), normalization='constant')
        1 + O(q^50)

    TESTS:

    Test that :issue:`5102` is fixed::

        sage: eisenstein_series_qexp(10, 30, GF(17))
        15 + q + 3*q^2 + 15*q^3 + 7*q^4 + 13*q^5 + 11*q^6 + 11*q^7 + 15*q^8 + 7*q^9 + 5*q^10 + 7*q^11 + 3*q^12 + 14*q^13 + 16*q^14 + 8*q^15 + 14*q^16 + q^17 + 4*q^18 + 3*q^19 + 6*q^20 + 12*q^21 + 4*q^22 + 12*q^23 + 4*q^24 + 4*q^25 + 8*q^26 + 14*q^27 + 9*q^28 + 6*q^29 + O(q^30)

    This shows that the bug reported at :issue:`8291` is fixed::

        sage: eisenstein_series_qexp(26, 10, GF(13))
        7 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + 12*q^6 + 8*q^7 + 2*q^8 + O(q^10)

    We check that the function behaves properly over finite-characteristic base rings::

        sage: eisenstein_series_qexp(12, 5, K = Zmod(691), normalization='integral')
        566*q + 236*q^2 + 286*q^3 + 194*q^4 + O(q^5)
        sage: eisenstein_series_qexp(12, 5, K = Zmod(691), normalization='constant')
        Traceback (most recent call last):
        ...
        ValueError: The numerator of -B_k/(2*k) (=691) must be invertible in the ring Ring of integers modulo 691
        sage: eisenstein_series_qexp(12, 5, K = Zmod(691), normalization='linear')
        q + 667*q^2 + 252*q^3 + 601*q^4 + O(q^5)

        sage: eisenstein_series_qexp(12, 5, K = Zmod(2), normalization='integral')
        1 + O(q^5)
        sage: eisenstein_series_qexp(12, 5, K = Zmod(2), normalization='constant')
        1 + O(q^5)
        sage: eisenstein_series_qexp(12, 5, K = Zmod(2), normalization='linear')
        Traceback (most recent call last):
        ...
        ValueError: The denominator of -B_k/(2*k) (=65520) must be invertible in the ring Ring of integers modulo 2

    AUTHORS:

    - William Stein: original implementation

    - Craig Citro (2007-06-01): rewrote for massive speedup

    - Martin Raum (2009-08-02): port to cython for speedup

    - David Loeffler (2010-04-07): work around an integer overflow when `k` is large

    - David Loeffler (2012-03-15): add options for alternative normalizations
      (motivated by :issue:`12043`)
    """
    # we use this to prevent computation if it would fail anyway.
    if k <= 0 or k % 2 == 1:
        raise ValueError("k must be positive and even")

    a0 = - bernoulli(k) / (2*k)

    if normalization == 'linear':
        a0den = a0.denominator()
        try:
            a0fac = K(1/a0den)
        except ZeroDivisionError:
            raise ValueError("The denominator of -B_k/(2*k) (=%s) must be invertible in the ring %s" % (a0den, K))
    elif normalization == 'constant':
        a0num = a0.numerator()
        try:
            a0fac = K(1/a0num)
        except ZeroDivisionError:
            raise ValueError("The numerator of -B_k/(2*k) (=%s) must be invertible in the ring %s" % (a0num, K))
    elif normalization == 'integral':
        a0fac = None
    else:
        raise ValueError("Normalization (=%s) must be one of 'linear', 'constant', 'integral'" % normalization)

    R = PowerSeriesRing(K, var)
    if K == QQ and normalization == 'linear':
        ls = Ek_ZZ(k, prec)
        # The following is *dramatically* faster than doing the more natural
        # "R(ls)" would be:
        E = ZZ[var](ls, prec=prec, check=False).change_ring(QQ)
        if len(ls) > 0:
            E._unsafe_mutate(0, a0)
        return R(E, prec)
        # The following is an older slower alternative to the above three lines:
        # return a0fac*R(eisenstein_series_poly(k, prec).list(), prec=prec, check=False)
    else:
        # This used to work with check=False, but that can only be regarded as
        # an improbable lucky miracle. Enabling checking is a noticeable speed
        # regression; the morally right fix would be to expose FLINT's
        # fmpz_poly_to_nmod_poly command (at least for word-sized N).
        if a0fac is not None:
            return a0fac*R(eisenstein_series_poly(k, prec).list(), prec=prec, check=True)
        else:
            return R(eisenstein_series_poly(k, prec).list(), prec=prec, check=True)


def __common_minimal_basering(chi, psi):
    """
    Find the smallest basering over which chi and psi are valued, and
    return new chi and psi valued in that ring.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(1)[0], DirichletGroup(1)[0])
        (Dirichlet character modulo 1 of conductor 1,
         Dirichlet character modulo 1 of conductor 1)
        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(3).0, DirichletGroup(5).0)
        (Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1,
         Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4)
        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(12).0, DirichletGroup(36).0)
        (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
         Dirichlet character modulo 36 of conductor 4 mapping 19 |--> -1, 29 |--> 1)
    """
    chi = chi.minimize_base_ring()
    psi = psi.minimize_base_ring()
    n = lcm(chi.base_ring().zeta().multiplicative_order(),
                  psi.base_ring().zeta().multiplicative_order())
    if n <= 2:
        K = QQ
    else:
        K = CyclotomicField(n)
    chi = chi.change_ring(K)
    psi = psi.change_ring(K)
    return chi, psi


def __find_eisen_chars(character, k):
    r"""
    Find all triples `(\psi_1, \psi_2, t)` that give rise to an Eisenstein series of the given weight and character.

    EXAMPLES::

        sage: sage.modular.modform.eis_series.__find_eisen_chars(DirichletGroup(36).0, 4)
        []

        sage: pars =  sage.modular.modform.eis_series.__find_eisen_chars(DirichletGroup(36).0, 5)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, 1), 1),
        ((1, 1), (-1, 1), 3),
        ((1, 1), (-1, 1), 9),
        ((1, -1), (-1, -1), 1),
        ((-1, 1), (1, 1), 1),
        ((-1, 1), (1, 1), 3),
        ((-1, 1), (1, 1), 9),
        ((-1, -1), (1, -1), 1)]
    """
    from sage.misc.verbose import verbose

    N = character.modulus()
    if character.is_trivial():
        if k % 2:
            return []
        char_inv = ~character
        V = [(character, char_inv, t) for t in divisors(N) if t > 1]
        if k != 2:
            V.insert(0, (character, char_inv, 1))
        if is_squarefree(N):
            return V
        # Now include all pairs (chi,chi^(-1)) such that cond(chi)^2 divides N:
        # TODO: Optimize -- this is presumably way too hard work below.
        G = DirichletGroup(N)
        for chi in G:
            if not chi.is_trivial():
                f = chi.conductor()
                if N % (f**2) == 0:
                    chi = chi.minimize_base_ring()
                    chi_inv = ~chi
                    for t in divisors(N//(f**2)):
                        V.insert(0, (chi, chi_inv, t))
        return V

    eps = character
    if eps(-1) != (-1)**k:
        return []
    eps = eps.maximize_base_ring()
    G = eps.parent()

    # Find all pairs chi, psi such that:
    #
    #  (1) cond(chi)*cond(psi) divides the level, and
    #
    #  (2) chi*psi == eps, where eps is the nebentypus character of self.
    #
    # See [Miyake, Modular Forms] Lemma 7.1.1.

    C = {}

    t0 = cputime()

    for e in G:
        m = Integer(e.conductor())
        if m in C:
            C[m].append(e)
        else:
            C[m] = [e]

    verbose("Enumeration with conductors.", t0)

    params = []
    for L in divisors(N):
        verbose("divisor %s" % L)
        if L not in C:
            continue
        GL = C[L]
        for R in divisors(N // L):
            if R not in C:
                continue
            GR = C[R]
            for chi in GL:
                for psi in GR:
                    if chi*psi == eps:
                        chi0, psi0 = __common_minimal_basering(chi, psi)
                        for t in divisors(N//(R*L)):
                            if k != 1 or (psi0, chi0, t) not in params:
                                params.append((chi0, psi0, t))
    return params


def __find_eisen_chars_gammaH(N, H, k):
    r"""
    Find all triples `(\psi_1, \psi_2, t)` that give rise to an Eisenstein series of weight `k` on
    `\Gamma_H(N)`.

    EXAMPLES::

        sage: pars =  sage.modular.modform.eis_series.__find_eisen_chars_gammaH(15, [2], 5)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, -1), 1), ((-1, 1), (1, -1), 1), ((1, -1), (-1, 1), 1), ((-1, -1), (1, 1), 1)]
    """
    params = []
    for chi in DirichletGroup(N):
        if all(chi(h) == 1 for h in H):
            params += __find_eisen_chars(chi, k)
    return params


def __find_eisen_chars_gamma1(N, k):
    r"""
    Find all triples `(\psi_1, \psi_2, t)` that give rise to an Eisenstein series of weight `k` on
    `\Gamma_1(N)`.

    EXAMPLES::

        sage: pars = sage.modular.modform.eis_series.__find_eisen_chars_gamma1(12, 4)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (1, 1), 1),
        ((1, 1), (1, 1), 2),
        ((1, 1), (1, 1), 3),
        ((1, 1), (1, 1), 4),
        ((1, 1), (1, 1), 6),
        ((1, 1), (1, 1), 12),
        ((1, 1), (-1, -1), 1),
        ((-1, -1), (1, 1), 1),
        ((-1, 1), (1, -1), 1),
        ((1, -1), (-1, 1), 1)]

        sage: pars =  sage.modular.modform.eis_series.__find_eisen_chars_gamma1(12, 5)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, 1), 1),
        ((1, 1), (-1, 1), 3),
        ((-1, 1), (1, 1), 1),
        ((-1, 1), (1, 1), 3),
        ((1, 1), (1, -1), 1),
        ((1, 1), (1, -1), 2),
        ((1, 1), (1, -1), 4),
        ((1, -1), (1, 1), 1),
        ((1, -1), (1, 1), 2),
        ((1, -1), (1, 1), 4)]
    """
    pairs = []
    s = (-1)**k
    G = DirichletGroup(N)
    E = list(G)
    parity = [c(-1) for c in E]
    for i in range(len(E)):
        for j in range(i,len(E)):
            if parity[i]*parity[j] == s and N % (E[i].conductor()*E[j].conductor()) == 0:
                chi, psi = __common_minimal_basering(E[i], E[j])
                if k != 1:
                    pairs.append((chi, psi))
                    if i != j:
                        pairs.append((psi, chi))
                else:
                    # if weight is 1 then (chi, psi) and (chi, psi) are the
                    # same form
                    if psi.is_trivial() and not chi.is_trivial():
                        # need to put the trivial character first to get the L-value right
                        pairs.append((psi, chi))
                    else:
                        pairs.append((chi, psi))

    triples = []
    for chi, psi in pairs:
        c_chi = chi.conductor()
        c_psi = psi.conductor()
        D = divisors(N // (c_chi * c_psi))
        if k == 2 and chi.is_trivial() and psi.is_trivial():
            D.remove(1)
        chi, psi = __common_minimal_basering(chi, psi)
        triples.extend((chi, psi, t) for t in D)

    return triples


def eisenstein_series_lseries(weight, prec=53,
                              max_imaginary_part=0,
                              max_asymp_coeffs=40):
    r"""
    Return the `L`-series of the weight `2k` Eisenstein series `E_{2k}`
    on `\SL_2(\ZZ)`.

    This actually returns an interface to Tim Dokchitser's program
    for computing with the `L`-series of the Eisenstein series.
    See :class:`~sage.lfunctions.dokchitser.Dokchitser`.

    INPUT:

    - ``weight`` -- even integer

    - ``prec`` -- integer (bits precision)

    - ``max_imaginary_part`` -- real number

    - ``max_asymp_coeffs`` -- integer

    OUTPUT: the `L`-series of the Eisenstein series. This can be
    evaluated at argument `s`, or have
    :meth:`~sage.lfunctions.dokchitser.Dokchitser.derivative` called, etc.

    EXAMPLES:

    We compute with the `L`-series of `E_{16}` and then `E_{20}`::

        sage: L = eisenstein_series_lseries(16)
        sage: L(1)
        -0.291657724743874
        sage: L.derivative(1)
        0.0756072194360656
        sage: L = eisenstein_series_lseries(20)
        sage: L(2)
        -5.02355351645998

    Now with higher precision::

        sage: L = eisenstein_series_lseries(20, prec=200)
        sage: L(2)
        -5.0235535164599797471968418348135050804419155747868718371029
    """
    f = eisenstein_series_qexp(weight, prec)
    from sage.lfunctions.all import Dokchitser
    j = weight
    L = Dokchitser(conductor=1,
                   gammaV=[0, 1],
                   weight=j,
                   eps=(-1)**Integer(j // 2),
                   poles=[j],
                   # Using a string for residues is a hack but it works well
                   # since this will make PARI/GP compute sqrt(pi) with the
                   # right precision.
                   residues='[sqrt(Pi)*(%s)]' % ((-1)**Integer(j // 2) * bernoulli(j) / j),
                   prec=prec)

    s = 'coeff = %s;' % f.list()
    L.init_coeffs('coeff[k+1]', pari_precode=s,
                  max_imaginary_part=max_imaginary_part,
                  max_asymp_coeffs=max_asymp_coeffs)
    L.check_functional_equation()
    L.rename('L-series associated to the weight %s Eisenstein series %s on SL_2(Z)' % (j, f))
    return L


def compute_eisenstein_params(character, k):
    r"""
    Compute and return a list of all parameters `(\chi,\psi,t)` that
    define the Eisenstein series with given character and weight `k`.

    Only the parity of `k` is relevant (unless k = 1, which is a slightly different case).

    If ``character`` is an integer `N`, then the parameters for
    `\Gamma_1(N)` are computed instead.  Then the condition is that
    `\chi(-1)*\psi(-1) =(-1)^k`.

    If ``character`` is a list of integers, the parameters for `\Gamma_H(N)` are
    computed, where `H` is the subgroup of `(\ZZ/N\ZZ)^\times` generated by the
    integers in the given list.

    EXAMPLES::

        sage: sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(30)(1), 3)
        []

        sage: pars =  sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(30)(1), 4)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (1, 1), 1),
        ((1, 1), (1, 1), 2),
        ((1, 1), (1, 1), 3),
        ((1, 1), (1, 1), 5),
        ((1, 1), (1, 1), 6),
        ((1, 1), (1, 1), 10),
        ((1, 1), (1, 1), 15),
        ((1, 1), (1, 1), 30)]

        sage: pars = sage.modular.modform.eis_series.compute_eisenstein_params(15, 1)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, 1), 1),
        ((1, 1), (-1, 1), 5),
        ((1, 1), (1, zeta4), 1),
        ((1, 1), (1, zeta4), 3),
        ((1, 1), (-1, -1), 1),
        ((1, 1), (1, -zeta4), 1),
        ((1, 1), (1, -zeta4), 3),
        ((-1, 1), (1, -1), 1)]

        sage: sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(15).0, 1)
        [(Dirichlet character modulo 15 of conductor 1 mapping 11 |--> 1, 7 |--> 1, Dirichlet character modulo 15 of conductor 3 mapping 11 |--> -1, 7 |--> 1, 1),
        (Dirichlet character modulo 15 of conductor 1 mapping 11 |--> 1, 7 |--> 1, Dirichlet character modulo 15 of conductor 3 mapping 11 |--> -1, 7 |--> 1, 5)]

        sage: len(sage.modular.modform.eis_series.compute_eisenstein_params(GammaH(15, [4]), 3))
        8
    """
    if isinstance(character, (int, Integer)):
        return __find_eisen_chars_gamma1(character, k)
    elif isinstance(character, GammaH_class):
        return __find_eisen_chars_gammaH(character.level(), character._generators_for_H(), k)
    else:
        return __find_eisen_chars(character, k)
