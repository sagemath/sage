r"""
Dimensions of spaces of Siegel modular forms

Let `\Gamma_2` denote the integer symplectic group `\operatorname{Sp}(4, \ZZ)`.

This file is about the dimensions of the spaces of Siegel modular forms and
Siegel cusp forms for the group `\Gamma_2`
with respect to the representation `\det^k \operatorname{Sym}_j`
of `\operatorname{GL(2)}`.

They are denoted `A_{k,j}(\Gamma_2)` and `S_{k,j}(\Gamma_2)`.

Author:

Alex Ghitza

REFERENCES:

- [BFvdG2014] Bergstroem-Faber-van der Geer, *Siegel modular forms of
    degree three and the cohomology of local systems* Sel. Math. New Ser.
    (2014) 20:83-124.

- [Ibu1] Ibukiyama *Lifting conjectures from vector valued Siegel
    modular forms of degree two*

- [Ibu2] Ibukiyama *A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences*

- [IbWa] Ibukiyama and Wakatsuki

- [Peter] Petersen *Cohomology of
    local systems on the moduli of principally polarized abelian surfaces*

- [vdG] van der Geer, *Siegel modular forms*

- [Tsu1983] Tsushima, *An explicit dimension formula for the spaces
  of generalized automorphic forms with respect to Sp(2,Z)*
  Proc. Japan Acad. 59 Ser A (1983).

- [Take] Takemori *Structure theorems for vector valued Siegel
    modular forms of degree 2 and weight det^k otimes Sym(10)*

"""
from sage.modular.dims import dimension_cusp_forms, dimension_modular_forms
from sage.rings.integer_ring import ZZ
from sage.rings.lazy_series_ring import LazyPowerSeriesRing
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ


def dimension_cusp_forms_sp4z(k, J):
    r"""
    Return the dimension of the space of cusps forms `S_{k,j}(\Gamma_2)`.

    This uses Tsushima's formula from Theorem 4 of
    'An explicit dimension formula
    for the spaces of generalized automorphic forms with respect to Sp(2,Z)'
    Proc. Japan Acad. 59 Ser A (1983).

    Tsushima proves the correctness of his formula for (j = 0 and  k >= 4)
    or (j >= 1 and k >= 5), but Bergstroem-Faber-van der Geer prove that it
    holds for (j >= 0 and k >= 4), see page 97 of 'Siegel modular forms of
    degree three and the cohomology of local systems' Sel. Math. New Ser.
    (2014) 20:83-124.

    EXAMPLES::

        sage: dimension_cusp_forms_sp4z(5, 64)
        38
    """
    J = ZZ(J)
    if J < 0:
        raise ValueError("j cannot be negative")
    if J % 2:
        return ZZ.zero()
    j = J // 2

    k = ZZ(k)
    if k < 4:
        raise ValueError("not implemented for k < 4")

    res = (2 * j + 1) * (k - 2) * (2 * j + k - 1) * (2 * j + 2 * k - 3) / 2**7 / 3**3 / 5
    res -= (2 * j + 1) * (2 * j + 2 * k - 3) / 2**5 / 3**2
    res += (2 * j + 1) / 2**4 / 3
    res += (-1)**k * (7 * (k - 2) * (2 * j + k - 1) / 2**7 / 3**2 - (2 * j + 2 * k - 3) / 2**4 / 3 + ZZ(3) / 2**5)
    res += (-1)**j * (5 * (2 * j + 2 * k - 3) / 2**7 / 3 - ZZ.one() / 2**3)
    res += (-1)**(k + j) * (2 * j + 1) / 2**7

    i = CyclotomicField(4).gen()
    rho = CyclotomicField(3).gen()
    omega = CyclotomicField(5).gen()
    sigma = CyclotomicField(12).gen()

    res += (i**k * (i * (2 * j + k - 1) / 2**6 / 3 - i / 2**4)).trace()
    res += ((-1)**k * i**j / 2**5 * (i + 1)).trace()
    res += (i**k * (-1)**j * ((k - 2) / 2**6 / 3 - ZZ.one() / 2**4)).trace()
    res += ((-i)**k * i**j / 2**5 * (i + 1)).trace()
    res += ((-1)**k * rho**j / 3**3 * (rho + 1)).trace()
    res += (rho**k * rho**j / 2**2 / 3**4 * (2 * rho + 1) * (2 * j + 1)).trace()
    res += - (rho**k * (-rho)**j / 2**2 / 3**2 * (2 * rho + 1)).trace()
    res += ((-rho)**k * rho**j / 3**3).trace()
    res += (rho**j * ((1 - rho) / 2 / 3**4 * (2 * j + 2 * k - 3) - (1 - rho) / 2 / 3**2)).trace()
    res += (rho**k * ((rho + 2) / 2**3 / 3**4 * (2 * j + k - 1) - (5 * rho + 6) / 2**2 / 3**3)).trace()
    res += - ((-rho)**k * ((rho + 2) / 2**3 / 3**3 * (2 * j + k - 1) - (rho + 2) / 2**2 / 3**2)).trace()
    res += (rho**k * (rho**2)**j * ((1 - rho) / 2**3 / 3**4 * (k - 2) + (rho - 5) / 2**2 / 3**3)).trace()
    res += ((-rho)**k * (rho**2)**j * ((1 - rho) / 2**3 / 3**3 * (k - 2) - (1 - rho) / 2**2 / 3**2)).trace()
    res += (omega**k * (omega**4)**j / 5**2).trace()
    res += - (omega**k * (omega**3)**j / 5**2 * omega**2).trace()
    res += ((sigma**7)**k * (-1)**j / 2**3 / 3**2 * (sigma**2 + 1)).trace()
    res += - ((sigma**7)**k * (sigma**8)**j / 2**3 / 3**2 * (sigma + sigma**3)).trace()
    return ZZ(res)


def generating_series_cusp_forms_sp4z_wt3():
    """
    From page 45 of Petersen 'Cohomology of
    local systems on the moduli of principally polarized abelian surfaces'

    EXAMPLES::

        sage: generating_series_cusp_forms_sp4z_wt3()
        x^36 + x^42 + O(x^43)
    """
    x = LazyPowerSeriesRing(QQ, 'x').gen()
    denom = (1 - x**6) * (1 - x**8) * (1 - x**10) * (1 - x**12)
    return x**36 / denom


def generating_series_modular_forms_sp4z_wt4():
    """
    From page 8 of Ibukiyama 'Lifting conjectures from vector valued Siegel
    modular forms of degree two'

    EXAMPLES::

        sage: generating_series_modular_forms_sp4z_wt4()
        x^8 + x^12 + x^14 + O(x^15)
    """
    x = LazyPowerSeriesRing(QQ, 'x').gen()
    num = x**8 + x**12 - x**18 - x**20 - x**22 + x**28 + x**30 + x**32
    denom = (1 - x**6) * (1 - x**8) * (1 - x**10) * (1 - x**12)
    return num / denom


def generating_series_cusp_forms_sp4z_wt4():
    """
    From page 8 of Ibukiyama 'Lifting conjectures from vector valued Siegel
    modular forms of degree two'

    EXAMPLES::

        sage: generating_series_cusp_forms_sp4z_wt4()
        x^24 + x^28 + x^30 + O(x^31)
    """
    x = LazyPowerSeriesRing(QQ, 'x').gen()
    num = x**24 * (1 + x**4 + x**8 - x**10)
    denom = (1 - x**6) * (1 - x**8) * (1 - x**10) * (1 - x**12)
    return num / denom


def generating_series_cusp_forms_sp4z_wt5():
    r"""
    Formula for `S_{5,j}(\Gamma_2)`.

    From page 115 of Ibukiyama "A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences"

    EXAMPLES::

        sage: generating_series_cusp_forms_sp4z_wt5()
        x^18 + x^20 + 2*x^24 + O(x^25)
    """
    x = LazyPowerSeriesRing(QQ, 'x').gen()
    num = x**18 + x**20 + x**24
    denom = (1 - x**6) * (1 - x**8) * (1 - x**10) * (1 - x**12)
    return num / denom


def generating_series_cusp_forms_sp4z_wt7():
    r"""
    Formula for `S_{7,j}(\Gamma_2)`.

    From page 115 of Ibukiyama "A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences"

    EXAMPLES::

        sage: generating_series_cusp_forms_sp4z_wt7()
        x^12 + x^14 + x^16 + 2*x^18 + O(x^19)
    """
    x = LazyPowerSeriesRing(QQ, 'x').gen()
    num = x**12 + x**14 + x**16 + x**18 + x**20
    denom = (1 - x**6) * (1 - x**8) * (1 - x**10) * (1 - x**12)
    return num / denom


def generating_series_cusp_forms_sp4z_odd_wt(k, j):
    r"""
    Formula for `S_{2j+3,2k-6}(\Gamma_2)`.

    From page 115 of Ibukiyama "A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences"

    Note that there is a shift in the parameters j and k in Ibukiyama's formula
    """
    R = LazyPowerSeriesRing(QQ, ['x', 'y'])
    x, y = R.gens()
    denom = (1 - y**2) * (1 - y**3) * (1 - y**5) * (1 - y**6)
    denom *= (1 - x**3) * (1 - x**4) * (1 - x**5) * (1 - x**6)

    num = 1  # TODO, numerator in Table 1 of page 116
    # this needs work, the numerator seems to only be given for a particular
    # range in Ibukiyama
    return num / denom


def generating_series_modular_forms_sp4z_j10_even_wt():
    r"""
    Formula for `S_{k,10}(\Gamma_2)`.

    From Lemma 7.1 of Takemori 'Structure theorems for vector valued Siegel
    modular forms of degree 2 and weight det**k otimes Sym(10)'

    EXAMPLES::

        sage: generating_series_modular_forms_sp4z_j10_even_wt()
        y^6 + y^8 + 3*y^10 + 4*y^12 + O(y^13)
    """
    y = LazyPowerSeriesRing(QQ, 'y').gen()
    num = y**6 + y**8 + 2 * y**10 + 2 * y**12 + 3 * y**14 + 2 * y**16 + y**18 + y**20 - y**24 - y**26
    denom = (1 - y**4) * (1 - y**6) * (1 - y**10) * (1 - y**12)
    return num / denom


def generating_series_modular_forms_sp4z_j10_odd_wt():
    """
    From Lemma 7.1 of Takemori 'Structure theorems for vector valued Siegel
    modular forms of degree 2 and weight det**k otimes Sym(10)'

    EXAMPLES::

        sage: generating_series_modular_forms_sp4z_j10_odd_wt()
        y^9 + y^11 + 2*y^13 + 5*y^15 + O(y^16)
    """
    R = LazyPowerSeriesRing(QQ, 'y')
    y = R.gen()
    num = R([0] * 9 + [1, 0, 1, 0, 1, 0, 3, 0, 3, 0,
                       2, 0, 1, 0, 1, 0, 0, 0, -1, 0, -1])
    denom = (1 - y**4) * (1 - y**6) * (1 - y**10) * (1 - y**12)
    return num / denom


def generating_series_modular_forms_sp4z_j0_even_wt():
    """
    From Section 9 of van der Geer 'Siegel modular forms'

    EXAMPLES::

        sage: generating_series_modular_forms_sp4z_j0_even_wt()
        1 + y^4 + y^6 + O(y^7)
    """
    y = LazyPowerSeriesRing(QQ, 'y').gen()
    denom = (1 - y**4) * (1 - y**6) * (1 - y**10) * (1 - y**12)
    return 1 / denom


def generating_series_modular_forms_sp4z_j0_odd_wt():
    """
    From Section 9 of van der Geer 'Siegel modular forms'

    EXAMPLES::

        sage: generating_series_modular_forms_sp4z_j0_odd_wt()
        y^35 + y^39 + y^41 + O(y^42)
    """
    y = LazyPowerSeriesRing(QQ, 'y').gen()
    denom = (1 - y**4) * (1 - y**6) * (1 - y**10) * (1 - y**12)
    return y**35 / denom


def generating_function_cusp_forms_sp4z_k(k):
    """
    EXAMPLES::

        sage: generating_function_cusp_forms_sp4z_k(7)
        ?
    """
    R = LazyPowerSeriesRing(QQ, 'x')
    x = R.gen()
    denom = (1 - x**6) * (1 - x**8) * (1 - x**10) * (1 - x**12)
    return R(lambda j: dimension_cusp_forms_sp4z(k, j)) * denom


def generating_function_numerator_cusp_forms_sp4z_j(j, prec=None):

    def generating_function_denominator_cusp_forms_sp4z_j():
        x = PolynomialRing(ZZ).gen()
        return (1 - x**4) * (1 - x**6) * (1 - x**10) * (1 - x**12)

    denom = generating_function_denominator_cusp_forms_sp4z_j(j)
    if prec is None:
        prec = denom.degree()
    lst = [dimension_cusp_forms_sp4z(k, j) for k in range(prec)]
    R = PowerSeriesRing(QQ)
    res = (R(lst) * denom).add_bigoh(prec)
    return res.polynomial().change_ring(ZZ)


def dimension_modular_forms_sp4z(k, j):
    r"""
    Return the dimension of `A_{k,j}(\Gamma_2)`.

    EXAMPLES::

        sage: dimension_modular_forms_sp4z(5, 24)
        2
    """
    j = ZZ(j)
    if j < 0:
        raise ValueError("j cannot be negative")

    k = ZZ(k)
    if k < 0 or (k == 0 and j > 0):
        return ZZ.zero()

    if k % 2 or k == 2:
        return dimension_cusp_forms_sp4z(k, j)

    # if we are here then k is even and k >= 4
    res = dimension_cusp_forms_sp4z(k, j) + dimension_cusp_forms(1, k + j)
    if j == 0:
        # there's one Siegel Eisenstein form
        res += 1
    return res


def dimension_cusp_forms_Gamma_e(k, j):
    r"""
    What is `\Gamma_e` ?

    From Theorem 6.2 in Ibukiyama and Wakatsuki

    EXAMPLES::

        sage: dimension_cusp_forms_Gamma_e(5, 26)
        9
    """
    def mods(ts, ell, m):
        return ZZ(ts[(m % ell)])

    j2 = j // 2

    H1e = (j + 1) * (k - 2) * (j + k - 1) * (j + 2 * k - 3) / (2**6 * 3**3 * 5)
    H1u = - (j + 1) * (j + 2 * k - 3) / 2**6 / 3**2 + (j + 1) / 2**4 / 3
    H1 = H1e + H1u

    H2e = (-1)**k * (j + k - 1) * (k - 2) / 2**6 / 3**2
    H2qu = - (-1)**k * (j + 2 * k - 3) / 2**4 / 3 + 3 * ZZ(-1)**k / 2**6
    H2 = H2e + H2qu

    H3e = 0
    H3qu = - mods([(-1)**(j2), -1, (-1)**(j2 + 1), 1], 4, k) / 2**3 + mods([1, (-1)**(j2), -1, (-1)**(j2 + 1)], 4, k) / 2**4
    H3 = H3e + H3qu

    H4e = (mods([j + k - 1, -(j + k - 1), 0], 3, k) + mods([k - 2, 0, -(k - 2)], 3, j + k)) / 2**2 / 3**3
    H4qu = - (mods([1, -1, 0], 3, k) + mods([1, 0, -1], 3, j + k)) / 2**2 / 3**2 - (mods([0, -1, -1], 3, k) + mods([1, 1, 0], 3, j + k)) / 3**2
    H4 = H4e + H4qu

    H5e = (mods([-(j + k - 1), -(j + k - 1), 0, j + k - 1, j + k - 1, 0], 6, k) + mods([k - 2, 0, -(k - 2), -(k - 2), 0, k - 2], 6, j + k)) / 2**2 / 3**2
    H5qu = - (mods([-1, -1, 0, 1, 1, 0], 6, k) + mods([1, 0, -1, -1, 0, 1], 6, j + k)) / 2**2 / 3
    H5 = H5e + H5qu

    H6e = (-1)**j2 * (j + 2 * k - 3) / 2**6 + (-1)**(j2 + k) * (j + 1) / 2**6
    H6qu = - ZZ(-1)**j2 / 2**3
    H6 = H6e + H6qu

    H7e = (j + 2 * k - 3) / 3**3 * mods([1, -1, 0], 3, j) + (j + 1) * mods([0, 1, -1], 3, j + 2 * k) / 2 / 3**3
    H7qu = - mods([1, -1, 0], 3, j) / 2 / 3
    H7 = H7e + H7qu

    if j % 6 == 0:
        h9 = mods([1, 0, 0, -1, 0, 0], 6, k)
    elif j % 6 == 2:
        h9 = mods([-1, 1, 0, 1, -1, 0], 6, k)
    elif j % 6 == 4:
        h9 = mods([0, -1, 0, 0, 1, 0], 6, k)
    H9 = h9 / 2 / 3**2

    if j % 10 == 0:
        h10 = mods([1, 0, 0, -1, 0], 5, k)
    elif j % 10 == 2:
        h10 = mods([-1, 1, 0, 0, 0], 5, k)
    elif j % 10 == 4:
        h10 = 0
    elif j % 10 == 6:
        h10 = mods([0, 0, 0, 1, -1], 5, k)
    elif j % 10 == 8:
        h10 = mods([0, -1, 0, 0, 1], 5, k)
    H10 = h10 * 2 / 5

    if j % 8 == 0:
        h11 = mods([1, 1, -1, -1], 4, k)
    elif j % 8 == 2:
        h11 = mods([-1, 1, 1, -1], 4, k)
    elif j % 8 == 4:
        h11 = mods([-1, -1, 1, 1], 4, k)
    elif j % 8 == 6:
        h11 = mods([1, -1, -1, 1], 4, k)
    H11 = h11 / 2**3

    return H1 + H2 + H3 + H4 + H5 + H6 + H7 + H9 + H10 + H11


def dimension_cusp_forms_sp4z_sgn(k, j):
    """
    What is this ?
    """
    d1 = dimension_cusp_forms_Gamma_e(k, j)
    d2 = dimension_cusp_forms_sp4z(k, j)
    return d1 - d2


def dimension_V(k, j):
    """
    What is this ?
    """
    if not k % 2:
        d1 = dimension_cusp_forms(1, k + j / 2)
        d2 = sum(dimension_cusp_forms(1, k + j - 2 * a) * dimension_cusp_forms(1, k + 2 * a) for a in range(j / 2 + 1))
    else:
        d1 = - dimension_cusp_forms(1, k + j / 2)
        d2 = sum(dimension_cusp_forms(1, k - 1 + j - 2 * a) * dimension_cusp_forms(1, k + 1 + 2 * a) for a in range(j / 2))

    return (d1 + d2) // 2


def dimension_V2(k, j):
    """
    What is this ?
    """
    d1 = (-1)**k * dimension_cusp_forms(1, k + j / 2)
    d2 = sum(dimension_cusp_forms(1, k + j - a) * dimension_cusp_forms(1, k + a) for a in range(j + 1))
    return (d1 + d2) // 2


def dimension_W(k, j):
    """
    What is this ?
    """
    d1 = (-1)**(k + 1) * dimension_modular_forms(1, k + j / 2 - 6)
    d2 = sum(dimension_modular_forms(1, k + j - a - 6) * dimension_modular_forms(1, k + a - 6) for a in range(j + 1))
    return (d1 + d2) // 2
