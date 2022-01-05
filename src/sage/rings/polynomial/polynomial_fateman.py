"Polynomial multiplication by Kronecker substitution"

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
################################################################################

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ

# Faster than Sage's
from math import log as pylog
from math import ceil as pyceil


def _mul_fateman_to_int2(f_list, g_list):
    """
    Convert a polynomial to an integer by evaluating it
    INPUT: p, a list of integers
    OUTPUT: padding
    """
    max_coeff_f = max([abs(i) for i in f_list])
    max_coeff_g = max([abs(i) for i in g_list])
    b = (1+min(len(f_list),len(g_list)))*max_coeff_f*max_coeff_g
    return int(pyceil(pylog(b, 2)))


def _mul_fateman_to_poly(number,padding):
    """
    Converts a number to a polynomial, according
    to a padding
    OUTPUT: a list containing the coefficient of
    a polynomial of degree len(list)

    """
    coeffs = []
    flag=0
    append = coeffs.append
    if number < 0:
        number = -number
        flag=1

    while number > 0:
        r =  number%(1<<padding)
        number = (number-r) >> padding
        if r > (1<<(padding-1)):
            r -= 1<<padding
            number+=1
        append(r)

    if flag==1:
        return [-c for c in coeffs]
    return coeffs

def _mul_fateman_mul(f,g):
    """
    Multiply 2 polynomials
    """

    f=f.change_ring(QQ)
    g=g.change_ring(QQ)

    f_list = f.list()
    g_list = g.list()

    # If these polynomials have real
    # coefficients, convert them to
    # rational coefficients.
    # Note: no precision is lost in this
    # direction

    fgcd = f_list[0].content(f_list)
    ggcd = g_list[0].content(g_list)

    # Need to change ring to ZZ
    z_poly_f=(f*fgcd.denominator()).change_ring(ZZ)
    z_poly_g=(g*ggcd.denominator()).change_ring(ZZ)

    div = 1/(fgcd.denominator()*ggcd.denominator())

    z_poly_f_list = z_poly_f.coefficients(sparse=False)
    z_poly_g_list = z_poly_g.coefficients(sparse=False)
    padding = _mul_fateman_to_int2(z_poly_f_list,z_poly_g_list)

    n_f = z_poly_f(1<<padding)
    n_g = z_poly_g(1<<padding)

    if div == 1:
        return _mul_fateman_to_poly(n_f*n_g,padding)
    #return to_poly(n_f*n_g,padding)
    else:
        l=_mul_fateman_to_poly(n_f*n_g,padding)
        return [QQ(i*div) for i in l]

