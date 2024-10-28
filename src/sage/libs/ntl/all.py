r"""
Victor Shoup's NTL C++ Library

Sage provides an interface to Victor Shoup's C++ library NTL.
Features of this library include *incredibly fast* arithmetic with
polynomials and asymptotically fast factorization of polynomials.
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.libs.ntl.ntl_ZZ import (
    ntl_setSeed,
    ntl_ZZ as ZZ,
    randomBnd as ZZ_random,
    randomBits as ZZ_random_bits)

from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext as ZZ_pContext

from sage.libs.ntl.ntl_ZZ_p import (
    ntl_ZZ_p as ZZ_p,
    ntl_ZZ_p_random_element as ZZ_p_random)

from sage.libs.ntl.ntl_ZZX import (
    ntl_ZZX as ZZX,
    zero_ZZX, one_ZZX)

from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX as ZZ_pX

from sage.libs.ntl.ntl_ZZ_pEContext import ntl_ZZ_pEContext as ZZ_pEContext

from sage.libs.ntl.ntl_ZZ_pE import ntl_ZZ_pE as ZZ_pE

from sage.libs.ntl.ntl_ZZ_pEX import ntl_ZZ_pEX as ZZ_pEX

from sage.libs.ntl.ntl_lzz_pContext import ntl_zz_pContext as zz_pContext

from sage.libs.ntl.ntl_lzz_p import ntl_zz_p as zz_p

from sage.libs.ntl.ntl_lzz_pX import ntl_zz_pX as zz_pX

from sage.libs.ntl.ntl_mat_ZZ import ntl_mat_ZZ as mat_ZZ

from sage.libs.ntl.ntl_GF2 import ntl_GF2 as GF2

from sage.libs.ntl.ntl_GF2X import (
    ntl_GF2X as GF2X,
    GF2XHexOutput)

from sage.libs.ntl.ntl_GF2EContext import ntl_GF2EContext as GF2EContext

from sage.libs.ntl.ntl_GF2E import (
    ntl_GF2E as GF2E,
    ntl_GF2E_random as GF2E_random)

from sage.libs.ntl.ntl_GF2EX import ntl_GF2EX as GF2EX

from sage.libs.ntl.ntl_mat_GF2E import ntl_mat_GF2E as mat_GF2E

from sage.libs.ntl.ntl_mat_GF2 import ntl_mat_GF2 as mat_GF2
