"""
`p`-adic Flat Polynomials
"""

# ****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense, Polynomial
from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
from sage.rings.infinity import infinity
from cypari2.gen import Gen as pari_gen
import sage.rings.padics.misc


class Polynomial_padic_flat(Polynomial_generic_dense, Polynomial_padic):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False, absprec=None):
        """
        TESTS:

        Check that :issue:`13620` has been fixed::

            sage: K = ZpFM(3)
            sage: R.<t> = K[]
            sage: R(R.zero())
            0
        """
        if x is None:
            Polynomial_generic_dense.__init__(self, parent, x, check, is_gen, construct)
            return
        R = parent.base_ring()
        if isinstance(x, sage.rings.fraction_field_element.FractionFieldElement):
            if x.denominator() != 1:
                raise TypeError("denominator must be 1")
            else:
                x = x.numerator()
        if isinstance(x, Polynomial):
            if x.base_ring() is R:
                x = list(x.list())
            else:
                x = [R(a) for a in x.list()]
        elif isinstance(x, list):
            if check:
                x = [R(a) for a in x]
        elif isinstance(x, dict):
            if check:
                m = infinity
                zero = R(0)
                n = max(x) if x else 0
                v = [zero] * (n + 1)
                for i, z in x.items():
                    v[i] = R(z)
                    m = min(m, v[i].precision_absolute())
                x = v
            else:
                m = sage.rings.padics.misc.min(a.precision_absolute()
                                               for a in x.values())
            if absprec is not None:
                m = min(m, absprec)
            Polynomial_generic_dense.__init__(self, parent, x, absprec=m)
            return
        elif isinstance(x, pari_gen):
            x = [R(w) for w in x.list()]
        else:
            x = [R(x)]
        if absprec is None:
            absprec = infinity
        m = min([a.precision_absolute() for a in x] + [absprec])
        Polynomial_generic_dense.__init__(self, parent, x, absprec=m)
