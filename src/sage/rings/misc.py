"Miscellaneous utilities"

###########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.structure.sequence import Sequence
from sage.categories.fields import Fields
_Fields = Fields()

def composite_field(K, L):
    """
    Return a canonical field that contains both $K$ and $L$, if possible.
    Otherwise, raise a ValueError.

    INPUT:
        K -- field
        L -- field
    OUTPUT:
        field

    EXAMPLES:
        sage: composite_field(QQ,QQbar)
        Algebraic Field
        sage: composite_field(QQ,QQ[sqrt(2)])
        Number Field in sqrt2 with defining polynomial x^2 - 2
        sage: composite_field(QQ,QQ)
        Rational Field
        sage: composite_field(QQ,GF(7))
        Traceback (most recent call last):
        ...
        ValueError: unable to find a common field
    """
    C = Sequence([K(0), L(0)]).universe()
    if C not in _Fields:
        raise ValueError("unable to find a common field")
    return C
