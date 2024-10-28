# sage.doctest: needs sage.libs.flint
r"""
`q`-expansion of `j`-invariant
"""
from .eis_series import eisenstein_series_qexp
from .vm_basis import delta_qexp
from sage.rings.rational_field import QQ


def j_invariant_qexp(prec=10, K=QQ):
    r"""
    Return the `q`-expansion of the `j`-invariant to
    precision ``prec`` in the field `K`.

    .. SEEALSO::

        If you want to evaluate (numerically) the `j`-invariant at
        certain points, see the special function :func:`elliptic_j`.

    .. WARNING::

        Stupid algorithm -- we divide by Delta, which is slow.

    EXAMPLES::

        sage: j_invariant_qexp(4)
        q^-1 + 744 + 196884*q + 21493760*q^2 + 864299970*q^3 + O(q^4)
        sage: j_invariant_qexp(2)
        q^-1 + 744 + 196884*q + O(q^2)
        sage: j_invariant_qexp(100, GF(2))
        q^-1 + q^7 + q^15 + q^31 + q^47 + q^55 + q^71 + q^87 + O(q^100)
    """
    if prec <= -1:
        raise ValueError("the prec must be nonnegative.")
    prec += 2
    g6 = -504*eisenstein_series_qexp(6, prec, K=QQ)
    Delta = delta_qexp(prec).change_ring(QQ)
    j = (g6*g6) * (~Delta) + 1728
    if K != QQ:
        return j.change_ring(K)
    else:
        return j


# NOTE: this needs to be sped up.  The pari code src/basemath/trans3.c is
# faster.
