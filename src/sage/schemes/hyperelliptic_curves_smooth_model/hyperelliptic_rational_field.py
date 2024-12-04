import sage.rings.abc
from sage.rings.padics.factory import Qp as pAdicField
from sage.schemes.hyperelliptic_curves_smooth_model import hyperelliptic_generic


class HyperellipticCurveSmoothModel_rational_field(
    hyperelliptic_generic.HyperellipticCurveSmoothModel_generic
):
    def __init__(self, projective_model, f, h, genus):
        super().__init__(projective_model, f, h, genus)

    def matrix_of_frobenius(self, p, prec=20):
        """
        Compute the matrix of Frobenius on Monsky-Washnitzer cohomology using
        the `p`-adic field with precision ``prec``.

        This function is essentially a wrapper function of
        :meth:`sage.schemes.hyperelliptic_curves.monsky_washnitzer.matrix_of_frobenius_hyperelliptic`.

        INPUT:

        - ``p`` (prime integer or pAdic ring / field ) -- if ``p`` is an integer,
          constructs a ``pAdicField`` with ``p`` to compute the matrix of
          Frobenius, otherwise uses the supplied pAdic ring or field.

        - ``prec`` (optional) -- if ``p`` is an prime integer, the `p`-adic
          precision of the coefficient ring constructed

        EXAMPLES::

            sage: K = pAdicField(5, prec=3)
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 2*x + 3)
            sage: H.matrix_of_frobenius(K)
            [            4*5 + O(5^3)       5 + 2*5^2 + O(5^3) 2 + 3*5 + 2*5^2 + O(5^3)     2 + 5 + 5^2 + O(5^3)]
            [      3*5 + 5^2 + O(5^3)             3*5 + O(5^3)             4*5 + O(5^3)         2 + 5^2 + O(5^3)]
            [    4*5 + 4*5^2 + O(5^3)     3*5 + 2*5^2 + O(5^3)       5 + 3*5^2 + O(5^3)     2*5 + 2*5^2 + O(5^3)]
            [            5^2 + O(5^3)       5 + 4*5^2 + O(5^3)     4*5 + 3*5^2 + O(5^3)             2*5 + O(5^3)]

        You can also pass directly a prime `p` with to construct a pAdic field with precision
        ``prec``::

            sage: H.matrix_of_frobenius(3, prec=2)
            [        O(3^2)     3 + O(3^2)         O(3^2)         O(3^2)]
            [    3 + O(3^2)         O(3^2)         O(3^2) 2 + 3 + O(3^2)]
            [  2*3 + O(3^2)         O(3^2)         O(3^2)    3^-1 + O(3)]
            [        O(3^2)         O(3^2)     3 + O(3^2)         O(3^2)]
        """
        import monsky_washnitzer  # TODO fix this

        if isinstance(p, (sage.rings.abc.pAdicField, sage.rings.abc.pAdicRing)):
            K = p
        else:
            K = pAdicField(p, prec)
        frob_p, _ = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(
            self.change_ring(K)
        )
        return frob_p

    def lseries(self, prec=53):
        """
        Return the L-series of this hyperelliptic curve of genus 2.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: C = HyperellipticCurveSmoothModel(x^2+x, x^3+x^2+1)
            sage: C.lseries()
            PARI L-function associated to Hyperelliptic Curve
            over Rational Field defined by y^2 + (x^3 + x^2 + 1)*y = x^2 + x
        """
        from sage.lfunctions.pari import LFunction

        L = LFunction(lfun_genus2(self), prec=prec)
        L.rename("PARI L-function associated to %s" % self)
        return L


#
#
# TODO: this is from `sage/src/sage/lfunctions/pari.py` but it uses the old method
# so I have copy-pasted it here for now with the new curve...
#
#


def lfun_genus2(C):
    """
    Return the L-function of a curve of genus 2.

    INPUT:

    - ``C`` -- hyperelliptic curve of genus 2

    Currently, the model needs to be minimal at 2.

    This uses :pari:`lfungenus2`.

    EXAMPLES::

        sage: from sage.lfunctions.pari import LFunction
        sage: x = polygen(QQ, 'x')
        sage: C = HyperellipticCurveSmoothModel(x^5 + x + 1)
        sage: L = LFunction(lfun_genus2(C))
        ...
        sage: L(3)
        0.965946926261520

        sage: C = HyperellipticCurveSmoothModel(x^2 + x, x^3 + x^2 + 1)
        sage: L = LFunction(lfun_genus2(C))
        sage: L(2)
        0.364286342944359

    TESTS::

        sage: x = polygen(QQ, 'x')
        sage: H = HyperellipticCurveSmoothModel(x^15 + x + 1)
        sage: L = LFunction(lfun_genus2(H))
        Traceback (most recent call last):
        ...
        ValueError: curve must be hyperelliptic of genus 2
    """
    from sage.libs.pari import pari
    from sage.schemes.hyperelliptic_curves_smooth_model import hyperelliptic_g2

    if not isinstance(C, hyperelliptic_g2.HyperellipticCurveSmoothModel_g2):
        raise ValueError("curve must be hyperelliptic of genus 2")
    P, Q = C.hyperelliptic_polynomials()
    return pari.lfungenus2(P) if not Q else pari.lfungenus2([P, Q])
