from sage.misc.cachefunc import cached_method
from sage.schemes.hyperelliptic_curves_smooth_model import (
    hyperelliptic_generic,
    invariants,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_finite_field import (
    HyperellipticCurveSmoothModel_finite_field,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_padic_field import (
    HyperellipticCurveSmoothModel_padic_field,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_rational_field import (
    HyperellipticCurveSmoothModel_rational_field,
)

"""
TODO List

### Invariants

There seems to be massive redundancy in having these methods and the ones imported into
invariants. I think we should fix this by putting the methods themseleves into this class.
"""


class HyperellipticCurveSmoothModel_g2(
    hyperelliptic_generic.HyperellipticCurveSmoothModel_generic
):
    def is_odd_degree(self):
        """
        Return ``True`` if the curve is an odd degree model.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^5 - x^4 + 3
            sage: HyperellipticCurveSmoothModel(f).is_odd_degree()
            True
        """
        f, h = self.hyperelliptic_polynomials()
        df = f.degree()
        if h.degree() < 3:
            return df % 2 == 1
        elif df < 6:
            return False
        else:
            a0 = f.leading_coefficient()
            c0 = h.leading_coefficient()
            return (c0**2 + 4 * a0) == 0

    @cached_method
    def jacobian(self):
        from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_g2_generic import (
            HyperellipticJacobian_g2_generic,
        )

        return HyperellipticJacobian_g2_generic(self)

    # -----------------------------------
    # Genus Two invariant computations
    #
    # TODO: should we move logic from functions in `invariants()`
    # into this file rather than import them
    # -----------------------------------

    def clebsch_invariants(self):
        r"""
        Return the Clebsch invariants `(A, B, C, D)` of Mestre, p 317, [Mes1991]_.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.invariants`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^5 - x^4 + 3
            sage: HyperellipticCurveSmoothModel(f).clebsch_invariants()
            (0, -2048/375, -4096/25, -4881645568/84375)
            sage: HyperellipticCurveSmoothModel(f(2*x)).clebsch_invariants()
            (0, -8388608/375, -1073741824/25, -5241627016305836032/84375)

            sage: HyperellipticCurveSmoothModel(f, x).clebsch_invariants()
            (-8/15, 17504/5625, -23162896/140625, -420832861216768/7119140625)
            sage: HyperellipticCurveSmoothModel(f(2*x), 2*x).clebsch_invariants()
            (-512/15, 71696384/5625, -6072014209024/140625, -451865844002031331704832/7119140625)

        TESTS::

            sage: # optional - magma
            sage: magma(HyperellipticCurveSmoothModel(f)).ClebschInvariants()
            [ 0, -2048/375, -4096/25, -4881645568/84375 ]
            sage: magma(HyperellipticCurveSmoothModel(f(2*x))).ClebschInvariants()
            [ 0, -8388608/375, -1073741824/25, -5241627016305836032/84375 ]
            sage: magma(HyperellipticCurveSmoothModel(f, x)).ClebschInvariants()
            [ -8/15, 17504/5625, -23162896/140625, -420832861216768/7119140625 ]
            sage: magma(HyperellipticCurveSmoothModel(f(2*x), 2*x)).ClebschInvariants()
            [ -512/15, 71696384/5625, -6072014209024/140625, -451865844002031331704832/7119140625 ]
        """
        f, h = self.hyperelliptic_polynomials()
        return invariants.clebsch_invariants(4 * f + h**2)

    def igusa_clebsch_invariants(self):
        r"""
        Return the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.invariants`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^5 - x + 2
            sage: HyperellipticCurveSmoothModel(f).igusa_clebsch_invariants()
            (-640, -20480, 1310720, 52160364544)
            sage: HyperellipticCurveSmoothModel(f(2*x)).igusa_clebsch_invariants()
            (-40960, -83886080, 343597383680, 56006764965979488256)

            sage: HyperellipticCurveSmoothModel(f, x).igusa_clebsch_invariants()
            (-640, 17920, -1966656, 52409511936)
            sage: HyperellipticCurveSmoothModel(f(2*x), 2*x).igusa_clebsch_invariants()
            (-40960, 73400320, -515547070464, 56274284941110411264)

        TESTS::

            sage: # optional - magma
            sage: magma(HyperellipticCurveSmoothModel(f)).IgusaClebschInvariants()
            [ -640, -20480, 1310720, 52160364544 ]
            sage: magma(HyperellipticCurveSmoothModel(f(2*x))).IgusaClebschInvariants()
            [ -40960, -83886080, 343597383680, 56006764965979488256 ]
            sage: magma(HyperellipticCurveSmoothModel(f, x)).IgusaClebschInvariants()
            [ -640, 17920, -1966656, 52409511936 ]
            sage: magma(HyperellipticCurveSmoothModel(f(2*x), 2*x)).IgusaClebschInvariants()
            [ -40960, 73400320, -515547070464, 56274284941110411264 ]
        """
        f, h = self.hyperelliptic_polynomials()
        return invariants.igusa_clebsch_invariants(4 * f + h**2)

    def absolute_igusa_invariants_wamelen(self):
        r"""
        Return the three absolute Igusa invariants used by van Wamelen [Wam1999]_.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: HyperellipticCurveSmoothModel(x^5 - 1).absolute_igusa_invariants_wamelen()
            (0, 0, 0)
            sage: HyperellipticCurveSmoothModel((x^5 - 1)(x - 2), (x^2)(x - 2)).absolute_igusa_invariants_wamelen()
            (0, 0, 0)
        """
        f, h = self.hyperelliptic_polynomials()
        return invariants.absolute_igusa_invariants_wamelen(4 * f + h**2)

    def absolute_igusa_invariants_kohel(self):
        r"""
        Return the three absolute Igusa invariants used by Kohel [KohECHIDNA]_.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.invariants`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: HyperellipticCurveSmoothModel(x^5 - 1).absolute_igusa_invariants_kohel()
            (0, 0, 0)
            sage: HyperellipticCurveSmoothModel(x^5 - x + 1, x^2).absolute_igusa_invariants_kohel()
            (-1030567/178769, 259686400/178769, 20806400/178769)
            sage: HyperellipticCurveSmoothModel((x^5 - x + 1)(3*x + 1), (x^2)(3*x + 1)).absolute_igusa_invariants_kohel()
            (-1030567/178769, 259686400/178769, 20806400/178769)
        """
        f, h = self.hyperelliptic_polynomials()
        return invariants.absolute_igusa_invariants_kohel(4 * f + h**2)


class HyperellipticCurveSmoothModel_g2_padic_field(
    HyperellipticCurveSmoothModel_g2, HyperellipticCurveSmoothModel_padic_field
):
    pass


class HyperellipticCurveSmoothModel_g2_finite_field(
    HyperellipticCurveSmoothModel_g2, HyperellipticCurveSmoothModel_finite_field
):
    pass


class HyperellipticCurveSmoothModel_g2_rational_field(
    HyperellipticCurveSmoothModel_g2, HyperellipticCurveSmoothModel_rational_field
):
    pass
