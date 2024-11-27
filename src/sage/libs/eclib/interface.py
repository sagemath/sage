r"""
Sage interface to Cremona's ``eclib`` library (also known as ``mwrank``)

This is the Sage interface to John Cremona's ``eclib`` C++ library for
arithmetic on elliptic curves.  The classes defined in this module
give Sage interpreter-level access to some of the functionality of
``eclib``.  For most purposes, it is not necessary to directly use these
classes. Instead, one can create an
:class:`EllipticCurve <sage.schemes.elliptic_curves.constructor.EllipticCurve>`
and call methods that are implemented using this module.

.. NOTE::

   This interface is a direct library-level interface to ``eclib``,
   including the 2-descent program ``mwrank``.

TESTS:

Check that ``eclib`` is imported as needed::

    sage: [k for k in sys.modules if k.startswith("sage.libs.eclib")]
    []
    sage: EllipticCurve('11a1').mwrank_curve()
    y^2 + y = x^3 - x^2 - 10 x - 20
    sage: [k for k in sys.modules if k.startswith("sage.libs.eclib")]
    ['...']
"""
import sys
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import IntegerRing

from .mwrank import _Curvedata, _two_descent, _mw, parse_point_list

class mwrank_EllipticCurve(SageObject):
    r"""
    The :class:`mwrank_EllipticCurve` class represents an elliptic
    curve using the ``Curvedata`` class from ``eclib``, called here an 'mwrank
    elliptic curve'.

    Create the mwrank elliptic curve with invariants
    ``ainvs``, which is a list of 5 or less *integers* `a_1`,
    `a_2`, `a_3`, `a_4`, and `a_5`.

    If strictly less than 5 invariants are given, then the *first*
    ones are set to 0, so, e.g., ``[3,4]`` means `a_1=a_2=a_3=0` and
    `a_4=3`, `a_5=4`.

    INPUT:

    - ``ainvs`` -- list or tuple list of 5 or less integers, the
      coefficients of a nonsingular Weierstrass equation

    - ``verbose`` -- boolean (default: ``False``); verbosity flag.  If ``True``,
      then all Selmer group computations will be verbose.

    EXAMPLES:

    We create the elliptic curve `y^2 + y = x^3 + x^2 - 2x`::

        sage: e = mwrank_EllipticCurve([0, 1, 1, -2, 0])
        sage: e.ainvs()
        [0, 1, 1, -2, 0]

    This example illustrates that omitted `a`-invariants default to `0`::

        sage: e = mwrank_EllipticCurve([3, -4])
        sage: e
        y^2 = x^3 + 3 x - 4
        sage: e.ainvs()
        [0, 0, 0, 3, -4]

    The entries of the input list are coerced to :class:`int`.
    If this is impossible, then an error is raised::

        sage: e = mwrank_EllipticCurve([3, -4.8]); e
        Traceback (most recent call last):
        ...
        TypeError: ainvs must be a list or tuple of integers.

    When you enter a singular model you get an exception::

        sage: e = mwrank_EllipticCurve([0, 0])
        Traceback (most recent call last):
        ...
        ArithmeticError: Invariants (= 0,0,0,0,0) do not describe an elliptic curve.
    """

    def __init__(self, ainvs, verbose=False):
        r"""
        Create the mwrank elliptic curve with invariants
        ``ainvs``, which is a list of 5 or less *integers* `a_1`,
        `a_2`, `a_3`, `a_4`, and `a_5`.

        See the docstring of this class for full documentation.

        EXAMPLES:

        We create the elliptic curve `y^2 + y = x^3 + x^2 - 2x`::

            sage: e = mwrank_EllipticCurve([0, 1, 1, -2, 0])
            sage: e.ainvs()
            [0, 1, 1, -2, 0]
        """

        ainvs = list(ainvs)
        if len(ainvs) > 5:
            raise TypeError("ainvs must have length at most 5")

        # Pad ainvs on the beginning by 0's, so e.g. [a4, a6] works
        ainvs = [0] * (5 - len(ainvs)) + ainvs

        # Convert each entry to an int
        try:
            a_int = [IntegerRing()(x) for x in ainvs]
        except (TypeError, ValueError):
            raise TypeError("ainvs must be a list or tuple of integers.")
        self.__ainvs = a_int
        self.__curve = _Curvedata(a_int[0], a_int[1], a_int[2],
                                  a_int[3], a_int[4])

        if verbose:
            self.__verbose = True
        else:
            self.__verbose = False

        # place holders
        self.__saturate = -2  # not yet saturated
        self.__descent = None

    def __reduce__(self):
        r"""
        Standard Python function used in pickling.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.__reduce__()
            (<class 'sage.libs.eclib.interface.mwrank_EllipticCurve'>, ([0, 0, 1, -7, 6], False))
        """
        return mwrank_EllipticCurve, (self.__ainvs, self.__verbose)

    def set_verbose(self, verbose):
        """
        Set the verbosity of printing of output by the :meth:`two_descent()` and
        other functions.

        INPUT:

        - ``verbose`` -- integer; if positive, print lots of output when
          doing 2-descent

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.saturate() # no output
            sage: E.gens()
            [[0, -1, 1]]

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.set_verbose(1)
            sage: E.saturate() # tol 1e-10
            Basic pair: I=48, J=-432
            disc=255744
            2-adic index bound = 2
            By Lemma 5.1(a), 2-adic index = 1
            2-adic index = 1
            One (I,J) pair
            Looking for quartics with I = 48, J = -432
            Looking for Type 2 quartics:
            Trying positive a from 1 up to 1 (square a first...)
            (1,0,-6,4,1)        --trivial
            Trying positive a from 1 up to 1 (...then non-square a)
            Finished looking for Type 2 quartics.
            Looking for Type 1 quartics:
            Trying positive a from 1 up to 2 (square a first...)
            (1,0,0,4,4) --nontrivial...(x:y:z) = (1 : 1 : 0)
            Point = [0:0:1]
                height = 0.0511114082399688402358
            Rank of B=im(eps) increases to 1 (The previous point is on the egg)
            Exiting search for Type 1 quartics after finding one which is globally soluble.
            Mordell rank contribution from B=im(eps) = 1
            Selmer  rank contribution from B=im(eps) = 1
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            Searching for points (bound = 8)...done:
              found points which generate a subgroup of rank 1
              and regulator 0.0511114082399688402358
            Processing points found during 2-descent...done:
              now regulator = 0.0511114082399688402358
            Saturating (with bound = -1)...done:
              points were already saturated.
        """
        self.__verbose = verbose

    def _curve_data(self):
        r"""
        Return the underlying :class:`_Curvedata` class for this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-1,0])
            sage: E._curve_data()
            [0,0,1,-1,0]
            b2 = 0       b4 = -2         b6 = 1  b8 = -1
            c4 = 48             c6 = -216
            disc = 37   (# real components = 2)
            #torsion not yet computed
        """
        return self.__curve

    def ainvs(self):
        r"""
        Return the `a`-invariants of this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-1,0])
            sage: E.ainvs()
            [0, 0, 1, -1, 0]
        """
        return self.__ainvs

    def isogeny_class(self, verbose=False):
        r"""
        Return the isogeny class of this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.isogeny_class()
            ([[0, -1, 1, 0, 0], [0, -1, 1, -10, -20], [0, -1, 1, -7820, -263580]], [[0, 5, 0], [5, 0, 5], [0, 5, 0]])
        """
        return self.__curve.isogeny_class(verbose)

    def __repr__(self):
        r"""
        Return the string representation of this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.__repr__()
            'y^2 + y = x^3 - x^2'
        """
        a1, a2, a3, a4, a6 = self.__ainvs
        # we do not assume a1, a2, a3 are reduced to {0,1}, {-1,0,1}, {0,1}

        def coeff(a):
            return ''.join([" +" if a > 0 else " -",
                            " " + str(abs(a)) if abs(a) > 1 else ""])

        return ''.join(['y^2',
                        ' '.join([coeff(a1), 'xy']) if a1 else '',
                        ' '.join([coeff(a3), 'y']) if a3 else '',
                        ' = x^3',
                        ' '.join([coeff(a2), 'x^2']) if a2 else '',
                        ' '.join([coeff(a4), 'x']) if a4 else '',
                        ' '.join([" +" if a6 > 0 else " -", str(abs(a6))]) if a6 else ''])

    def two_descent(self, verbose=True, selmer_only=False, first_limit=20,
                    second_limit=8, n_aux=-1, second_descent=True):
        r"""
        Compute 2-descent data for this curve.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``); print what mwrank is doing

        - ``selmer_only`` -- boolean (default: ``False``); ``selmer_only`` switch

        - ``first_limit`` -- integer (default: 20); bound on `|x|+|z|` in
          quartic point search

        - ``second_limit`` -- integer (default: 8); bound on
          `\log \max(|x|,|z|)`, i.e. logarithmic

        - ``n_aux`` -- integer (default: -1); (only relevant for general
          2-descent when 2-torsion trivial) number of primes used for
          quartic search.  ``n_aux=-1`` causes default (8) to be used.
          Increase for curves of higher rank.

        - ``second_descent`` -- boolean (default: ``True``); (only relevant
          for curves with 2-torsion, where mwrank uses descent via
          2-isogeny) flag determining whether or not to do second
          descent.  *Default strongly recommended.*

        OUTPUT: nothing

        TESTS:

        See :issue:`7992`::

            sage: EllipticCurve([0, prod(prime_range(10))]).mwrank_curve().two_descent()
            Basic pair: I=0, J=-5670
            disc=-32148900
            2-adic index bound = 2
            2-adic index = 2
            Two (I,J) pairs
            Looking for quartics with I = 0, J = -5670
            Looking for Type 3 quartics:
            Trying positive a from 1 up to 5 (square a first...)
            Trying positive a from 1 up to 5 (...then non-square a)
            (2,0,-12,19,-6)     --nontrivial...(x:y:z) = (2 : 4 : 1)
            Point = [-2488:-4997:512]
                height = 6.46767239...
            Rank of B=im(eps) increases to 1
            Trying negative a from -1 down to -3
            Finished looking for Type 3 quartics.
            Looking for quartics with I = 0, J = -362880
            Looking for Type 3 quartics:
            Trying positive a from 1 up to 20 (square a first...)
            Trying positive a from 1 up to 20 (...then non-square a)
            Trying negative a from -1 down to -13
            Finished looking for Type 3 quartics.
            Mordell rank contribution from B=im(eps) = 1
            Selmer  rank contribution from B=im(eps) = 1
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: EllipticCurve([0, prod(prime_range(100))]).mwrank_curve().two_descent()
            Traceback (most recent call last):
            ...
            RuntimeError: A 2-descent did not complete successfully.

        Calling this method twice does not cause a segmentation fault
        (see :issue:`10665`)::

            sage: E = EllipticCurve([1, 1, 0, 0, 528])
            sage: E.two_descent(verbose=False)
            True
            sage: E.two_descent(verbose=False)
            True
        """
        first_limit = int(first_limit)
        second_limit = int(second_limit)
        n_aux = int(n_aux)
        second_descent = int(second_descent)    # convert from bool to (int) 0 or 1
        self.__descent = _two_descent()
        self.__descent.do_descent(self.__curve,
                                  verbose,
                                  selmer_only,
                                  first_limit,
                                  second_limit,
                                  n_aux,
                                  second_descent)
        if not self.__descent.ok():
            raise RuntimeError("A 2-descent did not complete successfully.")
        self.__saturate = -2  # not yet saturated

    def __two_descent_data(self):
        r"""
        Return the 2-descent data for this elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E._mwrank_EllipticCurve__two_descent_data()
            <sage.libs.eclib.mwrank._two_descent object at ...>
        """
        if self.__descent is None:
            self.two_descent(self.__verbose)
        return self.__descent

    def conductor(self):
        """
        Return the conductor of this curve, computed using Cremona's
        implementation of Tate's algorithm.

        .. NOTE::

           This is independent of PARI's.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([1, 1, 0, -6958, -224588])
            sage: E.conductor()
            2310
        """
        return self.__curve.conductor()

    def rank(self):
        """
        Return the rank of this curve, computed using :meth:`two_descent()`.

        In general this may only be a lower bound for the rank; an
        upper bound may be obtained using the function :meth:`rank_bound()`.
        To test whether the value has been proved to be correct, use
        the method :meth:`certain()`.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.rank()
            0
            sage: E.certain()
            True

        ::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.rank()
            0
            sage: E.certain()
            False
        """
        return self.__two_descent_data().getrank()

    def rank_bound(self):
        """
        Return an upper bound for the rank of this curve, computed
        using :meth:`two_descent()`.

        If the curve has no 2-torsion, this is equal to the 2-Selmer
        rank.  If the curve has 2-torsion, the upper bound may be
        smaller than the bound obtained from the 2-Selmer rank minus
        the 2-rank of the torsion, since more information is gained
        from the 2-isogenous curve or curves.

        EXAMPLES:

        The following is the curve 960D1, which has rank 0,
        but Sha of order 4::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.rank_bound()
            0
            sage: E.rank()
            0

        In this case the rank was computed using a second descent,
        which is able to determine (by considering a 2-isogenous
        curve) that Sha is nontrivial.  If we deliberately stop the
        second descent, the rank bound is larger::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.two_descent(second_descent = False, verbose=False)
            sage: E.rank_bound()
            2

        In contrast, for the curve 571A, also with rank 0 and Sha
        of order 4, we only obtain an upper bound of 2::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.rank_bound()
            2

        In this case the value returned by :meth:`rank()` is only a
        lower bound in general (though this is correct)::

            sage: E.rank()
            0
            sage: E.certain()
            False
        """
        return self.__two_descent_data().getrankbound()

    def selmer_rank(self):
        r"""
        Return the rank of the 2-Selmer group of the curve.

        EXAMPLES:

        The following is the curve 960D1, which has rank 0, but Sha of
        order 4.  The 2-torsion has rank 2, and the Selmer rank is 3::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.selmer_rank()
            3

        Nevertheless, we can obtain a tight upper bound on the rank
        since a second descent is performed which establishes the
        2-rank of Sha::

            sage: E.rank_bound()
            0

        To show that this was resolved using a second descent, we do
        the computation again but turn off ``second_descent``::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.two_descent(second_descent = False, verbose=False)
            sage: E.rank_bound()
            2

        For the curve 571A, also with rank 0 and Sha of order 4,
        but with no 2-torsion, the Selmer rank is strictly greater
        than the rank::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.selmer_rank()
            2
            sage: E.rank_bound()
            2

        In cases like this with no 2-torsion, the rank upper bound is
        always equal to the 2-Selmer rank.  If we ask for the rank,
        all we get is a lower bound::

            sage: E.rank()
            0
            sage: E.certain()
            False
        """
        return self.__two_descent_data().getselmer()

    def regulator(self):
        r"""
        Return the regulator of the saturated Mordell-Weil group.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.regulator()
            0.05111140823996884
        """
        self.saturate()
        if not self.certain():
            raise RuntimeError("Unable to saturate Mordell-Weil group.")
        R = self.__two_descent_data().regulator()
        return float(R)

    def saturate(self, bound=-1, lower=2):
        """
        Compute the saturation of the Mordell-Weil group.

        INPUT:

        - ``bound`` -- integer (default: -1); if `-1`, saturate at *all*
          primes by computing a bound on the saturation index,
          otherwise saturate at all primes up to the minimum of
          ``bound`` and the saturation index bound

        - ``lower`` -- integer (default: 2); only saturate at primes not
          less than this

        EXAMPLES:

        Since the 2-descent automatically saturates at primes up to
        20, further saturation often has no effect::

            sage: E = mwrank_EllipticCurve([0, 0, 0, -1002231243161, 0])
            sage: E.gens()
            [[-1001107, -4004428, 1]]
            sage: E.saturate()
            sage: E.gens()
            [[-1001107, -4004428, 1]]

        Check that :issue:`18031` is fixed::

            sage: E = EllipticCurve([0,-1,1,-266,968])
            sage: Q1 = E([-1995,3674,125])
            sage: Q2 = E([157,1950,1])
            sage: E.saturation([Q1,Q2])
            ([(1 : -27 : 1), (157 : 1950 : 1)], 3, 0.801588644684981)
        """
        bound = int(bound)
        if self.__saturate < bound:
            self.__two_descent_data().saturate(bound, lower)
            self.__saturate = bound

    def gens(self) -> list:
        """
        Return a list of the generators for the Mordell-Weil group.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.gens()
            [[0, -1, 1]]
        """
        self.saturate()
        return parse_point_list(self.__two_descent_data().getbasis())

    def certain(self):
        r"""
        Return ``True`` if the last :meth:`two_descent()` call provably correctly
        computed the rank.  If :meth:`two_descent()` hasn't been
        called, then it is first called by :meth:`certain()`
        using the default parameters.

        The result is ``True`` if and only if the results of the methods
        :meth:`rank()` and :meth:`rank_bound()` are equal.

        EXAMPLES:

        A 2-descent does not determine `E(\QQ)` with certainty
        for the curve `y^2 + y = x^3 - x^2 - 120x - 2183`::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -120, -2183])
            sage: E.two_descent(False)
            ...
            sage: E.certain()
            False
            sage: E.rank()
            0

        The previous value is only a lower bound; the upper bound is greater::

            sage: E.rank_bound()
            2

        In fact the rank of `E` is actually 0 (as one could see by
        computing the `L`-function), but Sha has order 4 and the
        2-torsion is trivial, so mwrank cannot conclusively
        determine the rank in this case.
        """
        return bool(self.__two_descent_data().getcertain())

    # def fullmw(self):
    #     return self.__two_descent_data().getfullmw()

    def CPS_height_bound(self):
        r"""
        Return the Cremona-Prickett-Siksek height bound.

        This is a
        floating point number `B` such that if `P` is a point on the
        curve, then the naive logarithmic height `h(P)` is less than
        `B+\hat{h}(P)`, where `\hat{h}(P)` is the canonical height of
        `P`.

        .. warning::

           We assume the model is minimal!

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 0, -1002231243161, 0])
            sage: E.CPS_height_bound()
            14.163198527061496
            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.CPS_height_bound()
            0.0
        """
        return self.__curve.cps_bound()

    def silverman_bound(self):
        r"""
        Return the Silverman height bound.  This is a floating point
        number `B` such that if `P` is a point on the curve, then the
        naive logarithmic height `h(P)` is less than `B+\hat{h}(P)`,
        where `\hat{h}(P)` is the canonical height of `P`.

        .. warning::

           We assume the model is minimal!

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 0, -1002231243161, 0])
            sage: E.silverman_bound()
            18.29545210468247
            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.silverman_bound()
            6.284833369972403
        """
        return self.__curve.silverman_bound()


class mwrank_MordellWeil(SageObject):
    r"""
    The :class:`mwrank_MordellWeil` class represents a subgroup of a
    Mordell-Weil group.  Use this class to saturate a specified list
    of points on an :class:`mwrank_EllipticCurve`, or to search for
    points up to some bound.

    INPUT:

    - ``curve`` -- :class:`mwrank_EllipticCurve`; the underlying
      elliptic curve

    - ``verbose`` -- boolean (default: ``False``); verbosity flag (controls
      amount of output produced in point searches)

    - ``pp`` -- integer (default: 1); process points flag (if nonzero,
      the points found are processed, so that at all times only a
      `\ZZ`-basis for the subgroup generated by the points found
      so far is stored. If zero, no processing is done and all
      points found are stored).

    - ``maxr`` -- integer (default: 999); maximum rank (quit point
      searching once the points found generate a subgroup of this
      rank. Useful if an upper bound for the rank is already
      known).

    EXAMPLES::

        sage: E = mwrank_EllipticCurve([1,0,1,4,-6])
        sage: EQ = mwrank_MordellWeil(E)
        sage: EQ
        Subgroup of Mordell-Weil group: []
        sage: EQ.search(2)
        P1 = [0:1:0]     is torsion point, order 1
        P1 = [1:-1:1]    is torsion point, order 2
        P1 = [2:2:1]     is torsion point, order 3
        P1 = [9:23:1]    is torsion point, order 6

        sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
        sage: EQ = mwrank_MordellWeil(E)
        sage: EQ.search(2)
        P1 = [0:1:0]     is torsion point, order 1
        P1 = [-3:0:1]     is generator number 1
        ...
        P4 = [-91:804:343]       = -2*P1 + 2*P2 + 1*P3 (mod torsion)
        sage: EQ
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

    Example to illustrate the verbose parameter::

        sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
        sage: EQ = mwrank_MordellWeil(E, verbose=False)
        sage: EQ.search(1)
        sage: EQ
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

        sage: EQ = mwrank_MordellWeil(E, verbose=True)
        sage: EQ.search(1)
        P1 = [0:1:0]     is torsion point, order 1
        P1 = [-3:0:1]     is generator number 1
        saturating up to 20...Saturation index bound (for points of good reduction)  = 3
        Reducing saturation bound from given value 20 to computed index bound 3
        Tamagawa index primes are [ 2 ]...
        Checking saturation at [ 2 3 ]
        Checking 2-saturation
        Points were proved 2-saturated (max q used = 7)
        Checking 3-saturation
        Points were proved 3-saturated (max q used = 7)
        done
        P2 = [-2:3:1]     is generator number 2
        saturating up to 20...Saturation index bound (for points of good reduction)  = 4
        Reducing saturation bound from given value 20 to computed index bound 4
        Tamagawa index primes are [ 2 ]...
        Checking saturation at [ 2 3 ]
        Checking 2-saturation
        possible kernel vector = [1,1]
        This point may be in 2E(Q): [14:-52:1]
        ...and it is!
        Replacing old generator #1 with new generator [1:-1:1]
        Reducing index bound from 4 to 2
        Points have successfully been 2-saturated (max q used = 7)
        Index gain = 2^1
        done, index = 2.
        Gained index 2, new generators = [ [1:-1:1] [-2:3:1] ]
        P3 = [-14:25:8]   is generator number 3
        saturating up to 20...Saturation index bound (for points of good reduction)  = 3
        Reducing saturation bound from given value 20 to computed index bound 3
        Tamagawa index primes are [ 2 ]...
        Checking saturation at [ 2 3 ]
        Checking 2-saturation
        Points were proved 2-saturated (max q used = 11)
        Checking 3-saturation
        Points were proved 3-saturated (max q used = 13)
        done, index = 1.
        P4 = [-1:3:1]    = -1*P1 + -1*P2 + -1*P3 (mod torsion)
        P4 = [0:2:1]     = 2*P1 + 0*P2 + 1*P3 (mod torsion)
        P4 = [2:13:8]    = -3*P1 + 1*P2 + -1*P3 (mod torsion)
        P4 = [1:0:1]     = -1*P1 + 0*P2 + 0*P3 (mod torsion)
        P4 = [2:0:1]     = -1*P1 + 1*P2 + 0*P3 (mod torsion)
        P4 = [18:7:8]    = -2*P1 + -1*P2 + -1*P3 (mod torsion)
        P4 = [3:3:1]     = 1*P1 + 0*P2 + 1*P3 (mod torsion)
        P4 = [4:6:1]     = 0*P1 + -1*P2 + -1*P3 (mod torsion)
        P4 = [36:69:64]  = 1*P1 + -2*P2 + 0*P3 (mod torsion)
        P4 = [68:-25:64]         = -2*P1 + -1*P2 + -2*P3 (mod torsion)
        P4 = [12:35:27]  = 1*P1 + -1*P2 + -1*P3 (mod torsion)
        sage: EQ
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

    Example to illustrate the process points (``pp``) parameter::

        sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
        sage: EQ = mwrank_MordellWeil(E, verbose=False, pp=1)
        sage: EQ.search(1); EQ # generators only
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]
        sage: EQ = mwrank_MordellWeil(E, verbose=False, pp=0)
        sage: EQ.search(1); EQ # all points found
        Subgroup of Mordell-Weil group: [[-3:0:1], [-2:3:1], [-14:25:8], [-1:3:1], [0:2:1], [2:13:8], [1:0:1], [2:0:1], [18:7:8], [3:3:1], [4:6:1], [36:69:64], [68:-25:64], [12:35:27]]
    """

    def __init__(self, curve, verbose=True, pp=1, maxr=999):
        r"""
        Constructor for the :class:`mwrank_MordellWeil` class.

        See the docstring of this class for full documentation.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([1,0,1,4,-6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ
            Subgroup of Mordell-Weil group: []
        """
        if not isinstance(curve, mwrank_EllipticCurve):
            raise TypeError("curve (=%s) must be an mwrank_EllipticCurve" % curve)
        self.__curve = curve
        self.__verbose = verbose
        self.__pp = pp
        self.__maxr = maxr
        if verbose:
            verb = 1
        else:
            verb = 0
        self.__mw = _mw(curve._curve_data(), verb, pp, maxr)

    def __reduce__(self):
        r"""
        Standard Python function used in pickling.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.__reduce__()
            (<class 'sage.libs.eclib.interface.mwrank_MordellWeil'>, (y^2 + y = x^3 - 7 x + 6, True, 1, 999))
        """
        return mwrank_MordellWeil, (self.__curve, self.__verbose, self.__pp, self.__maxr)

    def __repr__(self):
        r"""
        String representation of this Mordell-Weil subgroup.

        OUTPUT:

        (string) String representation of this Mordell-Weil subgroup.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E, verbose=False)
            sage: EQ.__repr__()
            'Subgroup of Mordell-Weil group: []'
            sage: EQ.search(1)
            sage: EQ.__repr__()
            'Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]'
        """
        return "Subgroup of Mordell-Weil group: %s" % self.__mw

    def process(self, v, saturation_bound=0):
        """Process points in the list ``v``.

        This function allows one to add points to a :class:`mwrank_MordellWeil` object.

        INPUT:

        - ``v`` -- list of 3-tuples or lists of ints or Integers; a
          list of triples of integers, which define points on the
          curve

        - ``saturation_bound`` -- integer (default: 0); saturate at primes up to
          ``saturation_bound``, or at *all* primes if ``saturation_bound`` is
          -1. When ``saturation_bound`` is 0 (the default), do no saturation.

        OUTPUT:

        None.  But note that if the ``verbose`` flag is set, then there
        will be some output as a side-effect.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.gens()
            [[1, -1, 1], [-2, 3, 1], [-14, 25, 8]]
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1, -1, 1], [-2, 3, 1], [-14, 25, 8]])
            P1 = [1:-1:1]         is generator number 1
            P2 = [-2:3:1]         is generator number 2
            P3 = [-14:25:8]       is generator number 3

        ::

            sage: EQ.points()
            [[1, -1, 1], [-2, 3, 1], [-14, 25, 8]]

        Example to illustrate the saturation parameter ``saturation_bound``::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], saturation_bound=20)
            P1 = [1547:-2967:343]         is generator number 1
            ...
            Gained index 5, new generators = [ [-2:3:1] [-14:25:8] [1:-1:1] ]

            sage: EQ.points()
            [[-2, 3, 1], [-14, 25, 8], [1, -1, 1]]

        Here the processing was followed by saturation at primes up to
        20.  Now we prevent this initial saturation::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], saturation_bound=0)
            P1 = [1547:-2967:343]         is generator number 1
            P2 = [2707496766203306:864581029138191:2969715140223272]      is generator number 2
            P3 = [-13422227300:-49322830557:12167000000]          is generator number 3
            sage: EQ.points()
            [[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]]
            sage: EQ.regulator()
            375.42920288254555
            sage: EQ.saturate(2)  # points were not 2-saturated
            saturating basis...Saturation index bound (for points of good reduction)  = 93
            Only p-saturating for p up to given value 2.
            The resulting points may not be p-saturated for p between this and the computed index bound 93
            Tamagawa index primes are [ 2 ]...
            Checking saturation at [ 2 ]
            Checking 2-saturation
            possible kernel vector = [1,0,0]
            This point may be in 2E(Q): [1547:-2967:343]
            ...and it is!
            Replacing old generator #1 with new generator [-2:3:1]
            Reducing index bound from 93 to 46
            Points have successfully been 2-saturated (max q used = 11)
            Index gain = 2^1
            done
            Gained index 2
            New regulator =  93.85730072
            (True, 2, '[ ]')
            sage: EQ.points()
            [[-2, 3, 1], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]]
            sage: EQ.regulator()
            93.85730072063639
            sage: EQ.saturate(3)  # points were not 3-saturated
            saturating basis...Saturation index bound (for points of good reduction)  = 46
            Only p-saturating for p up to given value 3.
            The resulting points may not be p-saturated for p between this and the computed index bound 46
            Tamagawa index primes are [ 2 ]...
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            possible kernel vector = [0,1,0]
            This point may be in 3E(Q): [2707496766203306:864581029138191:2969715140223272]
            ...and it is!
            Replacing old generator #2 with new generator [-14:25:8]
            Reducing index bound from 46 to 15
            Points have successfully been 3-saturated (max q used = 13)
            Index gain = 3^1
            done
            Gained index 3
            New regulator =  10.42858897
            (True, 3, '[ ]')
            sage: EQ.points()
            [[-2, 3, 1], [-14, 25, 8], [-13422227300, -49322830557, 12167000000]]
            sage: EQ.regulator()
            10.4285889689596
            sage: EQ.saturate(5)  # points were not 5-saturated
            saturating basis...Saturation index bound (for points of good reduction)  = 15
            Only p-saturating for p up to given value 5.
            The resulting points may not be p-saturated for p between this and the computed index bound 15
            Tamagawa index primes are [ 2 ]...
            Checking saturation at [ 2 3 5 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            Checking 5-saturation
            possible kernel vector = [0,0,1]
            This point may be in 5E(Q): [-13422227300:-49322830557:12167000000]
            ...and it is!
            Replacing old generator #3 with new generator [1:-1:1]
            Reducing index bound from 15 to 3
            Points have successfully been 5-saturated (max q used = 71)
            Index gain = 5^1
            done
            Gained index 5
            New regulator =  0.4171435588
            (True, 5, '[ ]')
            sage: EQ.points()
            [[-2, 3, 1], [-14, 25, 8], [1, -1, 1]]
            sage: EQ.regulator()
            0.417143558758384
            sage: EQ.saturate()   # points are now saturated
            saturating basis...Saturation index bound (for points of good reduction)  = 3
            Tamagawa index primes are [ 2 ]...
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            done
            (True, 1, '[ ]')
        """
        if not isinstance(v, list):
            raise TypeError("v (=%s) must be a list" % v)
        saturation_bound = int(saturation_bound)
        for P in v:
            if not isinstance(P, (list, tuple)) or len(P) != 3:
                raise TypeError("v (=%s) must be a list of 3-tuples (or 3-element lists) of ints" % v)
            self.__mw.process(P, saturation_bound)

    def regulator(self):
        """
        Return the regulator of the points in this subgroup of
        the Mordell-Weil group.

        .. NOTE::

           ``eclib`` can compute the regulator to arbitrary precision,
           but the interface currently returns the output as a ``float``.

        OUTPUT:

        (float) The regulator of the points in this subgroup.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.regulator()
            1.0

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.regulator()
            0.417143558758384
        """
        return self.__mw.regulator()

    def rank(self):
        """
        Return the rank of this subgroup of the Mordell-Weil group.

        OUTPUT: integer

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.rank()
            0

        A rank 3 example::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.rank()
            0
            sage: EQ.regulator()
            1.0

        The preceding output is correct, since we have not yet tried
        to find any points on the curve either by searching or
        2-descent::

            sage: EQ
            Subgroup of Mordell-Weil group: []

        Now we do a very small search::

            sage: EQ.search(1)
            P1 = [0:1:0]         is torsion point, order 1
            P1 = [-3:0:1]         is generator number 1
            saturating up to 20...Checking 2-saturation
            ...
            P4 = [12:35:27]      = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ
            Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]
            sage: EQ.rank()
            3
            sage: EQ.regulator()
            0.417143558758384

        We do in fact now have a full Mordell-Weil basis.
        """
        return self.__mw.rank()

    def saturate(self, max_prime=-1, min_prime=2):
        r"""Saturate this subgroup of the Mordell-Weil group.

        INPUT:

        - ``max_prime`` -- integer (default: -1); if `-1`, an
          upper bound is computed for the primes at which the subgroup
          may not be saturated, and saturation is performed for all
          primes up to this bound.  Otherwise, the bound used is the
          minimum of ``max_prime`` and the computed bound.

        - ``min_prime`` -- integer (default: 2); only do saturation at
          primes no less than this.  (For example, if the points have
          been found via :meth:`two_descent` they should already be
          2-saturated so a value of 3 is appropriate.)

        OUTPUT:

        (3-tuple) (``ok``, ``index``, ``unsatlist``) where:

        - ``ok`` -- boolean; ``True`` if and only if the saturation was
          provably successful at all primes attempted.  If the default
          was used for ``max_prime``, then ``True`` indicates that the
          subgroup is saturated at *all* primes.

        - ``index`` -- integer; the index of the group generated by the
          original points in their saturation

        - ``unsatlist`` -- list of ints list of primes at which
          saturation could not be proved or achieved

        .. NOTE::

          In versions up to v20190909, ``eclib`` used floating point
          methods based on elliptic logarithms to divide points, and
          did not compute the precision necessary, which could cause
          failures. Since v20210310, ``eclib`` uses exact method based
          on division polynomials, which should mean that such
          failures does not happen.

        .. NOTE::

           We emphasize that if this function returns ``True`` as the
           first return argument (``ok``), and if the default was used
           for the parameter ``max_prime``, then the points in the
           basis after calling this function are saturated at *all*
           primes, i.e., saturating at the primes up to ``max_prime``
           are sufficient to saturate at all primes.  Note that the
           function computes an upper bound for the index of
           saturation, and does no work for primes greater than this
           even if ``max_prime`` is larger.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)

        We initialise with three points which happen to be 2, 3 and 5
        times the generators of this rank 3 curve.  To prevent
        automatic saturation at this stage we set the parameter
        ``sat`` to 0 (which is in fact the default)::

            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], saturation_bound=0)
            P1 = [1547:-2967:343]         is generator number 1
            P2 = [2707496766203306:864581029138191:2969715140223272]      is generator number 2
            P3 = [-13422227300:-49322830557:12167000000]          is generator number 3
            sage: EQ
            Subgroup of Mordell-Weil group: [[1547:-2967:343], [2707496766203306:864581029138191:2969715140223272], [-13422227300:-49322830557:12167000000]]
            sage: EQ.regulator()
            375.42920288254555

        Now we saturate at `p=2`, and gain index 2::

            sage: EQ.saturate(2)  # points were not 2-saturated
            saturating basis...Saturation index bound (for points of good reduction) = 93
            Only p-saturating for p up to given value 2.
            ...
            Gained index 2
            New regulator =  93.857...
            (True, 2, '[ ]')
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [2707496766203306:864581029138191:2969715140223272], [-13422227300:-49322830557:12167000000]]
            sage: EQ.regulator()
            93.85730072063639

        Now we saturate at `p=3`, and gain index 3::

            sage: EQ.saturate(3)  # points were not 3-saturated
            saturating basis...Saturation index bound (for points of good reduction) = 46
            Only p-saturating for p up to given value 3.
            ...
            Gained index 3
            New regulator =  10.428...
            (True, 3, '[ ]')
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [-14:25:8], [-13422227300:-49322830557:12167000000]]
            sage: EQ.regulator()
            10.4285889689596

        Now we saturate at `p=5`, and gain index 5::

            sage: EQ.saturate(5)  # points were not 5-saturated
            saturating basis...Saturation index bound (for points of good reduction) = 15
            Only p-saturating for p up to given value 5.
            ...
            Gained index 5
            New regulator =  0.417...
            (True, 5, '[ ]')
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [-14:25:8], [1:-1:1]]
            sage: EQ.regulator()
            0.417143558758384

        Finally we finish the saturation.  The output here shows that
        the points are now provably saturated at all primes::

            sage: EQ.saturate()   # points are now saturated
            saturating basis...Saturation index bound (for points of good reduction) = 3
            Tamagawa index primes are [ 2 ]...
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            done
            (True, 1, '[ ]')

        Of course, the :meth:`process()` function would have done all this
        automatically for us::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], saturation_bound=5)
            P1 = [1547:-2967:343]         is generator number 1
            ...
            Gained index 5, new generators = [ [-2:3:1] [-14:25:8] [1:-1:1] ]
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [-14:25:8], [1:-1:1]]
            sage: EQ.regulator()
            0.417143558758384

        But we would still need to use the :meth:`saturate()` function to
        verify that full saturation has been done::

            sage: EQ.saturate()
            saturating basis...Saturation index bound (for points of good reduction) = 3
            Tamagawa index primes are [ 2 ]...
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            done
            (True, 1, '[ ]')

        Note the output of the preceding command: it proves that the
        index of the points in their saturation is at most 3, then
        proves saturation at 2 and at 3, by reducing the points modulo
        all primes of good reduction up to 11, respectively 13.
        """
        ok, index, unsat = self.__mw.saturate(int(max_prime), int(min_prime))
        return bool(ok), int(str(index)), unsat

    def search(self, height_limit=18, verbose=False):
        r"""
        Search for new points, and add them to this subgroup of the
        Mordell-Weil group.

        INPUT:

        - ``height_limit``-- float (default: 18); search up to this
          logarithmic height

        .. NOTE::

          On 32-bit machines, this *must* be < 21.48 (`31\log(2)`) else
          `\exp(h_{\text{lim}}) > 2^{31}` and overflows.  On 64-bit machines, it
          must be *at most* 43.668  (`63\log(2)`) .  However, this bound is a logarithmic
          bound and increasing it by just 1 increases the running time
          by (roughly) `\exp(1.5)=4.5`, so searching up to even 20
          takes a very long time.

        .. NOTE::

           The search is carried out with a quadratic sieve, using
           code adapted from a version of Michael Stoll's
           ``ratpoints`` program.  It would be preferable to use a
           newer version of ``ratpoints``.

        - ``verbose`` -- boolean (default: ``False``)turn verbose operation on
          or off

        EXAMPLES:

        A rank 3 example, where a very small search is sufficient to
        find a Mordell-Weil basis::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.search(1)
            P1 = [0:1:0]         is torsion point, order 1
            P1 = [-3:0:1]         is generator number 1
            ...
            P4 = [12:35:27]      = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ
            Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

        In the next example, a search bound of 12 is needed to find a
        non-torsion point::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -18392, -1186248]) #1056g4
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.search(11); EQ
            P1 = [0:1:0]         is torsion point, order 1
            P1 = [161:0:1]       is torsion point, order 2
            Subgroup of Mordell-Weil group: []
            sage: EQ.search(12); EQ
            P1 = [0:1:0]         is torsion point, order 1
            P1 = [161:0:1]       is torsion point, order 2
            P1 = [4413270:10381877:27000]         is generator number 1
            ...
            Subgroup of Mordell-Weil group: [[4413270:10381877:27000]]
        """
        height_limit = float(height_limit)
        int_bits = sys.maxsize.bit_length()
        max_height_limit = int_bits * 0.693147  # log(2.0) = 0.693147 approx
        if height_limit >= max_height_limit:
            raise ValueError("The height limit must be < {} = {}log(2) on a {}-bit machine.".format(max_height_limit, int_bits, int_bits + 1))

        moduli_option = 0  # Use Stoll's sieving program... see strategies in ratpoints-1.4.c

        verbose = bool(verbose)
        self.__mw.search(height_limit, moduli_option, verbose)

    def points(self):
        """
        Return a list of the generating points in this Mordell-Weil
        group.

        OUTPUT: list of lists of length 3, each holding the
        primitive integer coordinates `[x,y,z]` of a generating
        point

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.search(1)
            P1 = [0:1:0]         is torsion point, order 1
            P1 = [-3:0:1]         is generator number 1
            ...
            P4 = [12:35:27]      = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ.points()
            [[1, -1, 1], [-2, 3, 1], [-14, 25, 8]]
        """
        return self.__mw.getbasis()
