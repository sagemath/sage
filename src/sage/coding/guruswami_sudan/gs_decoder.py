# sage.doctest: needs sage.modules sage.rings.finite_rings sage.symbolic
r"""
Guruswami-Sudan decoder for (Generalized) Reed-Solomon codes

REFERENCES:

- [GS1999]_

- [Nie2013]_

AUTHORS:

- Johan S. R. Nielsen, original implementation (see [Nie]_ for details)
- David Lucas, ported the original implementation in Sage
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2015 Johan S. R. Nielsen <jsrn@jsrn.dk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.arith.misc import integer_floor as floor
from sage.coding.grs_code import GeneralizedReedSolomonCode
from sage.rings.integer_ring import ZZ
from sage.coding.decoder import Decoder
from sage.coding.guruswami_sudan.interpolation import gs_interpolation_linalg, gs_interpolation_lee_osullivan
from sage.coding.guruswami_sudan.utils import (johnson_radius,
                                               gilt,
                                               solve_degree2_to_integer_range)
from sage.functions.other import floor
from sage.misc.functional import sqrt


def n_k_params(C, n_k):
    r"""
    Internal helper function for the :class:`GRSGuruswamiSudanDecoder` class for
    allowing to specify either a GRS code `C` or the length and dimensions `n,
    k` directly, in all the static functions.

    If neither `C` or `n,k` were specified to those functions, an appropriate
    error should be raised. Otherwise, `n, k` of the code or the supplied tuple
    directly is returned.

    INPUT:

    - ``C`` -- a GRS code or ``None``

    - ``n_k`` -- tuple `(n,k)` being length and dimension of a GRS code, or
      ``None``

    OUTPUT:

    - ``n_k`` -- tuple `(n,k)` being length and dimension of a GRS code

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.gs_decoder import n_k_params
        sage: n_k_params(None, (10, 5))
        (10, 5)
        sage: C = codes.GeneralizedReedSolomonCode(GF(11).list()[:10], 5)
        sage: n_k_params(C, None)
        (10, 5)
        sage: n_k_params(None,None)
        Traceback (most recent call last):
        ...
        ValueError: Please provide either the code or its length and dimension
        sage: n_k_params(C, (12, 2))
        Traceback (most recent call last):
        ...
        ValueError: Please provide only the code or its length and dimension
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        return C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        return n_k


def roth_ruckenstein_root_finder(p, maxd=None, precision=None):
    """
    Wrapper for Roth-Ruckenstein algorithm to compute the roots of a polynomial
    with coefficients in `F[x]`.

    TESTS::

        sage: from sage.coding.guruswami_sudan.gs_decoder import roth_ruckenstein_root_finder
        sage: R.<x> = GF(13)[]
        sage: S.<y> = R[]
        sage: p = (y - x^2 - x - 1) * (y + x + 1)
        sage: roth_ruckenstein_root_finder(p, maxd = 2)
        [12*x + 12, x^2 + x + 1]
    """
    gens = p.parent().gens()
    if len(gens) == 2:
        p = p.polynomial(gens[1])
    return p.roots(multiplicities=False, degree_bound=maxd, algorithm='Roth-Ruckenstein')


def alekhnovich_root_finder(p, maxd=None, precision=None):
    """
    Wrapper for Alekhnovich's algorithm to compute the roots of a polynomial
    with coefficients in `F[x]`.

    TESTS::

        sage: from sage.coding.guruswami_sudan.gs_decoder import alekhnovich_root_finder
        sage: R.<x> = GF(13)[]
        sage: S.<y> = R[]
        sage: p = (y - x^2 - x - 1) * (y + x + 1)
        sage: alekhnovich_root_finder(p, maxd = 2)
        [12*x + 12, x^2 + x + 1]
    """
    gens = p.parent().gens()
    if len(gens) == 2:
        p = p.polynomial(gens[1])
    return p.roots(multiplicities=False, degree_bound=maxd, algorithm='Alekhnovich')


class GRSGuruswamiSudanDecoder(Decoder):
    r"""
    The Guruswami-Sudan list-decoding algorithm for decoding Generalized
    Reed-Solomon codes.

    The Guruswami-Sudan algorithm is a polynomial time algorithm to decode
    beyond half the minimum distance of the code. It can decode up to the
    Johnson radius which is `n - \sqrt(n(n-d))`, where `n, d` is the length,
    respectively minimum distance of the RS code. See [GS1999]_ for more details.
    It is a list-decoder meaning that it returns a list of all closest codewords
    or their corresponding message polynomials. Note that the output of the
    ``decode_to_code`` and ``decode_to_message`` methods are therefore lists.

    The algorithm has two free parameters, the list size and the multiplicity,
    and these determine how many errors the method will correct: generally,
    higher decoding radius requires larger values of these parameters. To decode
    all the way to the Johnson radius, one generally needs values in the order
    of `O(n^2)`, while decoding just one error less requires just `O(n)`.

    This class has static methods for computing choices of parameters given the
    decoding radius or vice versa.

    The Guruswami-Sudan consists of two computationally intensive steps:
    Interpolation and Root finding, either of which can be completed in multiple
    ways. This implementation allows choosing the sub-algorithms among currently
    implemented possibilities, or supplying your own.

    INPUT:

    - ``code`` -- a code associated to this decoder

    - ``tau`` -- integer (default: ``None``); the number of errors one wants the
      Guruswami-Sudan algorithm to correct

    - ``parameters`` -- (default: ``None``) a pair of integers, where:

      - the first integer is the multiplicity parameter, and
      - the second integer is the list size parameter.

    - ``interpolation_alg`` -- (default: ``None``) the interpolation algorithm
      that will be used. The following possibilities are currently available:

      * ``'LinearAlgebra'`` -- uses a linear system solver.
      * ``'LeeOSullivan'`` -- uses Lee O'Sullivan method based on row reduction
        of a matrix
      * ``None`` -- one of the above will be chosen based on the size of the
        code and the parameters.

      You can also supply your own function to perform the interpolation. See
      NOTE section for details on the signature of this function.

    - ``root_finder`` -- (default: ``None``) the rootfinding algorithm that will
      be used. The following possibilities are currently available:

      * ``'Alekhnovich'`` -- uses Alekhnovich's algorithm.

      * ``'RothRuckenstein'`` -- uses Roth-Ruckenstein algorithm.

      * ``None`` -- one of the above will be chosen based on the size of the
        code and the parameters.

      You can also supply your own function to perform the interpolation. See
      NOTE section for details on the signature of this function.

    .. NOTE::

        One has to provide either ``tau`` or ``parameters``. If neither are given,
        an exception will be raised.

        If one provides a function as ``root_finder``, its signature has to be:
        ``my_rootfinder(Q, maxd=default_value, precision=default_value)``. `Q`
        will be given as an element of `F[x][y]`. The function must return the
        roots as a list of polynomials over a univariate polynomial ring. See
        :meth:`roth_ruckenstein_root_finder` for an example.

        If one provides a function as ``interpolation_alg``, its signature has
        to be: ``my_inter(interpolation_points, tau, s_and_l, wy)``. See
        :meth:`sage.coding.guruswami_sudan.interpolation.gs_interpolation_linalg`
        for an example.

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, tau=97); D
        Guruswami-Sudan decoder for [250, 70, 181] Reed-Solomon Code over GF(251)
         decoding 97 errors with parameters (1, 2)

    One can specify multiplicity and list size instead of ``tau``::

        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, parameters=(1,2)); D
        Guruswami-Sudan decoder for [250, 70, 181] Reed-Solomon Code over GF(251)
         decoding 97 errors with parameters (1, 2)

    One can pass a method as ``root_finder`` (works also for ``interpolation_alg``)::

        sage: from sage.coding.guruswami_sudan.gs_decoder import roth_ruckenstein_root_finder
        sage: rf = roth_ruckenstein_root_finder
        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, parameters=(1,2),
        ....:                                             root_finder=rf);  D
        Guruswami-Sudan decoder for [250, 70, 181] Reed-Solomon Code over GF(251)
         decoding 97 errors with parameters (1, 2)

    If one wants to use the native Sage algorithms for the root finding step,
    one can directly pass the string given in the ``Input`` block of this class.
    This works for ``interpolation_alg`` as well::


        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, parameters=(1,2),
        ....:                                             root_finder='RothRuckenstein'); D
        Guruswami-Sudan decoder for [250, 70, 181] Reed-Solomon Code over GF(251)
         decoding 97 errors with parameters (1, 2)

    Actually, we can construct the decoder from ``C`` directly::

        sage: D = C.decoder("GuruswamiSudan", tau=97); D
        Guruswami-Sudan decoder for [250, 70, 181] Reed-Solomon Code over GF(251)
         decoding 97 errors with parameters (1, 2)
    """

    ####################### static methods ###############################

    @staticmethod
    def parameters_given_tau(tau, C=None, n_k=None):
        r"""
        Return the smallest possible multiplicity and list size given the
        given parameters of the code and decoding radius.

        INPUT:

        - ``tau`` -- integer; number of errors one wants the Guruswami-Sudan
          algorithm to correct
        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a pair of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`

        OUTPUT:

        - ``(s, l)`` -- a pair of integers, where:

          - ``s`` is the multiplicity parameter, and
          - ``l`` is the list size parameter.

        .. NOTE::

            One should to provide either ``C`` or ``(n, k)``. If neither or both
            are given, an exception will be raised.

        EXAMPLES::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: tau, n, k = 97, 250, 70
            sage: GSD.parameters_given_tau(tau, n_k=(n, k))
            (1, 2)

        Another example with a bigger decoding radius::

            sage: tau, n, k = 118, 250, 70
            sage: GSD.parameters_given_tau(tau, n_k=(n, k))
            (47, 89)

        Choosing a decoding radius which is too large results in an errors::

            sage: tau = 200
            sage: GSD.parameters_given_tau(tau, n_k=(n, k))
            Traceback (most recent call last):
            ...
            ValueError: The decoding radius must be less than
            the Johnson radius (which is 118.66)
        """
        n,k = n_k_params(C, n_k)

        johnson = johnson_radius(n, n - k + 1)
        if tau >= johnson:
            raise ValueError("The decoding radius must be less than the Johnson radius (which is %.2f)"
                             % float(johnson))

        # We start with l=1 and check if a satisfiable s can be chosen. We keep
        # increasing l by 1 until this is the case. The governing equation is
        #   s*(s+1)/2 * n < (l+1)*s*(n-tau) - l*(l+1)/2*(k-1)
        # See [GS1999]_
        def try_l(l):
            (mins,maxs) = solve_degree2_to_integer_range(n, n-2*(l+1)*(n-tau), (k-1)*l*(l+1))
            if maxs > 0 and maxs >= mins:
                return max(1, mins)
            else:
                return None
        s, l = None, 0
        while s is None:
            l += 1
            s = try_l(l)

        return (s, l)

    @staticmethod
    def guruswami_sudan_decoding_radius(C=None, n_k=None, l=None, s=None):
        r"""
        Return the maximal decoding radius of the Guruswami-Sudan decoder and
        the parameter choices needed for this.

        If ``s`` is set but ``l`` is not it will return the best decoding radius using this ``s``
        alongside with the required ``l``. Vice versa for ``l``. If both are
        set, it returns the decoding radius given this parameter choice.

        INPUT:

        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a pair of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`
        - ``s`` -- integer (default: ``None``); the multiplicity parameter of
          Guruswami-Sudan algorithm
        - ``l`` -- integer (default: ``None``); the list size parameter

        .. NOTE::

            One has to provide either ``C`` or ``n_k``. If none or both are
            given, an exception will be raised.

        OUTPUT:

        - ``(tau, (s, l))`` -- where

          - ``tau`` is the obtained decoding radius, and

          - ``(s, l)`` are the multiplicity parameter and the list size
            parameter giving the radius

        EXAMPLES::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: n, k = 250, 70
            sage: GSD.guruswami_sudan_decoding_radius(n_k=(n, k))
            (118, (47, 89))

        One parameter can be restricted at a time::

            sage: n, k = 250, 70
            sage: GSD.guruswami_sudan_decoding_radius(n_k=(n, k), s=3)
            (109, (3, 5))
            sage: GSD.guruswami_sudan_decoding_radius(n_k=(n, k), l=7)
            (111, (4, 7))

        The function can also just compute the decoding radius given the parameters::

            sage: GSD.guruswami_sudan_decoding_radius(n_k=(n, k), s=2, l=6)
            (92, (2, 6))
        """
        n,k = n_k_params(C, n_k)

        def get_tau(s, l):
            "Return the decoding radius given this s and l"
            if s <= 0 or l <= 0:
                return -1
            return gilt(n - n/2*(s+1)/(l+1) - (k-1)/2*l/s)
        if l is None and s is None:
            tau = gilt(johnson_radius(n, n - k + 1))
            return (tau, GRSGuruswamiSudanDecoder.parameters_given_tau(tau, n_k=(n, k)))
        if l is not None and s is not None:
            return (get_tau(s,l), (s,l))

        # Either s or l is set, but not both. First a shared local function
        def find_integral_max(real_max, f):
            """Given a real (local) maximum of a function `f`, return that of
            the integers around ``real_max`` which gives the (local) integral
            maximum, and the value of at that point."""
            if real_max in ZZ:
                int_max = ZZ(real_max)
                return (int_max, f(int_max))
            else:
                x_f = floor(real_max)
                x_c = x_f + 1
                f_f, f_c = f(x_f), f(x_c)
                return (x_f, f_f) if f_f >= f_c else (x_c, f_c)

        if s is not None:
            # maximising tau under condition
            # n*(s+1 choose 2) < (ell+1)*s*(n-tau) - (ell+1 choose 2)*(k-1)
            # knowing n and s, we can just minimise
            # ( n*(s+1 choose 2) + (ell+1 choose 2)*(k-1) )/(ell+1)
            # Differentiating and setting to zero yields ell best choice:
            lmax = sqrt(n*s*(s+1.)/(k-1.)) - 1.
            #the best integral value will be
            (l,tau) = find_integral_max(lmax, lambda l: get_tau(s,l))
            #Note that we have not proven that this ell is minimal in integral
            #sense! It just seems that this most often happens
            return (tau,(s,l))
        if l is not None:
            # Acquired similarly to when restricting s
            smax = sqrt((k-1.)/n*l*(l+1.))
            (s,tau) = find_integral_max(smax, lambda s: get_tau(s,l))
            return (tau, (s,l))

    @staticmethod
    def _suitable_parameters_given_tau(tau, C=None, n_k=None):
        r"""
        Return quite good multiplicity and list size parameters for the code
        parameters and the decoding radius.

        These parameters are not guaranteed to be the best ones possible
        for the provided ``tau``, but arise from easily-evaluated closed
        expressions and are very good approximations of the best ones.

        See [Nie2013]_ pages 53-54, proposition 3.11 for details.

        INPUT:

        - ``tau`` -- integer; number of errors one wants the Guruswami-Sudan
          algorithm to correct
        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a pair of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`

        OUTPUT:

        - ``(s, l)`` -- a pair of integers, where:

          - ``s`` is the multiplicity parameter, and
          - ``l`` is the list size parameter.

        .. NOTE::

            One has to provide either ``C`` or ``(n, k)``. If neither or both
            are given, an exception will be raised.

        EXAMPLES:


        The following is an example where the parameters are optimal::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: tau = 98
            sage: n, k = 250, 70
            sage: GSD._suitable_parameters_given_tau(tau, n_k=(n, k))
            (2, 3)
            sage: GSD.parameters_given_tau(tau, n_k=(n, k))
            (2, 3)

        This is an example where they are not::

            sage: tau = 97
            sage: n, k = 250, 70
            sage: GSD._suitable_parameters_given_tau(tau, n_k=(n, k))
            (2, 3)
            sage: GSD.parameters_given_tau(tau, n_k=(n, k))
            (1, 2)

        We can provide a GRS code instead of `n` and `k` directly::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: GSD._suitable_parameters_given_tau(tau, C=C)
            (2, 3)

        Another one with a bigger ``tau``::

            sage: GSD._suitable_parameters_given_tau(118, C=C)
            (47, 89)
        """
        n,k = n_k_params(C, n_k)
        w = k - 1
        atau = n - tau
        smin = tau * w / (atau ** 2 - n * w)
        s = floor(1 + smin)
        D = (s - smin) * (atau ** 2 - n * w) * s + (w**2) / 4
        l = floor(atau / w * s + 0.5 - sqrt(D)/w)
        return (s, l)

    @staticmethod
    def gs_satisfactory(tau, s, l, C=None, n_k=None):
        r"""
        Return whether input parameters satisfy the governing equation of
        Guruswami-Sudan.

        See [Nie2013]_ page 49, definition 3.3 and proposition 3.4 for details.

        INPUT:

        - ``tau`` -- integer; number of errors one expects Guruswami-Sudan algorithm
          to correct
        - ``s`` -- integer; multiplicity parameter of Guruswami-Sudan algorithm
        - ``l`` -- integer; list size parameter
        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`

        .. NOTE::

            One has to provide either ``C`` or ``(n, k)``. If none or both are
            given, an exception will be raised.

        EXAMPLES::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: tau, s, l = 97, 1, 2
            sage: n, k = 250, 70
            sage: GSD.gs_satisfactory(tau, s, l, n_k=(n, k))
            True

        One can also pass a GRS code::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: GSD.gs_satisfactory(tau, s, l, C=C)
            True

        Another example where ``s`` and ``l`` does not satisfy the equation::

            sage: tau, s, l = 118, 47, 80
            sage: GSD.gs_satisfactory(tau, s, l, n_k=(n, k))
            False

        If one provides both ``C`` and ``n_k`` an exception is returned::

            sage: tau, s, l = 97, 1, 2
            sage: n, k = 250, 70
            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: GSD.gs_satisfactory(tau, s, l, C=C, n_k=(n, k))
            Traceback (most recent call last):
            ...
            ValueError: Please provide only the code or its length and dimension

        Same if one provides none of these::

            sage: GSD.gs_satisfactory(tau, s, l)
            Traceback (most recent call last):
            ...
            ValueError: Please provide either the code or its length and dimension
        """
        n,k = n_k_params(C, n_k)
        return l > 0 and s > 0 and n * s * (s+1) < (l+1) * (2*s*(n-tau) - (k-1) * l)

    ####################### decoder itself ###############################
    def __init__(self, code, tau=None, parameters=None, interpolation_alg=None, root_finder=None):
        r"""
        TESTS:

        If neither ``tau`` nor ``parameters`` is given, an exception is returned::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = GSD(C)
            Traceback (most recent call last):
            ...
            ValueError: Specify either tau or parameters

        If one provides something else than one of the allowed strings or a method as ``interpolation_alg``,
        an exception is returned::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = GSD(C, tau=97, interpolation_alg=42)
            Traceback (most recent call last):
            ...
            ValueError: Please provide a method or one of the allowed strings for interpolation_alg

        Same thing for ``root_finder``::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = GSD(C, tau=97, root_finder='FortyTwo')
            Traceback (most recent call last):
            ...
            ValueError: Please provide a method or one of the allowed strings for root_finder

        If one provides a full set of parameters (tau, s and l) which are not satisfactory, an
        error message is returned::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = GSD(C, tau=142, parameters=(1, 2))
            Traceback (most recent call last):
            ...
            ValueError: Impossible parameters for the Guruswami-Sudan algorithm

        If ``code`` is not a GRS code, an error is raised::

            sage: C  = codes.random_linear_code(GF(11), 10, 4)
            sage: GSD(C, tau=2)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a generalized Reed-Solomon code
        """
        if not isinstance(code, GeneralizedReedSolomonCode):
            raise ValueError("code has to be a generalized Reed-Solomon code")
        n, k = code.length(), code.dimension()
        if tau and parameters:
            if not GRSGuruswamiSudanDecoder.gs_satisfactory(tau, parameters[0], parameters[1], C=code):
                raise ValueError("Impossible parameters for the Guruswami-Sudan algorithm")
            self._tau, self._s, self._ell = tau, parameters[0], parameters[1]
        elif tau:
            self._tau = tau
            self._s, self._ell = GRSGuruswamiSudanDecoder.parameters_given_tau(tau, n_k=(n, k))
        elif parameters:
            self._s = parameters[0]
            self._ell = parameters[1]
            (self._tau,_) = GRSGuruswamiSudanDecoder.guruswami_sudan_decoding_radius(C=code, s=self._s, l=self._ell)
        else:
            raise ValueError("Specify either tau or parameters")
        if callable(interpolation_alg):
            self._interpolation_alg = interpolation_alg
        elif interpolation_alg is None or interpolation_alg == "LeeOSullivan":
            self._interpolation_alg = gs_interpolation_lee_osullivan
        elif interpolation_alg == "LinearAlgebra":
            self._interpolation_alg = gs_interpolation_linalg
        else:
            raise ValueError("Please provide a method or one of the allowed strings for interpolation_alg")
        if callable(root_finder):
            self._root_finder = root_finder
        elif root_finder == "RothRuckenstein":
            self._root_finder = roth_ruckenstein_root_finder
        elif root_finder is None or root_finder == "Alekhnovich":
            self._root_finder = alekhnovich_root_finder
        else:
            raise ValueError("Please provide a method or one of the allowed strings for root_finder")
        super().__init__(code, code.ambient_space(), "EvaluationPolynomial")

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D
            Guruswami-Sudan decoder for [250, 70, 181] Reed-Solomon Code over GF(251) decoding 97 errors with parameters (1, 2)
        """
        return "Guruswami-Sudan decoder for %s decoding %s errors with parameters %s" % (self.code(), self.decoding_radius(), (self.multiplicity(), self.list_size()))

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: latex(D)
            \textnormal{Guruswami-Sudan decoder for } [250, 70, 181] \textnormal{ Reed-Solomon Code over } \Bold{F}_{251}\textnormal{ decoding }97\textnormal{ errors with parameters }(1, 2)
        """
        return "\\textnormal{Guruswami-Sudan decoder for } %s\\textnormal{ decoding }%s\\textnormal{ errors with parameters }%s" % (self.code()._latex_(), self.decoding_radius(), (self.multiplicity(), self.list_size()))

    def __eq__(self, other):
        r"""
        Test equality between GRSGuruswamiSudanDecoder objects.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D1 = C.decoder("GuruswamiSudan", tau=97)
            sage: D2 = C.decoder("GuruswamiSudan", tau=97)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSGuruswamiSudanDecoder)\
                and self.code() == other.code()\
                and self.decoding_radius() == other.decoding_radius()\
                and self.multiplicity() == other.multiplicity()\
                and self.list_size() == other.list_size()\
                and self.interpolation_algorithm() == other.interpolation_algorithm()\
                and self.rootfinding_algorithm() == other.rootfinding_algorithm()

    def interpolation_algorithm(self):
        r"""
        Return the interpolation algorithm that will be used.

        Remember that its signature has to be:
        ``my_inter(interpolation_points, tau, s_and_l, wy)``.
        See :meth:`sage.coding.guruswami_sudan.interpolation.gs_interpolation_linalg`
        for an example.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D.interpolation_algorithm()
            <function gs_interpolation_lee_osullivan at 0x...>
        """
        return self._interpolation_alg

    def rootfinding_algorithm(self):
        r"""
        Return the rootfinding algorithm that will be used.

        Remember that its signature has to be:
        ``my_rootfinder(Q, maxd=default_value, precision=default_value)``.
        See :meth:`roth_ruckenstein_root_finder`
        for an example.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D.rootfinding_algorithm()
            <function alekhnovich_root_finder at 0x...>
        """
        return self._root_finder

    def parameters(self):
        r"""
        Return the multiplicity and list size parameters of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D.parameters()
            (1, 2)
        """
        return (self._s, self._ell)

    def multiplicity(self):
        r"""
        Return the multiplicity parameter of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D.multiplicity()
            1
        """
        return self._s

    def list_size(self):
        r"""
        Return the list size parameter of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D.list_size()
            2
        """
        return self._ell

    def decode_to_message(self, r):
        r"""
        Decode ``r`` to the list of polynomials whose encoding by
        :meth:`self.code()` is within Hamming distance
        :meth:`self.decoding_radius` of ``r``.

        INPUT:

        - ``r`` -- a received word, i.e. a vector in `F^n` where `F` and `n` are
          the base field respectively length of :meth:`self.code`

        EXAMPLES::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: C = codes.GeneralizedReedSolomonCode(GF(17).list()[:15], 6)
            sage: D = GSD(C, tau=5)
            sage: F.<x> = GF(17)[]
            sage: m = 13*x^4 + 7*x^3 + 10*x^2 + 14*x + 3
            sage: c = D.connected_encoder().encode(m)
            sage: r = vector(GF(17), [3,13,12,0,0,7,5,1,8,11,15,12,14,7,10])
            sage: (c-r).hamming_weight()
            5
            sage: messages = D.decode_to_message(r)
            sage: len(messages)
            2
            sage: m in messages
            True

        TESTS:

        If one has provided a method as a ``root_finder`` or a ``interpolation_alg`` which
        does not fit the allowed signature, an exception will be raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(17).list()[:15], 6)
            sage: D = GSD(C, tau=5, root_finder=next_prime)
            sage: F.<x> = GF(17)[]
            sage: m = 9*x^5 + 10*x^4 + 9*x^3 + 7*x^2 + 15*x + 2
            sage: c = D.connected_encoder().encode(m)
            sage: r = vector(GF(17), [3,1,4,2,14,1,0,4,13,12,1,16,1,13,15])
            sage: m in D.decode_to_message(r)
            Traceback (most recent call last):
            ...
            ValueError: The provided root-finding algorithm has a wrong signature.
            See the documentation of `...rootfinding_algorithm()` for details
        """
        return [self.connected_encoder().unencode(c) for c in self.decode_to_code(r)]

    def decode_to_code(self, r):
        r"""
        Return the list of all codeword within radius :meth:`self.decoding_radius` of the received word `r`.

        INPUT:

        - ``r`` -- a received word, i.e. a vector in `F^n` where `F` and `n` are
          the base field respectively length of :meth:`self.code`

        EXAMPLES::

            sage: GSD = codes.decoders.GRSGuruswamiSudanDecoder
            sage: C = codes.GeneralizedReedSolomonCode(GF(17).list()[:15], 6)
            sage: D = GSD(C, tau=5)
            sage: c = vector(GF(17), [3,13,12,0,0,7,5,1,8,11,1,9,4,12,14])
            sage: c in C
            True
            sage: r = vector(GF(17), [3,13,12,0,0,7,5,1,8,11,15,12,14,7,10])
            sage: r in C
            False
            sage: codewords = D.decode_to_code(r)
            sage: len(codewords)
            2
            sage: c in codewords
            True

        TESTS:

        Check that :issue:`21347` is fixed::

            sage: C = codes.GeneralizedReedSolomonCode(GF(13).list()[:10], 3)
            sage: D = GSD(C, tau=4)
            sage: c = vector(GF(13), [6, 8, 2, 1, 5, 1, 2, 8, 6, 9])
            sage: e = vector(GF(13), [1, 0, 0, 1, 1, 0, 0, 1, 0, 1])
            sage: D.decode_to_code(c+e)
            []
        """
        C = self.code()
        k = C.dimension()
        alphas = C.evaluation_points()
        colmults = C.column_multipliers()
        s = self.multiplicity()
        l = self.list_size()
        tau = self.decoding_radius()
        # SETUP INTERPOLATION PROBLEM
        wy = k-1
        points = [(alphas[i], r[i]/colmults[i]) for i in range(len(alphas))]
        # SOLVE INTERPOLATION
        try:
            Q = self.interpolation_algorithm()(points, tau, (s,l), wy)
        except TypeError:
            raise ValueError("The provided interpolation algorithm has a wrong signature. See the documentation of `codes.decoders.GRSGuruswamiSudanDecoder.interpolation_algorithm()` for details")
        # EXAMINE THE FACTORS AND CONVERT TO CODEWORDS
        try:
            polynomials = self.rootfinding_algorithm()(Q, maxd=wy)
        except TypeError:
            raise ValueError("The provided root-finding algorithm has a wrong signature. See the documentation of `codes.decoders.GRSGuruswamiSudanDecoder.rootfinding_algorithm()` for details")
        if not polynomials:
            return []

        E = self.connected_encoder()
        codewords = [ E.encode(f) for f in polynomials]
        # Root-finding might find spurious roots. Return only the ones which give nearby codewords
        return [ c for c in codewords if (r - c).hamming_weight() <= tau ]

    def decoding_radius(self):
        r"""
        Return the maximal number of errors that ``self`` is able to correct.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau=97)
            sage: D.decoding_radius()
            97

        An example where tau is not one of the inputs to the constructor::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", parameters=(2,4))
            sage: D.decoding_radius()
            105
        """
        return self._tau


####################### types ###############################

GeneralizedReedSolomonCode._registered_decoders["GuruswamiSudan"] = GRSGuruswamiSudanDecoder
GRSGuruswamiSudanDecoder._decoder_type = {"list-decoder", "always-succeed", "hard-decision"}
