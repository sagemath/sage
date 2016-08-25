r"""
Generalized Reed-Solomon code

Given `n` different evaluation points `\alpha_1, \dots, \alpha_n` from some
finite field `F`, and `n` column multipliers `\beta_1, \dots, \beta_n`, the
corresponding GRS code of dimension `k` is the set:

.. math::

    \{ (\beta_1 f(\alpha_1), \ldots, \beta_n f(\alpha_n)  \mid  f \in F[x], \deg f < k \}

This file contains the following elements:

    - :class:`GeneralizedReedSolomonCode`, the class for GRS codes
    - :class:`GRSEvaluationVectorEncoder`, an encoder with a vectorial message space
    - :class:`GRSEvaluationPolynomialEncoder`, an encoder with a polynomial message space
    - :class:`GRSBerlekampWelchDecoder`, a decoder which corrects errors using Berlekamp-Welch algorithm
    - :class:`GRSGaoDecoder`, a decoder which corrects errors using Gao algorithm
    - :class:`GRSErrorErasureDecoder`, a decoder which corrects both errors and erasures
    - :class:`GRSKeyEquationSyndromeDecoder`, a decoder which corrects errors using the key equation on syndrome polynomials
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix, diagonal_matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.categories.cartesian_product import cartesian_product
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method
from copy import copy
from .linear_code import AbstractLinearCode
from .encoder import Encoder
from .decoder import Decoder, DecodingError
from sage.rings.arith import xgcd
from sage.misc.misc_c import prod
from sage.functions.other import binomial, floor, sqrt
from sage.calculus.var import var
from sage.misc.functional import symbolic_sum
from sage.rings.integer_ring import ZZ

class GeneralizedReedSolomonCode(AbstractLinearCode):
    r"""
    Representation of a Generalized Reed-Solomon code.

    INPUT:

    - ``evaluation_points`` -- A list of distinct elements of some finite field `F`.

    - ``dimension`` -- The dimension of the resulting code.

    - ``column_multipliers`` -- (default: ``None``) List of non-zero elements of `F`.
      All column multipliers are set to 1 if default value is kept.

    EXAMPLES:

    A Reed-Solomon code can be constructed by taking all non-zero elements of
    the field as evaluation points, and specifying no column multipliers::

        sage: F = GF(7)
        sage: evalpts = [F(i) for i in range(1,7)]
        sage: C = codes.GeneralizedReedSolomonCode(evalpts,3)
        sage: C
        [6, 3, 4] Generalized Reed-Solomon Code over Finite Field of size 7

    More generally, the following is a GRS code where the evaluation points are
    a subset of the field and includes zero::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: C
        [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    It is also possible to specify the column multipliers::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: colmults = F.list()[1:n+1]
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, colmults)
        sage: C
        [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, evaluation_points, dimension, column_multipliers=None):
        r"""
        TESTS:

        If the evaluation points are not from a finite field, it raises an error::

            sage: C = codes.GeneralizedReedSolomonCode([1,2,3], 1)
            Traceback (most recent call last):
            ...
            ValueError: Evaluation points must be in a finite field (and Integer Ring is not one)

        If the evaluation points are not from the same finite field, it raises an error::

            sage: F2, F3 = GF(2) , GF(3)
            sage: C = codes.GeneralizedReedSolomonCode([F2.zero(),F2.one(),F3(2)], 1)
            Traceback (most recent call last):
            ...
            ValueError: Failed converting all evaluation points to the same field (unable to find a common ring for all elements)

        If the column multipliers cannot be converted into the finite are not from a finite field, or cannot be not in the same
        finite field as the evaluation points, it raises an error::

            sage: F = GF(59)
            sage: F2 = GF(61)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, [.3]*n )
            Traceback (most recent call last):
            ...
            ValueError: Failed converting all evaluation points and column multipliers to the same field (unable to find a common ring for all elements)

            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F2.list()[1:n+1])
            Traceback (most recent call last):
            ...
            ValueError: Failed converting all evaluation points and column multipliers to the same field (unable to find a common ring for all elements)

        The number of column multipliers is checked as well::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[1:n])
            Traceback (most recent call last):
            ...
            ValueError: There must be the same number of evaluation points as column multipliers

        It is not allowed to have 0 as a column multiplier::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[:n])
            Traceback (most recent call last):
            ...
            ValueError: All column multipliers must be non-zero

        And all the evaluation points must be different. Note that they should
        be different after converting into the same field::

            sage: F = GF(5)
            sage: C = codes.GeneralizedReedSolomonCode([ F(0), 1, 2, 3, 5 ], 3)
            Traceback (most recent call last):
            ...
            ValueError: All evaluation points must be different

        The dimension is not allowed to exceed the length::

            sage: F = GF(59)
            sage: n, k = 40, 100
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            Traceback (most recent call last):
            ...
            ValueError: The dimension must be a positive integer at most the length of the code.
        """
        if column_multipliers:
            if len(evaluation_points) != len(column_multipliers):
                raise ValueError("There must be the same number of evaluation points as column multipliers");
            try:
                common_points = vector(list(evaluation_points) + list(column_multipliers))
                F = common_points.base_ring()
                self._evaluation_points = common_points[:len(evaluation_points)]
                self._column_multipliers = common_points[len(evaluation_points):]
            except (TypeError, ValueError) as e:
                raise ValueError("Failed converting all evaluation points and column multipliers to the same field (%s)" % e.message)
        else:
            try:
                self._evaluation_points = vector(evaluation_points)
                F = self._evaluation_points.base_ring()
                self._column_multipliers = vector(F, [F.one()] * len(self._evaluation_points))
            except (TypeError, ValueError) as e:
                raise ValueError("Failed converting all evaluation points to the same field (%s)" % e.message)

        if F.is_finite() == False or F.is_field() == False:
            raise ValueError("Evaluation points must be in a finite field (and %s is not one)" % F)
        super(GeneralizedReedSolomonCode, self).__init__(F,
                len(self._evaluation_points), "EvaluationVector", "Gao")

        if dimension not in ZZ or dimension > self._length or dimension < 1:
            raise ValueError("The dimension must be a positive integer at most the length of the code.")
        self._dimension = dimension

        if F.zero() in self._column_multipliers:
            raise ValueError("All column multipliers must be non-zero")
        if len(self._evaluation_points) != len(set(self._evaluation_points)):
            raise ValueError("All evaluation points must be different")

    def __eq__(self, other):
        r"""
        Tests equality between Generalized Reed-Solomon Code objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C1.__eq__(C2)
            True
        """
        return isinstance(other, GeneralizedReedSolomonCode) \
                and self.base_field() == other.base_field() \
                and self.length() == other.length() \
                and self.dimension() == other.dimension() \
                and self.evaluation_points() == other.evaluation_points() \
                and self.column_multipliers() == other.column_multipliers()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C
            [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "[%s, %s, %s] Generalized Reed-Solomon Code over %s"\
                % (self.length(), self.dimension(),\
                self.minimum_distance(), self.base_field())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: latex(C)
            [40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "[%s, %s, %s] \\textnormal{ Generalized Reed-Solomon Code over } %s"\
                % (self.length(), self.dimension() ,self.minimum_distance(),
                self.base_field()._latex_())

    def minimum_distance(self):
        r"""
        Returns the minimum distance of ``self``. Since a GRS code is always
        Maximum-Distance-Separable (MDS), this returns ``C.length() -
        C.dimension() + 1``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.minimum_distance()
            29
        """
        return self.length() - self.dimension() + 1

    def evaluation_points(self):
        r"""
        Returns the evaluation points of ``self`` as a vector.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.evaluation_points()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        """
        return self._evaluation_points

    def column_multipliers(self):
        r"""
        Returns the column multipliers of ``self`` as a vector.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.column_multipliers()
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        """
        return self._column_multipliers

    @cached_method
    def multipliers_product(self):
        r"""
        Returns the component-wise product of the column multipliers of ``self``
        with the column multipliers of the dual GRS code.

        This is a simple Cramer's rule-like expression on the evaluation points
        of ``self``. Recall that the column multipliers of the dual GRS code is
        also the column multipliers of the parity check matrix of ``self``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.multipliers_product()
            [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        a = self.evaluation_points()
        one = self.base_ring().one()
        return [ one/prod([ a[i] - a[h] for h in range(0, len(a)) if h != i ])
                    for i in range(0,len(a)) ]

    @cached_method
    def parity_column_multipliers(self):
        r"""
        Returns the list of column multipliers of ``self``'s parity check matrix.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.parity_column_multipliers()
            [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        n = self.length()
        col_mults = self.column_multipliers()
        etas = self.multipliers_product()
        return [ etas[i]/col_mults[i] for i in range(n) ]

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.parity_check_matrix()
            [10  9  8  7  6  5  4  3  2  1]
            [ 0  9  5 10  2  3  2 10  5  9]
            [ 0  9 10  8  8  4  1  4  7  4]
            [ 0  9  9  2 10  9  6  6  1  3]
            [ 0  9  7  6  7  1  3  9  8  5]
        """
        return self.dual_code().generator_matrix()

    @cached_method
    def dual_code(self):
        r"""
        Returns the dual code of ``self``, which is also a GRS code.

        EXAMPLES::

            sage: F =  GF(59)
            sage: colmults = [ F.random_element() for i in range(40) ]
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:40], 12, colmults)
            sage: Cd = C.dual_code(); Cd
            [40, 28, 13] Generalized Reed-Solomon Code over Finite Field of size 59

        The dual code of the dual code is the original code::

            sage: C == Cd.dual_code()
            True
        """
        col_mults = self.parity_column_multipliers()
        return GeneralizedReedSolomonCode(self.evaluation_points(), self.length() - self.dimension(), col_mults)

    def covering_radius(self):
        r"""
        Returns the covering radius of ``self``.

        The covering radius of a linear code `C` is the smallest
        number `r` s.t. any element of the ambient space of `C` is at most at
        distance `r` to `C`.

        As GRS codes are Maximum Distance Separable codes (MDS), their covering
        radius is always `d-1`, where `d` is the minimum distance. This is
        opposed to random linear codes where the covering radius is
        computationally hard to determine.

        EXAMPLES::

            sage: F = GF(2^8, 'a')
            sage: n, k = 256, 100
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.covering_radius()
            156
        """
        return self.length() - self.dimension()

    @cached_method
    def weight_distribution(self):
        r"""
        Returns the list whose `i`'th entry is the number of words of weight `i`
        in ``self``.

        Computing the weight distribution for a GRS code is very fast. Note that
        for random linear codes, it is computationally hard.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.weight_distribution()
            [1, 0, 0, 0, 0, 0, 2100, 6000, 29250, 61500, 62200]
        """
        d = self.minimum_distance()
        n = self.length()
        q = self.base_ring().order()
        s = var('s')
        wd = [1] + [0] * (d - 1)
        for i in range(d, n+1):
            tmp = binomial(n, i) * (q - 1)
            wd.append(tmp * symbolic_sum(binomial(i-1, s) * (-1) ** s * q ** (i - d - s), s, 0, i-d))
        return wd

    def weight_enumerator(self):
        r"""
        Returns the polynomial whose coefficient to `x^i` is the number of codewords of weight `i` in ``self``.

        Computing the weight enumerator for a GRS code is very fast. Note that
        for random linear codes, it is computationally hard.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.weight_enumerator()
            62200*x^10 + 61500*x^9 + 29250*x^8 + 6000*x^7 + 2100*x^6 + 1
        """
        PolRing = ZZ['x']
        x = PolRing.gen()
        s = var('s')
        wd = self.weight_distribution()
        d = self.minimum_distance()
        n = self.length()
        w_en = PolRing(1)
        for i in range(n + 1 - d):
            w_en += wd[i + d] * x ** (i + d)
        return w_en

    def _punctured_form(self, points):
        r"""
        Returns a representation of self as a
        :class:`GeneralizedReedSolomonCode` punctured in ``points``.

        INPUT:

        - ``points`` -- a set of positions where to puncture ``self``

        EXAMPLES::

            sage: C_grs = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C_grs._punctured_form({4, 3})
            [38, 12, 27] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        if not isinstance(points, (Integer, int, set)):
            raise TypeError("points must be either a Sage Integer, a Python int, or a set")
        alphas = list(self.evaluation_points())
        col_mults = list(self.column_multipliers())
        n = self.length()
        punctured_alphas = []
        punctured_col_mults = []
        punctured_alphas = [alphas[i] for i in range(n) if i not in points]
        punctured_col_mults = [col_mults[i] for i in range(n) if i not in points]
        G = self.generator_matrix()
        G = G.delete_columns(list(points))
        dimension = G.rank()
        return GeneralizedReedSolomonCode(punctured_alphas, dimension, punctured_col_mults)

    def decode_to_message(self, r):
        r"""
        Decodes``r`` to an element in message space of ``self``

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be unencoded as is. In that case,
            if ``r`` is not a codeword, the output is unspecified.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 2))
            sage: C.decode_to_message(r)
            (3, 6, 6, 3, 1)
        """
        if self.length() == self.dimension():
            return self.encoder().unencode_nocheck(r)
        return vector(self.decoder().decode_to_message(r))






####################### encoders ###############################


####################### encoders ###############################


class GRSEvaluationVectorEncoder(Encoder):
    r"""
    Encoder for Generalized Reed-Solomon codes which encodes vectors into codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
        sage: E
        Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    Actually, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("EvaluationVector")
        sage: E
        Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E
            Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        super(GRSEvaluationVectorEncoder, self).__init__(code)
        self._R = code.base_field()['x']

    def __eq__(self, other):
        r"""
        Tests equality between GRSEvaluationVectorEncoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: D2 = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: D1.__eq__(D2)
            True
            sage: D1 is D2
            False
        """
        return isinstance(other, GRSEvaluationVectorEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E
            Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Evaluation vector-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: latex(E)
            \textnormal{Evaluation vector-style encoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Evaluation vector-style encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E.generator_matrix()
            [1 1 1 1 1 1 1 1 1 1]
            [0 1 2 3 4 5 6 7 8 9]
            [0 1 4 9 5 3 3 5 9 4]
            [0 1 8 5 9 4 7 2 6 3]
            [0 1 5 4 3 9 9 3 4 5]
        """
        C = self.code()
        alphas = C.evaluation_points()
        col_mults = C.column_multipliers()
        return matrix(C.base_field(), C.dimension(), C.length(), lambda i,j: col_mults[j] * alphas[j]**i)








class GRSEvaluationPolynomialEncoder(Encoder):
    r"""
    Encoder for Generalized Reed-Solomon codes which uses evaluation of
    polynomials to obtain codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
        sage: E
        Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        sage: E.message_space()
        Univariate Polynomial Ring in x over Finite Field of size 59

    Actually, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("EvaluationPolynomial")
        sage: E
        Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: E
            Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        super(GRSEvaluationPolynomialEncoder, self).__init__(code)
        self._R = code.base_field()['x']

    def __eq__(self, other):
        r"""
        Tests equality between GRSEvaluationPolynomialEncoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: D2 = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: D1 is D2
            False
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSEvaluationPolynomialEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: E
            Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Evaluation polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: latex(E)
            \textnormal{Evaluation polynomial-style encoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Evaluation polynomial-style encoder for }%s" % self.code()._latex_()

    def encode(self, p):
        r"""
        Transforms the polynomial ``p`` into a codeword of :meth:`code`.

        INPUT:

        - ``p`` -- A polynomial from the message space of ``self`` of degree
          less than ``self.code().dimension()``.

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: Fx.<x> = F[]
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: p = x^2 + 3*x + 10
            sage: c = E.encode(p); c
            (10, 3, 9, 6, 5, 6, 9, 3, 10, 8)
            sage: c in C
            True

        If a polynomial of too high degree is given, an error is raised::

            sage: p = x^10
            sage: E.encode(p)
            Traceback (most recent call last):
            ...
            ValueError: The polynomial to encode must have degree at most 4

        If ``p`` is not an element of the proper polynomial ring, an error is raised::

            sage: Qy.<y> = QQ[]
            sage: p = y^2 + 1
            sage: E.encode(p)
            Traceback (most recent call last):
            ...
            ValueError: The value to encode must be in Univariate Polynomial Ring in x over Finite Field of size 11
        """
        M = self.message_space()
        if p not in M:
            raise ValueError("The value to encode must be in %s" % M)
        C = self.code()
        if p.degree() >= C.dimension():
            raise ValueError("The polynomial to encode must have degree at most %s" % (C.dimension() - 1))
        alphas    = C.evaluation_points()
        col_mults = C.column_multipliers()
        c = vector(C.base_ring(), [col_mults[i]*p(alphas[i]) for i in range(C.length())])
        return c

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to the codeword ``c``.

        Use this method with caution: it does not check if ``c``
        belongs to the code, and if this is not the case, the output is
        unspecified. Instead, use :meth:`unencode`.

        INPUT:

        - ``c`` -- A codeword of :meth:`code`.

        OUTPUT:

        - An polynomial of degree less than ``self.code().dimension()``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: c = vector(F, (10, 3, 9, 6, 5, 6, 9, 3, 10, 8))
            sage: c in C
            True
            sage: p = E.unencode_nocheck(c); p
            x^2 + 3*x + 10
            sage: E.encode(p) == c
            True

        Note that no error is thrown if ``c`` is not a codeword, and that the
        result is undefined::

            sage: c = vector(F, (11, 3, 9, 6, 5, 6, 9, 3, 10, 8))
            sage: c in C
            False
            sage: p = E.unencode_nocheck(c); p
            6*x^4 + 6*x^3 + 2*x^2
            sage: E.encode(p) == c
            False

        """
        C = self.code()
        alphas    = C.evaluation_points()
        col_mults = C.column_multipliers()

        c = [c[i]/col_mults[i] for i in range(C.length())]
        points = [(alphas[i], c[i]) for i in range(C.dimension())]

        Pc = self._R.lagrange_polynomial(points)
        return Pc

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: E.message_space()
            Univariate Polynomial Ring in x over Finite Field of size 11
        """
        return self._R









####################### decoders ###############################

class GRSBerlekampWelchDecoder(Decoder):
    r"""
    Decoder for Generalized Reed-Solomon codes which uses Berlekamp-Welch
    decoding algorithm to correct errors in codewords.

    This algorithm recovers the error locator polynomial by solving a linear system.
    See [HJ04]_ pp. 51-52 for details.

    REFERENCES:

    .. [HJ04] Tom Hoeholdt and Joern Justesen, A Course In Error-Correcting Codes,
       EMS, 2004

    INPUT:

    - ``code`` -- A code associated to this decoder

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
        sage: D
        Berlekamp-Welch decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    Actually, we can construct the decoder from ``C`` directly::

        sage: D = C.decoder("BerlekampWelch")
        sage: D
        Berlekamp-Welch decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        TESTS:

        If ``code`` is not a GRS code, an error is raised::

            sage: C  = codes.RandomLinearCode(10, 4, GF(11))
            sage: codes.decoders.GRSBerlekampWelchDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a generalized Reed-Solomon code
        """
        if not isinstance(code, GeneralizedReedSolomonCode):
            raise ValueError("code has to be a generalized Reed-Solomon code")
        super(GRSBerlekampWelchDecoder, self).__init__(code, code.ambient_space(),
            "EvaluationPolynomial")

    def __eq__(self, other):
        r"""
        Tests equality between GRSBerlekampWelchDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D2 = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D1.__eq__(D2)
            True
            sage: D1 is D2
            False
        """
        return isinstance(other, GRSBerlekampWelchDecoder) \
                and self.code() == other.code()\
                and self.input_space() == other.input_space()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D
            Berlekamp-Welch decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Berlekamp-Welch decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: latex(D)
            \textnormal{Berlekamp Welch decoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Berlekamp Welch decoder for }%s"\
                % self.code()._latex_()

    def _decode_to_code_and_message(self, r):
        r"""
        Decodes ``r`` to an element in message space of ``self`` and its
        representation in the ambient space of the code associated to ``self``.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - ``(c, f)`` -- ``c`` is the representation of ``r`` decoded in the ambient
          space of the associated code of ``self``, ``f`` its representation in
          the message space of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: c_dec, f_dec = D._decode_to_code_and_message(y)
            sage: f_dec == D.connected_encoder().unencode(c)
            True
            sage: c_dec == c
            True
        """
        C = self.code()
        if r not in C.ambient_space():
            raise ValueError("The word to decode has to be in the ambient space of the code")
        n, k = C.length(), C.dimension()
        if n == k:
            return r, self.connected_encoder().unencode_nocheck(r)
        if r in C:
            return r, self.connected_encoder().unencode_nocheck(r)
        col_mults = C.column_multipliers()

        r_list = copy(r)
        r_list = [r[i]/col_mults[i] for i in range(0, C.length())]

        t  = (C.minimum_distance()-1) // 2
        l0 = n-1-t
        l1 = n-1-t-(k-1)
        S  = matrix(C.base_field(), n, l0+l1+2, lambda i,j :
                (C.evaluation_points()[i])**j if j<(l0+1)
                else r_list[i]*(C.evaluation_points()[i])**(j-(l0+1)))
        S  = S.right_kernel()
        S  = S.basis_matrix().row(0)
        R = C.base_field()['x']

        Q0 = R(S.list_from_positions(xrange(0, l0+1)))
        Q1 = R(S.list_from_positions(xrange(l0+1, l0+l1+2)))

        f, rem = (-Q0).quo_rem(Q1)
        if not rem.is_zero():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        if f not in R:
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        c = self.connected_encoder().encode(f)
        if (c - r).hamming_weight() > self.decoding_radius():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        return c, f

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to an element in message space of ``self``.

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be unencoded as is. In that case,
            if ``r`` is not a codeword, the output is unspecified.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: D.connected_encoder().unencode(c) == D.decode_to_message(y)
            True

        TESTS:

        If one tries to decode a word which is too far from any codeword, an exception is raised::

            sage: e = vector(F,[0, 0, 54, 23, 1, 0, 0, 0, 53, 21, 0, 0, 0, 34, 6, 11, 0, 0, 16, 0, 0, 0, 9, 0, 10, 27, 35, 0, 0, 0, 0, 46, 0, 0, 0, 0, 0, 0, 44, 0]); e.hamming_weight()
            15
            sage: D.decode_to_message(c + e)
            Traceback (most recent call last):
            ...
            DecodingError: Decoding failed because the number of errors exceeded the decoding radius

        If one tries to decode something which is not in the ambient space of the code,
        an exception is raised::

            sage: D.decode_to_message(42)
            Traceback (most recent call last):
            ...
            ValueError: The word to decode has to be in the ambient space of the code

        The bug detailed in :trac:`20340` has been fixed::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: c = C.random_element()
            sage: D = C.decoder("BerlekampWelch")
            sage: E = D.connected_encoder()
            sage: m = E.message_space().random_element()
            sage: c = E.encode(m)
            sage: D.decode_to_message(c) == m
            True
        """
        return self._decode_to_code_and_message(r)[1]

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``r`` and returns a codeword.

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be returned as is.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a vector of ``self.code()``

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: c == D.decode_to_code(y)
            True

        TESTS:

        If one tries to decode a word which is too far from any codeword, an exception is raised::

            sage: e = vector(F,[0, 0, 54, 23, 1, 0, 0, 0, 53, 21, 0, 0, 0, 34, 6, 11, 0, 0, 16, 0, 0, 0, 9, 0, 10, 27, 35, 0, 0, 0, 0, 46, 0, 0, 0, 0, 0, 0, 44, 0]); e.hamming_weight()
            15
            sage: D.decode_to_code(c + e)
            Traceback (most recent call last):
            ...
            DecodingError: Decoding failed because the number of errors exceeded the decoding radius

        If one tries to decode something which is not in the ambient space of the code,
        an exception is raised::

            sage: D.decode_to_code(42)
            Traceback (most recent call last):
            ...
            ValueError: The word to decode has to be in the ambient space of the code

        The bug detailed in :trac:`20340` has been fixed::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: c = C.random_element()
            sage: D = C.decoder("BerlekampWelch")
            sage: D.decode_to_code(c) == c
            True
        """
        return self._decode_to_code_and_message(r)[0]

    def decoding_radius(self):
        r"""
        Returns maximal number of errors that ``self`` can decode.

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D.decoding_radius()
            14
        """
        return (self.code().minimum_distance()-1)//2











class GRSGaoDecoder(Decoder):
    r"""
    Decoder for Generalized Reed-Solomon codes which uses Gao
    decoding algorithm to correct errors in codewords.

    Gao decoding algorithm uses early terminated extended Euclidean algorithm
    to find the error locator polynomial. See [G02]_ for details.

    REFERENCES:

    .. [G02] Shuhong Gao, A new algorithm for decoding Reed-Solomon Codes, January 31, 2002

    INPUT:

    - ``code`` -- The associated code of this decoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: D = codes.decoders.GRSGaoDecoder(C)
        sage: D
        Gao decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    Actually, we can construct the decoder from ``C`` directly::

        sage: D = C.decoder("Gao")
        sage: D
        Gao decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        TESTS:

        If ``code`` is not a GRS code, an error is raised::

            sage: C  = codes.RandomLinearCode(10, 4, GF(11))
            sage: codes.decoders.GRSGaoDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a generalized Reed-Solomon code
        """
        if not isinstance(code, GeneralizedReedSolomonCode):
            raise ValueError("code has to be a generalized Reed-Solomon code")
        super(GRSGaoDecoder, self).__init__(code, code.ambient_space(),
                "EvaluationPolynomial")

    def __eq__(self, other):
        r"""
        Tests equality of GRSGaoDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.decoders.GRSGaoDecoder(C)
            sage: D2 = codes.decoders.GRSGaoDecoder(C)
            sage: D1.__eq__(D2)
            True
            sage: D1 is D2
            False
        """
        return isinstance(other, GRSGaoDecoder) \
                and self.code() == other.code()\
                and self.input_space() == other.input_space()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: D
            Gao decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Gao decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: latex(D)
            \textnormal{Gao decoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Gao decoder for }%s" % self.code()._latex_()

    @cached_method
    def _polynomial_vanishing_at_alphas(self, PolRing):
        r"""
        Return the unique minimal-degree polynomial vanishing at all the evaluation points.

        INPUT:

        - ``PolRing`` -- polynomial ring of the output

        OUTPUT:

        - a polynomial over ``PolRing``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: P = PolynomialRing(F,'x')
            sage: D._polynomial_vanishing_at_alphas(P)
            x^10 + 10*x^9 + x^8 + 10*x^7 + x^6 + 10*x^5 + x^4 + 10*x^3 + x^2 + 10*x
        """
        alphas = self.code().evaluation_points()
        G = PolRing.one()
        x = PolRing.gen()
        for i in range(0, self.code().length()):
            G = G*(x-self.code().evaluation_points()[i])
        return G

    def _partial_xgcd(self, a, b, PolRing):
        r"""
        Performs an Euclidean algorithm on ``a`` and ``b`` until a remainder
        has degree less than `\frac{n+k}{2}`, `n` being the dimension of the
        code, `k` its dimension, and returns `(r, s)` such that in the step
        just before termination, `r = a\times s + b\times t`.

        INPUT:

        - ``a, b`` -- polynomials over ``PolRing``

        - ``PolRing`` -- polynomial ring of the output

        OUTPUT:

        - a tuple of polynomials

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: P = PolynomialRing(F,'x')
            sage: x = P.parameter()
            sage: a = 5*x^2 + 9*x + 8
            sage: b = 10*x^2 + 3*x + 5
            sage: D._partial_xgcd(a, b, P)
            (10*x^2 + 3*x + 5, 1)
        """
        stop = floor(self.code().dimension() + self.code().length()) // 2
        s = PolRing.one()
        prev_s = PolRing.zero()

        r = b
        prev_r = a
        while(r.degree() >= stop):
            q = prev_r.quo_rem(r)[0]
            (prev_r, r) = (r, prev_r - q * r)
            (prev_s, s) = (s, prev_s - q * s)

        return (r, s)

    def _decode_to_code_and_message(self, r):
        r"""
        Decodes ``r`` to an element in message space of ``self`` and its
        representation in the ambient space of the code associated to ``self``.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - ``(c, h)`` -- ``c`` is the representation of ``r`` decoded in the ambient
          space of the associated code of ``self``, ``h`` its representation in
          the message space of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: c_dec, h_dec = D._decode_to_code_and_message(y)
            sage: h_dec == D.connected_encoder().unencode(c)
            True
            sage: c_dec == c
            True
        """
        C = self.code()
        if r not in C.ambient_space():
            raise ValueError("The word to decode has to be in the ambient space of the code")
        alphas = C.evaluation_points()
        col_mults = C.column_multipliers()
        PolRing = C.base_field()['x']
        G = self._polynomial_vanishing_at_alphas(PolRing)
        n = C.length()

        if n == C.dimension() or r in C:
            return r, self.connected_encoder().unencode_nocheck(r)

        points = [(alphas[i], r[i]/col_mults[i]) for i in
                range(0, n)]
        R = PolRing.lagrange_polynomial(points)

        (Q1, Q0) = self._partial_xgcd(G, R, PolRing)

        h, rem = Q1.quo_rem(Q0)
        if not rem.is_zero():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        if h not in PolRing:
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        c = self.connected_encoder().encode(h)
        if (c - r).hamming_weight() > self.decoding_radius():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        return c, h

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to an element in message space of ``self``.

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be unencoded as is. In that case,
            if ``r`` is not a codeword, the output is unspecified.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: D.connected_encoder().unencode(c) == D.decode_to_message(y)
            True

        TESTS:

        If one tries to decode a word which is too far from any codeword, an exception is raised::

            sage: e = vector(F,[0, 0, 54, 23, 1, 0, 0, 0, 53, 21, 0, 0, 0, 34, 6, 11, 0, 0, 16, 0, 0, 0, 9, 0, 10, 27, 35, 0, 0, 0, 0, 46, 0, 0, 0, 0, 0, 0, 44, 0]); e.hamming_weight()
            15
            sage: D.decode_to_message(c + e)
            Traceback (most recent call last):
            ...
            DecodingError: Decoding failed because the number of errors exceeded the decoding radius

        If one tries to decode something which is not in the ambient space of the code,
        an exception is raised::

            sage: D.decode_to_message(42)
            Traceback (most recent call last):
            ...
            ValueError: The word to decode has to be in the ambient space of the code

        The bug detailed in :trac:`20340` has been fixed::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: c = C.random_element()
            sage: D = C.decoder("Gao")
            sage: E = D.connected_encoder()
            sage: m = E.message_space().random_element()
            sage: c = E.encode(m)
            sage: D.decode_to_message(c) == m
            True
        """
        return self._decode_to_code_and_message(r)[1]

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``r`` and returns a codeword.

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be returned as is.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a vector of ``self.code()``

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: c == D.decode_to_code(y)
            True

        TESTS:


        If one tries to decode a word which is too far from any codeword, an exception is raised::

            sage: e = vector(F,[0, 0, 54, 23, 1, 0, 0, 0, 53, 21, 0, 0, 0, 34, 6, 11, 0, 0, 16, 0, 0, 0, 9, 0, 10, 27, 35, 0, 0, 0, 0, 46, 0, 0, 0, 0, 0, 0, 44, 0]); e.hamming_weight()
            15
            sage: D.decode_to_code(c + e)
            Traceback (most recent call last):
            ...
            DecodingError: Decoding failed because the number of errors exceeded the decoding radius

        If one tries to decode something which is not in the ambient space of the code,
        an exception is raised::

            sage: D.decode_to_code(42)
            Traceback (most recent call last):
            ...
            ValueError: The word to decode has to be in the ambient space of the code

        The bug detailed in :trac:`20340` has been fixed::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: c = C.random_element()
            sage: D = C.decoder("Gao")
            sage: c = C.random_element()
            sage: D.decode_to_code(c) == c
            True
        """
        return self._decode_to_code_and_message(r)[0]

    def decoding_radius(self):
        r"""
        Return maximal number of errors that ``self`` can decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: D.decoding_radius()
            14
        """
        return (self.code().minimum_distance()-1)//2










class GRSErrorErasureDecoder(Decoder):
    r"""
    Decoder for Generalized Reed-Solomon codes which is able to correct both errors
    and erasures in codewords.

    INPUT:

    - ``code`` -- The associated code of this decoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: D = codes.decoders.GRSErrorErasureDecoder(C)
        sage: D
        Error-Erasure decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    Actually, we can construct the decoder from ``C`` directly::

        sage: D = C.decoder("ErrorErasure")
        sage: D
        Error-Erasure decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        TESTS:

        If ``code`` is not a GRS code, an error is raised::

            sage: C  = codes.RandomLinearCode(10, 4, GF(11))
            sage: codes.decoders.GRSErrorErasureDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a generalized Reed-Solomon code
        """
        if not isinstance(code, GeneralizedReedSolomonCode):
            raise ValueError("code has to be a generalized Reed-Solomon code")
        input_space = cartesian_product([code.ambient_space(),
                VectorSpace(GF(2), code.ambient_space().dimension())])
        super(GRSErrorErasureDecoder, self).__init__(code, input_space, "EvaluationVector")

    def __eq__(self, other):
        r"""
        Test equality of GRSErrorErasureDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D2 = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D1.__eq__(D2)
            True
            sage: D1 is D2
            False
        """
        return isinstance(other, GRSErrorErasureDecoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D
            Error-Erasure decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Error-Erasure decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: latex(D)
            \textnormal{Error-Erasure decoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Error-Erasure decoder for }%s"\
                % self.code()._latex_()

    def decode_to_message(self, word_and_erasure_vector):
        r"""
        Decode ``word_and_erasure_vector`` to an element in message space
        of ``self``

        INPUT:

        - word_and_erasure_vector -- a tuple whose:
          - first element is an element of the ambient space of the code
          - second element is a vector over GF(2) whose length is the same as the code's

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be unencoded as is.
            If the number of erasures is exactly `n - k`, where `n` is the
            length of the code associated to ``self`` and `k` its dimension,
            ``r`` will be returned as is.
            In either case, if ``r`` is not a codeword,
            the output is unspecified.

        INPUT:

        - ``word_and_erasure_vector`` -- a pair of vectors, where
          first element is a codeword of ``self`` and second element
          is a vector of GF(2) containing erasure positions

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: c = C.random_element()
            sage: n_era = randint(0, C.minimum_distance() - 2)
            sage: Chan = channels.ErrorErasureChannel(C.ambient_space(), D.decoding_radius(n_era), n_era)
            sage: y = Chan(c)
            sage: D.connected_encoder().unencode(c) == D.decode_to_message(y)
            True

        TESTS:

        If one tries to decode a word with too many erasures, it returns
        an exception::

            sage: Chan = channels.ErrorErasureChannel(C.ambient_space(), 0, C.minimum_distance() + 1)
            sage: y = Chan(c)
            sage: D.decode_to_message(y)
            Traceback (most recent call last):
            ...
            DecodingError: Too many erasures in the received word

        If one tries to decode something which is not in the ambient space of the code,
        an exception is raised::

            sage: D.decode_to_message((42, random_vector(GF(2), C.length())))
            Traceback (most recent call last):
            ...
            ValueError: The word to decode has to be in the ambient space of the code

        If one tries to pass an erasure_vector which is not a vector over GF(2) of the same length as code's,
        an exception is raised::

            sage: D.decode_to_message((C.random_element(), 42))
            Traceback (most recent call last):
            ...
            ValueError: The erasure vector has to be a vector over GF(2) of the same length as the code
        """
        C = self.code()
        word, erasure_vector = word_and_erasure_vector
        n, k = C.length(), C.dimension()
        if word not in C.ambient_space():
            raise ValueError("The word to decode has to be in the ambient space of the code")
        if not erasure_vector in VectorSpace(GF(2), n):
            raise ValueError("The erasure vector has to be a vector over GF(2) of the same length as the code")
        if erasure_vector.hamming_weight() >= self.code().minimum_distance():
            raise DecodingError("Too many erasures in the received word")

        punctured_word = vector(self.code().base_ring(), [word[i] for i in
            range(len(word)) if erasure_vector[i]!=1])
        C1_length = len(punctured_word)
        if C1_length == k:
            return self.connected_encoder().unencode_nocheck(word)
        C1_evaluation_points = [self.code().evaluation_points()[i] for i in
                range(n) if erasure_vector[i]!=1]
        C1_column_multipliers = [self.code().column_multipliers()[i] for i in
                range(n) if erasure_vector[i]!=1]
        C1 = GeneralizedReedSolomonCode(C1_evaluation_points, k,
                C1_column_multipliers)
        return C1.decode_to_message(punctured_word)

    def decoding_radius(self, number_erasures):
        r"""
        Return maximal number of errors that ``self`` can decode according
        to how many erasures it receives

        INPUT:

        - ``number_erasures`` -- the number of erasures when we try to decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D.decoding_radius(5)
            11

        If we receive too many erasures, it returns an exception as codeword will
        be impossible to decode::

            sage: D.decoding_radius(30)
            Traceback (most recent call last):
            ...
            ValueError: The number of erasures exceed decoding capability
        """
        diff = self.code().minimum_distance() - 1 - number_erasures
        if diff <= 0:
            raise ValueError("The number of erasures exceed decoding capability")
        else :
            return diff // 2








class GRSKeyEquationSyndromeDecoder(Decoder):
    r"""
    Decoder for Generalized Reed-Solomon codes which uses a
    Key equation decoding based on the syndrome polynomial to
    correct errors in codewords.

    This algorithm uses early terminated extended euclidean algorithm
    to solve the key equations, as described in [R06]_, pp. 183-195.

    REFERENCES:

        .. [R06] Ron Roth, Introduction to Coding Theory, Cambridge University Press, 2006

    INPUT:

    - ``code`` -- The associated code of this decoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
        sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
        sage: D
        Key equation decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    Actually, we can construct the decoder from ``C`` directly::

        sage: D = C.decoder("KeyEquationSyndrome")
        sage: D
        Key equation decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        TESTS::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Impossible to use this decoder over a GRS code which contains 0 amongst its evaluation points

        If ``code`` is not a GRS code, an error is raised::

            sage: C  = codes.RandomLinearCode(10, 4, GF(11))
            sage: codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a generalized Reed-Solomon code
        """
        if not isinstance(code, GeneralizedReedSolomonCode):
            raise ValueError("code has to be a generalized Reed-Solomon code")
        if code.base_field().zero() in code.evaluation_points():
            raise ValueError("Impossible to use this decoder over a GRS code which contains 0 amongst its evaluation points")
        super(GRSKeyEquationSyndromeDecoder, self).__init__(code, code.ambient_space(),
                "EvaluationVector")

    def __eq__(self, other):
        r"""
        Tests equality of GRSKeyEquationSyndromeDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D1 = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D2 = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D1.__eq__(D2)
            True
            sage: D1 is D2
            False
        """
        return isinstance(other, GRSKeyEquationSyndromeDecoder) \
                and self.code() == other.code()\
                and self.input_space() == other.input_space()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D
            Key equation decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Key equation decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: latex(D)
            \textnormal{Key equation decoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Key equation decoder for }%s" % self.code()._latex_()

    def _partial_xgcd(self, a, b, PolRing):
        r"""
        Performs an Euclidean algorithm on ``a`` and ``b`` until a remainder
        has degree less than `\frac{n+k}{2}`, `n` being the dimension of the
        code, `k` its dimension, and returns `(r, t)` such that in the step
        just before termination, `r = a\times s + b\times t`.

        INPUT:

        - ``a, b`` -- polynomials over ``PolRing``

        - ``PolRing`` -- polynomial ring of the output

        OUTPUT:

        - a tuple of polynomials

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: P = PolynomialRing(F,'x')
            sage: x = P.parameter()
            sage: a = 5*x^2 + 9*x + 8
            sage: b = 10*x^2 + 3*x + 5
            sage: D._partial_xgcd(a, b, P)
            (5, 8*x + 10)
        """
        prev_t = PolRing.zero()
        t = PolRing.one()

        prev_r = a
        r = b

        while(r.degree() >= t.degree()):
            q = prev_r.quo_rem(r)[0]
            prev_r, r = r, prev_r - q * r
            prev_t, t = t, prev_t - q * t

        return (r, t)

    def _syndrome(self, r):
        r"""
        Returns the coefficients of the syndrome polynomial of ``r``.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a list

        EXAMPLES::


            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 2))
            sage: D._syndrome(r)
            [1, 10, 1, 10, 1]
        """
        C = self.code()
        F = C.base_ring()
        S = []
        col_mults = C.parity_column_multipliers()
        alphas = C.evaluation_points()

        for l in range(C.minimum_distance() - 1):
            Sl = F.zero()
            for j in range(C.length()):
                Sl += r[j] * col_mults[j] * (alphas[j] ** l)
            S.append(Sl)

        return S

    def _forney_formula(self, error_evaluator, error_locator):
        r"""
        Returns the error vector computed through Forney's formula.

        INPUT:

        - ``error_evaluator``, ``error_locator`` -- two polynomials

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: R.<x> = F[]
            sage: evaluator, locator = R(10), R([10, 10])
            sage: D._forney_formula(evaluator, locator)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
        """
        C = self.code()
        alphas = C.evaluation_points()
        col_mults = C.parity_column_multipliers()
        ELPp = error_locator.derivative()
        F = C.base_ring()
        zero, one = F.zero(), F.one()
        e = []

        for i in range(C.length()):
            alpha_inv = one/alphas[i]
            if error_locator(alpha_inv) == zero:
                e.append(-alphas[i]/col_mults[i] * error_evaluator(alpha_inv)/ELPp(alpha_inv))
            else:
                e.append(zero)

        return vector(F, e)

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``r`` and returns a codeword.

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be returned as is.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a vector of ``self.code()``

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: c == D.decode_to_code(y)
            True

        TESTS:

        If one tries to decode a word with too many errors, it returns
        an exception::

            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius()+1)
            sage: y = Chan(c)
            sage: D.decode_to_message(y)
            Traceback (most recent call last):
            ...
            DecodingError: Decoding failed because the number of errors exceeded the decoding radius

        If one tries to decode something which is not in the ambient space of the code,
        an exception is raised::

            sage: D.decode_to_code(42)
            Traceback (most recent call last):
            ...
            ValueError: The word to decode has to be in the ambient space of the code
        """
        C = self.code()
        if r not in C.ambient_space():
            raise ValueError("The word to decode has to be in the ambient space of the code")
        F = C.base_field()
        PolRing = C.base_field()['x']
        x = PolRing.gen()

        if C.length() == C.dimension() or r in C:
            return r

        S = PolRing(self._syndrome(r))
        a = x ** (C.minimum_distance() - 1)

        (EEP, ELP) = self._partial_xgcd(a, S, PolRing)

        e = self._forney_formula(EEP, ELP)
        dec = r - e
        if dec not in C:
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        return dec

    def decode_to_message(self, r):
        r"""
        Decodes``r`` to an element in message space of ``self``

        .. NOTE::

            If the code associated to ``self`` has the same length as its
            dimension, ``r`` will be unencoded as is. In that case,
            if ``r`` is not a codeword, the output is unspecified.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: D.connected_encoder().unencode(c) == D.decode_to_message(y)
            True
        """
        C = self.code()
        if C.length() == C.dimension():
            return self.connected_encoder().unencode_nocheck(r)
        return super(GRSKeyEquationSyndromeDecoder, self).decode_to_message(r)

    def decoding_radius(self):
        r"""
        Return maximal number of errors that ``self`` can decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D.decoding_radius()
            14
        """
        return (self.code().minimum_distance()-1)//2


####################### registration ###############################

GeneralizedReedSolomonCode._registered_encoders["EvaluationVector"] = GRSEvaluationVectorEncoder
GeneralizedReedSolomonCode._registered_encoders["EvaluationPolynomial"] = GRSEvaluationPolynomialEncoder

GeneralizedReedSolomonCode._registered_decoders["BerlekampWelch"] = GRSBerlekampWelchDecoder
GRSBerlekampWelchDecoder._decoder_type = {"hard-decision", "unique", "always-succeed"}
GeneralizedReedSolomonCode._registered_decoders["Gao"] = GRSGaoDecoder
GRSGaoDecoder._decoder_type = {"hard-decision", "unique", "always-succeed"}
GeneralizedReedSolomonCode._registered_decoders["ErrorErasure"] = GRSErrorErasureDecoder
GRSErrorErasureDecoder._decoder_type = {"error-erasure", "unique", "always-succeed"}
GeneralizedReedSolomonCode._registered_decoders["KeyEquationSyndrome"] = GRSKeyEquationSyndromeDecoder
GRSKeyEquationSyndromeDecoder._decoder_type = {"hard-decision", "unique", "always-succeed"}
