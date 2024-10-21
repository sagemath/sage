# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Gabidulin Code

This module provides the :class:`~sage.coding.gabidulin.GabidulinCode`, which constructs
Gabidulin Codes that are the rank metric equivalent of Reed Solomon codes and are
defined as the evaluation codes of degree-restricted skew polynomials.

This module also provides :class:`~sage.coding.gabidulin.GabidulinPolynomialEvaluationEncoder`,
an encoder with a skew polynomial message space and :class:`~sage.coding.gabidulin.GabidulinVectorEvaluationEncoder`,
an encoder based on the generator matrix. It also provides a decoder
:class:`~sage.coding.gabidulin.GabidulinGaoDecoder` which corrects errors using
the Gao algorithm in the rank metric.

AUTHOR:

- Arpit Merchant (2016-08-16)
- Marketa Slukova (2019-08-19): initial version
"""
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder, DecodingError
from sage.coding.linear_rank_metric import AbstractLinearRankMetricCode
from sage.categories.fields import Fields


class GabidulinCode(AbstractLinearRankMetricCode):
    r"""
    A Gabidulin Code.

    DEFINITION:

    A linear Gabidulin code Gab[n, k] over `F_{q^m}` of length `n` (at most
    `m`) and dimension `k` (at most `n`) is the set of all codewords, that
    are the evaluation of a `q`-degree restricted skew polynomial `f(x)`
    belonging to the skew polynomial constructed over the base ring `F_{q^m}`
    and the twisting homomorphism `\sigma`.

    .. math::

        \{ \text{Gab[n, k]} = \big\{ (f(g_0) f(g_1) ... f(g_{n-1})) = f(\textbf{g}) : \text{deg}_{q}f(x) < k \big\} \}

    where the fixed evaluation points `g_0, g_1,..., g_{n-1}` are linearly
    independent over `F_{q^m}`.

    EXAMPLES:

    A Gabidulin Code can be constructed in the following way::

        sage: Fqm = GF(16)
        sage: Fq = GF(4)
        sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
        sage: C
        [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, dimension, sub_field=None,
            twisting_homomorphism=None, evaluation_points=None):
        r"""
        Representation of a Gabidulin Code.

        INPUT:

        - ``base_field`` -- finite field of order `q^m` where `q` is a prime power
          and `m` is an integer

        - ``length`` -- length of the resulting code

        - ``dimension`` -- dimension of the resulting code

        - ``sub_field`` -- (default: ``None``) finite field of order `q`
          which is a subfield of the ``base_field``. If not given, it is the
          prime subfield of the ``base_field``.

        - ``twisting_homomorphism`` -- (default: ``None``) homomorphism of the
          underlying skew polynomial ring. If not given, it is the Frobenius
          endomorphism on ``base_field``, which sends an element `x` to `x^{q}`.

        - ``evaluation_points`` -- (default: ``None``) list of elements
          `g_0, g_1,...,g_{n-1}` of the ``base_field`` that are linearly
          independent over the ``sub_field``. These elements form the first row
          of the generator matrix. If not specified, these are the `nth` powers
          of the generator of the ``base_field``.

        Both parameters ``sub_field`` and ``twisting_homomorphism`` are optional.
        Since they are closely related, here is a complete list of behaviours:

        - both ``sub_field`` and ``twisting_homomorphism`` given -- in this case
          we only check that given that ``twisting_homomorphism`` has a fixed
          field method, it returns ``sub_field``

        - only ``twisting_homomorphism`` given -- we set ``sub_field`` to be the
          fixed field of the ``twisting_homomorphism``. If such method does not
          exist, an error is raised.

        - only ``sub_field`` given -- we set ``twisting_homomorphism`` to be the
          Frobenius of the field extension

        - neither ``sub_field`` or ``twisting_homomorphism`` given -- we take
          ``sub_field`` to be the prime field of ``base_field`` and the
          ``twisting_homomorphism`` to be the Frobenius wrt. the prime field

        TESTS:

        If ``length`` is bigger than the degree of the extension, an error is
        raised::

            sage: C = codes.GabidulinCode(GF(64), 4, 3, GF(4))
            Traceback (most recent call last):
            ...
            ValueError: 'length' can be at most the degree of the extension, 3

        If the number of evaluation points is not equal to the length
        of the code, an error is raised::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5)
            sage: aa = Fqm.gen()
            sage: evals = [ aa^i for i in range(21) ]
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq, None, evals)
            Traceback (most recent call last):
            ...
            ValueError: the number of evaluation points should be equal to the length of the code

        If evaluation points are not linearly independent over the ``base_field``,
        an error is raised::

            sage: evals = [ aa*i for i in range(2) ]
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq, None, evals)
            Traceback (most recent call last):
            ...
            ValueError: the evaluation points provided are not linearly independent

        If an evaluation point does not belong to the ``base_field``, an error
        is raised::

            sage: a = GF(3).gen()
            sage: evals = [ a*i for i in range(2) ]
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq, None, evals)
            Traceback (most recent call last):
            ...
            ValueError: evaluation point does not belong to the 'base field'

        Given that both ``sub_field`` and ``twisting_homomorphism`` are specified
        and ``twisting_homomorphism`` has a fixed field method. If the fixed
        field of ``twisting_homomorphism`` is not ``sub_field``, an error is
        raised::

            sage: Fqm = GF(64)
            sage: Fq = GF(8)
            sage: twist = GF(64).frobenius_endomorphism(n=2)
            sage: C = codes.GabidulinCode(Fqm, 3, 2, Fq, twist)
            Traceback (most recent call last):
            ...
            ValueError: the fixed field of the twisting homomorphism has to be the relative field of the extension

        If ``twisting_homomorphism`` is given, but ``sub_field`` is not. In case
        ``twisting_homomorphism`` does not have a fixed field method, and error
        is raised::

            sage: Fqm.<z6> = GF(64)
            sage: sigma = Hom(Fqm, Fqm)[1]; sigma
            Ring endomorphism of Finite Field in z6 of size 2^6
                Defn: z6 |--> z6^2
            sage: C = codes.GabidulinCode(Fqm, 3, 2, None, sigma)
            Traceback (most recent call last):
            ...
            ValueError: if 'sub_field' is not given, the twisting homomorphism has to have a 'fixed_field' method
        """
        twist_fix_field = None
        have_twist = (twisting_homomorphism is not None)
        have_subfield = (sub_field is not None)

        if have_twist and have_subfield:
            try:
                twist_fix_field = twisting_homomorphism.fixed_field()[0]
            except AttributeError:
                pass
            if twist_fix_field and twist_fix_field.order() != sub_field.order():
                raise ValueError("the fixed field of the twisting homomorphism has to be the relative field of the extension")

        if have_twist and not have_subfield:
            if not twist_fix_field:
                raise ValueError("if 'sub_field' is not given, the twisting homomorphism has to have a 'fixed_field' method")
            else:
                sub_field = twist_fix_field

        if (not have_twist) and have_subfield:
            twisting_homomorphism = base_field.frobenius_endomorphism(n=sub_field.degree())

        if (not have_twist) and not have_subfield:
            sub_field = base_field.base_ring()
            twisting_homomorphism = base_field.frobenius_endomorphism()

        self._twisting_homomorphism = twisting_homomorphism

        super().__init__(base_field, sub_field, length, "VectorEvaluation", "Gao")

        if length > self.extension_degree():
            raise ValueError("'length' can be at most the degree of the extension, {}".format(self.extension_degree()))
        if evaluation_points is None:
            evaluation_points = [base_field.gen()**i for i in range(base_field.degree())][:length]
        else:
            if not len(evaluation_points) == length:
                raise ValueError("the number of evaluation points should be equal to the length of the code")
            for i in range(length):
                if evaluation_points[i] not in base_field:
                    raise ValueError("evaluation point does not belong to the 'base field'")
            basis = self.matrix_form_of_vector(vector(evaluation_points))
            if basis.rank() != length:
                raise ValueError("the evaluation points provided are not linearly independent")
        self._evaluation_points = evaluation_points
        self._dimension = dimension

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq); C
            [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)
        """
        R = self.base_field()
        S = self.sub_field()
        if R and S in Fields():
            return "[%s, %s, %s] linear Gabidulin code over GF(%s)/GF(%s)" % (self.length(), self.dimension(), self.minimum_distance(), R.cardinality(), S.cardinality())
        else:
            return "[%s, %s, %s] linear Gabidulin code over %s/%s" % (self.length(), self.dimension(), self.minimum_distance(), R, S)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq);
            sage: latex(C)
            [2, 2, 1] \textnormal{ linear Gabidulin code over } \Bold{F}_{2^{4}}/\Bold{F}_{2^{2}}
        """
        txt = "[%s, %s, %s] \\textnormal{ linear Gabidulin code over } %s/%s"
        return txt % (self.length(), self.dimension(), self.minimum_distance(),
                      self.base_field()._latex_(), self.sub_field()._latex_())

    def __eq__(self, other):
        """
        Test equality between Gabidulin Code objects.

        INPUT:

        - ``other`` -- another Gabidulin Code object

        OUTPUT: boolean

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C1 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: C2 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: C1.__eq__(C2)
            True

            sage: Fqmm = GF(64)
            sage: C3 = codes.GabidulinCode(Fqmm, 2, 2, Fq)
            sage: C3.__eq__(C2)
            False
        """
        return isinstance(other, GabidulinCode) \
            and self.base_field() == other.base_field() \
            and self.sub_field() == other.sub_field() \
            and self.length() == other.length() \
            and self.dimension() == other.dimension() \
            and self.evaluation_points() == other.evaluation_points()

    def twisting_homomorphism(self):
        r"""
        Return the twisting homomorphism of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 5, 3, Fq)
            sage: C.twisting_homomorphism()
            Frobenius endomorphism z20 |--> z20^(5^4) on Finite Field in z20 of size 5^20
        """
        return self._twisting_homomorphism

    def minimum_distance(self):
        r"""
        Return the minimum distance of ``self``.

        Since Gabidulin Codes are Maximum-Distance-Separable (MDS), this returns
         ``self.length() - self.dimension() + 1``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5)
            sage: C = codes.GabidulinCode(Fqm, 20, 15, Fq)
            sage: C.minimum_distance()
            6
        """
        return self.length() - self.dimension() + 1

    def parity_evaluation_points(self):
        r"""
        Return the parity evaluation points of ``self``.

        These form the first row of the parity check matrix of ``self``.

        EXAMPLES::

            sage: C = codes.GabidulinCode(GF(2^10), 5, 2)
            sage: list(C.parity_check_matrix().row(0)) == C.parity_evaluation_points() #indirect_doctest
            True
        """
        eval_pts = self.evaluation_points()
        n = self.length()
        k = self.dimension()
        sigma = self.twisting_homomorphism()

        coefficient_matrix = matrix(self.base_field(), n - 1, n,
                                    lambda i, j: (sigma**(-n + k + 1 + i))(eval_pts[j]))
        solution_space = coefficient_matrix.right_kernel()
        return list(solution_space.basis()[0])

    def dual_code(self):
        r"""
        Return the dual code `C^{\perp}` of ``self``, the code `C`,

        .. MATH::

            C^{\perp} = \{ v \in V\ |\ v\cdot c = 0,\ \forall c \in C \}.

        EXAMPLES::

            sage: C = codes.GabidulinCode(GF(2^10), 5, 2)
            sage: C1 = C.dual_code(); C1
            [5, 3, 3] linear Gabidulin code over GF(1024)/GF(2)
            sage: C == C1.dual_code()
            True
        """
        return GabidulinCode(self.base_field(), self.length(),
                             self.length() - self.dimension(),
                             self.sub_field(),
                             self.twisting_homomorphism(),
                             self.parity_evaluation_points())

    def parity_check_matrix(self):
        r"""
        Return the parity check matrix of ``self``.

        This is the generator matrix of the dual code of ``self``.

        EXAMPLES::

            sage: C = codes.GabidulinCode(GF(2^3), 3, 2)
            sage: C.parity_check_matrix()
            [        1        z3 z3^2 + z3]
            sage: C.parity_check_matrix() == C.dual_code().generator_matrix()
            True
        """
        return self.dual_code().generator_matrix()

    def evaluation_points(self):
        """
        Return the evaluation points of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: C.evaluation_points()
            [1, z20, z20^2, z20^3]
        """
        return self._evaluation_points

# ---------------------- encoders ------------------------------


class GabidulinVectorEvaluationEncoder(Encoder):

    def __init__(self, code):
        """
        This method constructs the vector evaluation encoder for
        Gabidulin Codes.

        INPUT:

        - ``code`` -- the associated code of this encoder

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E = codes.encoders.GabidulinVectorEvaluationEncoder(C)
            sage: E
            Vector evaluation style encoder for [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)

        Alternatively, we can construct the encoder from ``C`` directly::

            sage: E = C.encoder("VectorEvaluation")
            sage: E
            Vector evaluation style encoder for [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)

        TESTS:

        If the code is not a Gabidulin code, an error is raised::

            sage: C = codes.HammingCode(GF(4), 2)
            sage: E = codes.encoders.GabidulinVectorEvaluationEncoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a Gabidulin code
        """
        if not isinstance(code, GabidulinCode):
            raise ValueError("code has to be a Gabidulin code")
        super().__init__(code)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: E = codes.encoders.GabidulinVectorEvaluationEncoder(C); E
            Vector evaluation style encoder for [4, 4, 1] linear Gabidulin code over GF(95367431640625)/GF(625)
        """
        return "Vector evaluation style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: E = codes.encoders.GabidulinVectorEvaluationEncoder(C)
            sage: latex(E)
            \textnormal{Vector evaluation style encoder for } [4, 4, 1] \textnormal{ linear Gabidulin code over } \Bold{F}_{5^{20}}/\Bold{F}_{5^{4}}
        """
        return "\\textnormal{Vector evaluation style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        """
        Test equality between Gabidulin Generator Matrix Encoder objects.

        INPUT:

        - ``other`` -- another Gabidulin Generator Matrix Encoder

        OUTPUT: boolean

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C1 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E1 = codes.encoders.GabidulinVectorEvaluationEncoder(C1)
            sage: C2 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E2 = codes.encoders.GabidulinVectorEvaluationEncoder(C2)
            sage: E1.__eq__(E2)
            True

            sage: Fqmm = GF(64)
            sage: C3 = codes.GabidulinCode(Fqmm, 2, 2, Fq)
            sage: E3 = codes.encoders.GabidulinVectorEvaluationEncoder(C3)
            sage: E3.__eq__(E2)
            False
        """
        return isinstance(other, GabidulinVectorEvaluationEncoder) \
            and self.code() == other.code()

    def generator_matrix(self):
        """
        Return the generator matrix of ``self``.

        EXAMPLES::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 3, 3, Fq)
            sage: (list(C.generator_matrix().row(1))
            ....:   == [C.evaluation_points()[i]**(2**3) for i in range(3)])
            True
        """
        from functools import reduce
        C = self.code()
        eval_pts = C.evaluation_points()
        sigma = C.twisting_homomorphism()

        def create_matrix_elements(A, k, f):
            return reduce(lambda L, x: [x] +
                          [list(map(f, l)) for l in L], [A] * k, [])
        return matrix(C.base_field(), C.dimension(), C.length(),
                      create_matrix_elements(eval_pts, C.dimension(), sigma))


class GabidulinPolynomialEvaluationEncoder(Encoder):
    r"""
    Encoder for Gabidulin codes which uses evaluation of skew polynomials to
    obtain codewords.

    Let `C` be a Gabidulin code of length `n` and dimension `k` over some
    finite field `F = GF(q^m)`. We denote by `\alpha_i` its evaluations
    points, where `1 \leq i \leq n`. Let `p`, a skew polynomial of degree at
    most `k-1` in `F[x]`, be the message.

    The encoding of `m` will be the following codeword:

    .. MATH::

        (p(\alpha_1), \dots, p(\alpha_n)).

    TESTS:

    This module uses the following experimental feature.
    This test block is here only to trigger the experimental warning so it does not
    interferes with doctests::

        sage: Fqm = GF(2^9)
        sage: Fq = GF(2^3)
        sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
        sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
        sage: z9 = Fqm.gen()
        sage: p = (z9^6 + z9^2 + z9 + 1)*x + z9^7 + z9^5 + z9^4 + z9^2
        sage: vector(p.multi_point_evaluation(C.evaluation_points()))
        doctest:...: FutureWarning: This class/method/function is marked as experimental.
        It, its functionality or its interface might change without a formal deprecation.
        See https://github.com/sagemath/sage/issues/13215 for details.
        (z9^7 + z9^6 + z9^5 + z9^4 + z9 + 1, z9^6 + z9^5 + z9^3 + z9)

    EXAMPLES::

        sage: Fqm = GF(16)
        sage: Fq = GF(4)
        sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
        sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
        sage: E
        Polynomial evaluation style encoder for
         [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)

    Alternatively, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("PolynomialEvaluation")
        sage: E
        Polynomial evaluation style encoder for
         [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)
    """

    def __init__(self, code):
        r"""
        INPUT:

        - ``code`` -- the associated code of this encoder

        TESTS:

        If the code is not a Gabidulin code, an error is raised::

            sage: C = codes.HammingCode(GF(4), 2)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a Gabidulin code
        """
        if not isinstance(code, GabidulinCode):
            raise ValueError("code has to be a Gabidulin code")
        super().__init__(code)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C); E
            Polynomial evaluation style encoder for [4, 4, 1] linear Gabidulin code over GF(95367431640625)/GF(625)
        """
        return "Polynomial evaluation style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: latex(E)
            \textnormal{Polynomial evaluation style encoder for } [4, 4, 1] \textnormal{ linear Gabidulin code over } \Bold{F}_{5^{20}}/\Bold{F}_{5^{4}}
        """
        return "\\textnormal{Polynomial evaluation style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        """
        Test equality between Gabidulin Polynomial Evaluation Encoder objects.

        INPUT:

        - ``other`` -- another Gabidulin Polynomial Evaluation Encoder

        OUTPUT: boolean

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C1 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E1 = codes.encoders.GabidulinPolynomialEvaluationEncoder(C1)
            sage: C2 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E2 = codes.encoders.GabidulinPolynomialEvaluationEncoder(C2)
            sage: E1.__eq__(E2)
            True

            sage: Fqmm = GF(64)
            sage: C3 = codes.GabidulinCode(Fqmm, 2, 2, Fq)
            sage: E3 = codes.encoders.GabidulinPolynomialEvaluationEncoder(C3)
            sage: E3.__eq__(E2)
            False
        """
        return isinstance(other, GabidulinPolynomialEvaluationEncoder) \
            and self.code() == other.code()

    def message_space(self):
        r"""
        Return the message space of the associated code of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: E.message_space()
            Ore Polynomial Ring in x over Finite Field in z20 of size 5^20
             twisted by z20 |--> z20^(5^4)
        """
        C = self.code()
        return C.base_field()['x', C.twisting_homomorphism()]

    def encode(self, p, form='vector'):
        """
        Transform the polynomial ``p`` into a codeword of :meth:`code`.

        The output codeword can be represented as a vector or a matrix,
        depending on the ``form`` input.

        INPUT:

        - ``p`` -- a skew polynomial from the message space of ``self`` of degree
          less than ``self.code().dimension()``

        - ``form`` -- type parameter taking strings "vector" or "matrix"
          as values and converting the output codeword into the respective form
          (default: ``'vector'``)

        OUTPUT: a codeword corresponding to `p` in vector or matrix form

        EXAMPLES::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z9 = Fqm.gen()
            sage: p = (z9^6 + z9^2 + z9 + 1)*x + z9^7 + z9^5 + z9^4 + z9^2
            sage: codeword_vector = E.encode(p, "vector"); codeword_vector
            (z9^7 + z9^6 + z9^5 + z9^4 + z9 + 1, z9^6 + z9^5 + z9^3 + z9)
            sage: codeword_matrix = E.encode(p, "matrix"); codeword_matrix
            [           z3     z3^2 + z3]
            [           z3             1]
            [         z3^2 z3^2 + z3 + 1]

        TESTS:

        If the skew polynomial, `p`, has degree greater than or equal to the
        dimension of the code, an error is raised::

            sage: t = z9^4*x^2 + z9
            sage: codeword_vector = E.encode(t, "vector"); codeword_vector
            Traceback (most recent call last):
            ...
            ValueError: the skew polynomial to encode must have degree at most 1

        The skew polynomial, `p`, must belong to the message space of the code.
        Otherwise, an error is raised::

            sage: Fqmm = GF(2^12)
            sage: S.<x> = Fqmm['x', Fqmm.frobenius_endomorphism(n=3)]
            sage: q = S.random_element(degree=2)
            sage: codeword_vector = E.encode(q, "vector"); codeword_vector
            Traceback (most recent call last):
            ...
            ValueError: the message to encode must be in Ore Polynomial Ring in x over Finite Field in z9 of size 2^9 twisted by z9 |--> z9^(2^3)
        """
        C = self.code()
        M = self.message_space()
        if p not in M:
            raise ValueError("the message to encode must be in %s" % M)
        if p.degree() >= C.dimension():
            raise ValueError("the skew polynomial to encode must have degree at most %s" % (C.dimension() - 1))
        eval_pts = C.evaluation_points()
        codeword = p.multi_point_evaluation(eval_pts)
        if form == "vector":
            return vector(codeword)
        elif form == "matrix":
            return C.matrix_form_of_vector(vector(codeword))
        else:
            return ValueError("the argument 'form' takes only either 'vector' or 'matrix' as valid input")

    def unencode_nocheck(self, c):
        """
        Return the message corresponding to the codeword ``c``.

        Use this method with caution: it does not check if ``c``
        belongs to the code, and if this is not the case, the output is
        unspecified.

        INPUT:

        - ``c`` -- a codeword of :meth:`code`

        OUTPUT:

        - a skew polynomial of degree less than ``self.code().dimension()``

        EXAMPLES::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z9 = Fqm.gen()
            sage: p = (z9^6 + z9^4)*x + z9^2 + z9
            sage: codeword_vector = E.encode(p, "vector")
            sage: E.unencode_nocheck(codeword_vector)
            (z9^6 + z9^4)*x + z9^2 + z9
        """
        C = self.code()
        eval_pts = C.evaluation_points()
        values = [c[i] for i in range(len(c))]
        points = [(eval_pts[i], values[i]) for i in range(len(eval_pts))]
        p = self.message_space().lagrange_polynomial(points)
        return p


# ---------------------- decoders ------------------------------


class GabidulinGaoDecoder(Decoder):

    def __init__(self, code):
        r"""
        Gao style decoder for Gabidulin Codes.

        INPUT:

        - ``code`` -- the associated code of this decoder

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: D
            Gao decoder for [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)

        Alternatively, we can construct the encoder from ``C`` directly::

            sage: D = C.decoder("Gao")
            sage: D
            Gao decoder for [2, 2, 1] linear Gabidulin code over GF(16)/GF(4)

        TESTS:

        If the code is not a Gabidulin code, an error is raised::

            sage: C = codes.HammingCode(GF(4), 2)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a Gabidulin code
        """
        if not isinstance(code, GabidulinCode):
            raise ValueError("code has to be a Gabidulin code")
        super().__init__(code, code.ambient_space(), "PolynomialEvaluation")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C); D
            Gao decoder for [4, 4, 1] linear Gabidulin code over GF(95367431640625)/GF(625)
        """
        return "Gao decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5^4)
            sage: C = codes.GabidulinCode(Fqm, 4, 4, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: latex(D)
            \textnormal{Gao decoder for } [4, 4, 1] \textnormal{ linear Gabidulin code over } \Bold{F}_{5^{20}}/\Bold{F}_{5^{4}}
        """
        return "\\textnormal{Gao decoder for } %s" % self.code()._latex_()

    def __eq__(self, other) -> bool:
        """
        Test equality between Gabidulin Gao Decoder objects.

        INPUT:

        - ``other`` -- another Gabidulin Gao Decoder

        OUTPUT: boolean

        EXAMPLES::

            sage: Fqm = GF(16)
            sage: Fq = GF(4)
            sage: C1 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: D1 = codes.decoders.GabidulinGaoDecoder(C1)
            sage: C2 = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: D2 = codes.decoders.GabidulinGaoDecoder(C2)
            sage: D1.__eq__(D2)
            True

            sage: Fqmm = GF(64)
            sage: C3 = codes.GabidulinCode(Fqmm, 2, 2, Fq)
            sage: D3 = codes.decoders.GabidulinGaoDecoder(C3)
            sage: D3.__eq__(D2)
            False
        """
        return isinstance(other, GabidulinGaoDecoder) \
            and self.code() == other.code()

    def _partial_xgcd(self, a, b, d_stop):
        """
        Compute the partial gcd of `a` and `b` using the right linearized
        extended Euclidean algorithm up to the `d_stop` iterations. This
        is a private method for internal use only.

        INPUT:

        - ``a`` -- a skew polynomial

        - ``b`` -- another skew polynomial

        - ``d_stop`` -- the number of iterations for which the algorithm
          is to be run

        OUTPUT: ``r_c`` -- right linearized remainder of `a` and `b`

        - ``u_c`` -- right linearized quotient of `a` and `b`

        EXAMPLES::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z9 = Fqm.gen()
            sage: p = (z9^6 + z9^4)*x + z9^2 + z9
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_message(codeword_vector) #indirect_doctest
            sage: r
            (z9^6 + z9^4)*x + z9^2 + z9
        """
        S = self.message_space()
        if (a not in S) or (b not in S):
            raise ValueError("both the input polynomials must belong to %s" % S)
        if a.degree() < b.degree():
            raise ValueError("degree of first polynomial must be greater than or equal to degree of second polynomial")
        r_p = a
        r_c = b
        u_p = S.zero()
        u_c = S.one()
        v_p = u_c
        v_c = u_p

        while r_c.degree() >= d_stop:
            (q, r_c), r_p = r_p.right_quo_rem(r_c), r_c
            u_c, u_p = u_p - q * u_c, u_c
            v_c, v_p = v_p - q * v_c, v_c
        return r_c, u_c

    def _decode_to_code_and_message(self, r):
        """
        Return the decoded codeword and message (skew polynomial)
        corresponding to the received codeword `r`. This is a
        private method for internal use only.

        INPUT:

        - ``r`` -- received codeword

        OUTPUT:

        - the decoded codeword and decoded message corresponding to
          the received codeword `r`

        EXAMPLES::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z9 = Fqm.gen()
            sage: p = (z9^6 + z9^4)*x + z9^2 + z9
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_message(codeword_vector) #indirect doctest
            sage: r
            (z9^6 + z9^4)*x + z9^2 + z9
        """
        C = self.code()
        length = len(r)
        eval_pts = C.evaluation_points()
        S = self.message_space()

        if length == C.dimension() or r in C:
            return r, self.connected_encoder().unencode_nocheck(r)

        points = [(eval_pts[i], r[i]) for i in range(len(eval_pts))]
        # R = S.lagrange_polynomial(eval_pts, list(r))
        R = S.lagrange_polynomial(points)
        r_out, u_out = self._partial_xgcd(S.minimal_vanishing_polynomial(eval_pts),
                R, (C.length() + C.dimension()) // 2)
        quo, rem = r_out.left_quo_rem(u_out)
        if not rem.is_zero():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        if quo not in S:
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        c = self.connected_encoder().encode(quo)
        if C.rank_weight_of_vector(c - r) > self.decoding_radius():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        return c, quo

    def decode_to_code(self, r):
        """
        Return the decoded codeword corresponding to the
        received word `r`.

        INPUT:

        - ``r`` -- received codeword

        OUTPUT: the decoded codeword corresponding to the received codeword

        EXAMPLES::

            sage: Fqm = GF(3^20)
            sage: Fq = GF(3)
            sage: C = codes.GabidulinCode(Fqm, 5, 3, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z20 = Fqm.gen()
            sage: p = x
            sage: codeword_vector = E.encode(p, "vector")
            sage: codeword_vector
            (1, z20^3, z20^6, z20^9, z20^12)
            sage: l = list(codeword_vector)
            sage: l[0] = l[1] #make an error
            sage: D.decode_to_code(vector(l))
            (1, z20^3, z20^6, z20^9, z20^12)
        """
        return self._decode_to_code_and_message(r)[0]

    def decode_to_message(self, r):
        """
        Return the skew polynomial (message) corresponding to the
        received word `r`.

        INPUT:

        - ``r`` -- received codeword

        OUTPUT: the message corresponding to the received codeword

        EXAMPLES::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z9 = Fqm.gen()
            sage: p = (z9^6 + z9^4)*x + z9^2 + z9
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_message(codeword_vector)
            sage: r
            (z9^6 + z9^4)*x + z9^2 + z9
        """
        return self._decode_to_code_and_message(r)[1]

    def decoding_radius(self):
        """
        Return the decoding radius of the Gabidulin Gao Decoder.

        EXAMPLES::

            sage: Fqm = GF(5^20)
            sage: Fq = GF(5)
            sage: C = codes.GabidulinCode(Fqm, 20, 4, Fq)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: D.decoding_radius()
            8
        """
        return (self.code().minimum_distance() - 1) // 2

# ----------------------------- registration --------------------------------


GabidulinCode._registered_encoders["PolynomialEvaluation"] = GabidulinPolynomialEvaluationEncoder
GabidulinCode._registered_encoders["VectorEvaluation"] = GabidulinVectorEvaluationEncoder

GabidulinCode._registered_decoders["Gao"] = GabidulinGaoDecoder
