# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Generic structures for linear codes over the rank metric

Rank Metric
===========

In coding theory, the most common metric is the Hamming metric, where distance
between two codewords is  given by the number of positions in which they differ.
An alternative to this is the rank metric. Take two fields, `\GF{q}` and `\GF{q^m}`,
and define a code `C` to be a set of vectors of length `n` with entries from
`\GF{q^m}`. Let `c` be a codeword. We can represent it as an `m \times n` matrix
`M` over `\GF{q}`.

A detailed description on the relationship between the two representations can
be found in :meth:`sage.coding.linear_rank_metric.to_matrix_representation`
and :meth:`sage.coding.linear_rank_metric.from_matrix_representation`.

We can define a metric using the rank of the matrix representation of the
codewords. A distance between two codewords `a, b` is the rank of the matrix
representation of `a - b`. A weight of a codeword `c` is the rank of the matrix
representation of `c`.

This module allows representing rank metric codes which are linear over the
big field `\GF{q^m}`, i.e. the usual linearity condition when the codewords are
considered in vector form. One can also consider rank metric codes which are
only linear over `\GF{q}`, but these are not currently supported in SageMath.

Note that linear rank metric codes per the definition of this file are
mathematically just linear block codes, and so could be considered as a
:class:`sage.coding.linear_code.LinearCode`. However, since most of the
functionality of that class is specific to the Hamming metric, the two notions
are implemented as entirely different in SageMath. If you wish to investigate
Hamming-metric properties of a linear rank metric code ``C``, you can easily
convert it by calling ``C_hamm = LinearCode(C)``.

Linear Rank Metric Code and Gabidulin Codes
===========================================

The class :class:`sage.coding.linear_rank_metric.LinearRankMetricCode` is the
analog of :class:`sage.coding.linear_code.LinearCode`, i.e. it is a generator
matrix-based representation of a linear rank metric code without specific
knowledge on the structure of the code.

Gabidulin codes are the main family of structured linear rank metric codes.
These codes are the rank-metric analog of Reed-Solomon codes.

``AbstractLinearRankMetricCode``
--------------------------------

This is a base class designed to contain methods, features and parameters
shared by every linear rank metric code. For instance, generic algorithms for
computing the minimum distance, etc. Many of these algorithms are slow,
e.g. exponential in the code length. It also contains methods for swapping
between vector and matrix representation of elements.

``AbstractLinearCodeNoMetric`` is an abstract class for linear rank metric codes,
so any linear rank metric code  class should inherit from this class.
Also ``AbstractLinearCodeNoMetric`` should never itself be instantiated.

See :class:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode`
for details and examples.

``LinearRankMetricCode``
------------------------

This class is used to represent arbitrary and unstructured linear rank metric
codes. It mostly relies directly on generic methods provided by
``AbstractLinearRankMetricCode``, which means that basic operations on the code
(e.g. computation of the minimum distance) will use slow algorithms.

A ``LinearRankMetricCode`` is instantiated by providing a generator::

    sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
    sage: C = codes.LinearRankMetricCode(G, GF(4))
    sage: C
    [3, 2] linear rank metric code over GF(64)/GF(4)
    sage: C.generator_matrix()
    [1 1 0]
    [0 0 1]
    sage: c = vector(GF(64), (1, 1, 1))
    sage: c in C
    True

Further references
------------------

Read more about
:wikipedia:`rank metric and Gabidulin codes <Rank_error-correcting_code>`

AUTHORS:

- Marketa Slukova (2019-08-16): initial version
- Camille Garnier and Rubén Muñoz--Bertrand (2024-02-13): added rank_support_of_vector, and corrected the documentation

TESTS::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C = LinearCode(G)
    sage: C == loads(dumps(C))
    True
"""

# ****************************************************************************
#       Copyright (C) 2019 MARKETA SLUKOVA <em.slukova@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.fields import Fields
from sage.matrix.constructor import matrix
from sage.structure.element import Matrix, Vector
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity

from .linear_code_no_metric import AbstractLinearCodeNoMetric
from .linear_code import LinearCodeGeneratorMatrixEncoder
from .decoder import Decoder


def to_matrix_representation(v, sub_field=None, basis=None):
    r"""
    Return a matrix representation of ``v`` over ``sub_field`` in terms of
    ``basis``.

    Let `(b_1, b_2, \ldots, b_m)`, `b_i \in GF(q^m)`, be a basis of `\GF{q^m}` as
    a vector space over `\GF{q}`. Take an element `x \in \GF{q^m}`. We can write
    `x` as `x = u_1 b_1 + u_2 b_2 + \ldots u_m b_m`, where `u_i \in GF(q)`. This
    way we can represent an element from `\GF{q^m}` as a vector of length `m`
    over `\GF{q}`.

    Given a vector ``v`` of length `n` over some field `\GF{q^m}`, we can
    represent each entry as a vector of length `m`, yielding an `m \times n`
    matrix over ``sub_field``. In case ``sub_field`` is not given, we take the
    prime subfield `\GF{p}` of `\GF{q^m}`.

    INPUT:

    - ``v`` -- a vector over some field `\GF{q^m}`

    - ``sub_field`` -- (default: ``None``) a sub field of `\GF{q^m}`; if not
      specified, it is the prime subfield `\GF{p}` of `\GF{q^m}`

    - ``basis`` -- (default: ``None``) a basis of `\GF{q^m}` as a vector space over
      ``sub_field``. If not specified, the default basis is
      `1,\beta,\ldots,\beta^{m-1}` where `\beta` is the generator of the
      multiplicative group of `\GF{q^m}` given by Sage.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import to_matrix_representation
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: to_matrix_representation(a, GF(4))
        [1 1 1]
        [1 1 0]
        [0 0 0]

        sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
        sage: to_matrix_representation(m)
        Traceback (most recent call last):
        ...
        TypeError: Input must be a vector
    """
    if not isinstance(v, Vector):
        raise TypeError("Input must be a vector")
    base_field = v.base_ring()
    if not sub_field:
        sub_field = base_field.prime_subfield()
    n = v.length()
    m = base_field.degree()//sub_field.degree()
    extension, to_big_field, from_big_field = base_field.vector_space(sub_field, basis, map=True)
    return matrix(sub_field, m, n, lambda i, j: from_big_field(v[j])[i])


def from_matrix_representation(w, base_field=None, basis=None):
    r"""
    Return a vector representation of a matrix ``w`` over ``base_field`` in terms
    of ``basis``.

    Given an `m \times n` matrix over `\GF{q}` and some ``basis`` of `\GF{q^m}`
    over `\GF{q}`, we can represent each of its columns as an element of `\GF{q^m}`,
    yielding a vector of length `n` over `\GF{q}`.

    In case ``base_field`` is not given, we take `\GF{q^m}`, the field extension of
    `\GF{q}` of degree `m`, the number of rows of ``w``.

    INPUT:

    - ``w`` -- a matrix over some field `\GF{q}`

    - ``base_field`` -- (default: ``None``) an extension field of `\GF{q}`. If not
      specified, it is the field `\GF{q^m}`, where `m` is the number of rows of
      ``w``.

    - ``basis`` -- (default: ``None``) a basis of `\GF{q^m}` as a vector space over
      `\GF{q}`. If not specified, the default basis is
      `1,\beta,\ldots,\beta^{m-1}` where `\beta` is the generator of the
      multiplicative group of `\GF{q^m}` given by Sage.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import from_matrix_representation
        sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
        sage: from_matrix_representation(m)
        (z6 + 1, z6 + 1, 1)

        sage: v = vector(GF(4), (1, 0, 0))
        sage: from_matrix_representation(v)
        Traceback (most recent call last):
        ...
        TypeError: Input must be a matrix
    """
    if not isinstance(w, Matrix):
        raise TypeError("Input must be a matrix")
    sub_field = w.base_ring()
    if not base_field:
        base_field = sub_field.extension(w.nrows())
    extension, to_big_field, from_big_field = base_field.vector_space(sub_field, basis, map=True)
    v = [to_big_field(w.column(i)) for i in range(w.ncols())]
    return vector(v)


def rank_weight(c, sub_field=None, basis=None):
    r"""
    Return the rank of ``c`` as a matrix over ``sub_field``.

    If ``c`` is a vector over some field `\GF{q^m}`, the function converts it
    into a matrix over ``sub_field```.

    INPUT:

    - ``c`` -- a vector over some field `\GF{q^m}`; or a matrix over `\GF{q}`

    - ``sub_field`` -- (default: ``None``) a sub field of `\GF{q^m}`; if not
      specified, it is the prime subfield `\GF{p}` of `\GF{q^m}`

    - ``basis`` -- (default: ``None``) a basis of `\GF{q^m}` as a vector space over
      ``sub_field``. If not specified, the default basis is
      `1,\beta,\ldots,\beta^{m-1}` where `\beta` is the generator of the
      multiplicative group of `\GF{q^m}` given by Sage.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import rank_weight
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: rank_weight(a, GF(4))
        2
    """
    if isinstance(c, Vector):
        c = to_matrix_representation(c, sub_field, basis)
    return c.rank()


def rank_distance(a, b, sub_field=None, basis=None):
    r"""
    Return the rank of ``a`` - ``b`` as a matrix over ``sub_field``.

    Take two vectors ``a``, ``b`` over some field `\GF{q^m}`. This function
    converts them to matrices over `\GF{q}` and calculates the rank of their
    difference.

    If ``sub_field`` is not specified, we take the prime subfield `\GF{q}` of
    `\GF{q^m}`.

    INPUT:

    - ``a`` -- a vector over some field `\GF{q^m}`

    - ``b`` -- a vector over some field `\GF{q^m}`

    - ``sub_field`` -- (default: ``None``) a sub field of `\GF{q^m}`; if not
      specified, it is the prime subfield `\GF{p}` of `\GF{q^m}`

    - ``basis`` -- (default: ``None``) a basis of `\GF{q^m}` as a vector space over
      ``sub_field``. If not specified, the default basis is
      `1,\beta,\ldots,\beta^{m-1}` where `\beta` is the generator of the
      multiplicative group of `\GF{q^m}` given by Sage.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import rank_distance
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: b = vector(GF(64), (1, 0, 0))
        sage: rank_distance(a, b, GF(4))
        2

        sage: c = vector(GF(4), (1, 0, 0))
        sage: rank_distance(a, c, GF(4))
        Traceback (most recent call last):
        ...
        ValueError: The base field of (z6 + 1, z6 + 1, 1) and (1, 0, 0) has to be the same

        sage: d = Matrix(GF(64), (1, 0, 0))
        sage: rank_distance(a, d, GF(64))
        Traceback (most recent call last):
        ...
        TypeError: Both inputs have to be vectors

        sage: e = vector(GF(64), (1, 0))
        sage: rank_distance(a, e, GF(64))
        Traceback (most recent call last):
        ...
        ValueError: The length of (z6 + 1, z6 + 1, 1) and (1, 0) has to be the same
    """
    if not (a.base_ring() == b.base_ring()):
        raise ValueError("The base field of {} and {} has to be the same".format(a, b))
    if not (isinstance(a, Vector) and isinstance(b, Vector)):
        raise TypeError("Both inputs have to be vectors")
    if not len(a) == len(b):
        raise ValueError("The length of {} and {} has to be the same".format(a, b))

    a = to_matrix_representation(a, sub_field, basis)
    b = to_matrix_representation(b, sub_field, basis)
    return (a - b).rank()


class AbstractLinearRankMetricCode(AbstractLinearCodeNoMetric):
    r"""
    Abstract class for linear rank metric codes.

    This class contains methods that can be used on families of linear rank
    metric codes. Every linear rank metric code class should inherit from this
    abstract class.

    This class is intended for codes which are linear over the ``base_field``.

    Codewords of rank metric codes have two representations. They can either be
    written as a vector of length `n` over `\GF{q^m}`, or an `m \times n` matrix
    over `\GF{q}`. This implementation principally uses the vector representation.
    However, one can always get the matrix representation using the
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.to_matrix`
    method. To go back to a vector, use the
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.from_matrix`
    method.

    Instructions on how to make a new family of rank metric codes is analogous
    to making a new family of linear codes over the Hamming metric, instructions
    for which are in :class:`sage.coding.linear_code.AbstractLinearCode`. For an
    example on, see
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.__init__`

    .. WARNING::

        A lot of methods of the abstract class rely on the knowledge of a generator matrix.
        It is thus strongly recommended to set an encoder with a generator matrix implemented
        as a default encoder.
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, sub_field, length, default_encoder_name,
            default_decoder_name, basis=None):
        r"""
        Initialize mandatory parameters that every linear rank metric code has.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear rank metric code.
        The class :class:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode`
        should never be directly instantiated.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``sub_field`` -- the sub field of ``self``

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer),
          must be > 0 and at most the degree of the field extension

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``

        - ``basis`` -- (default: ``None``) a basis of `\GF{q^m}` as a vector space over
          ``sub_field``. If not specified, the default basis is
          `1,\beta,\ldots,\beta^{m-1}` where `\beta` is the generator of the
          multiplicative group of `\GF{q^m}` given by Sage.

        EXAMPLES:

        The following example demonstrates how to use subclass
        `AbstractLinearRankMetricCode` for representing a new family of rank
        metric codes. The example is a rank repetition code::

             sage: from sage.coding.linear_rank_metric import AbstractLinearRankMetricCode
             sage: class RankRepetitionCode(AbstractLinearRankMetricCode):
             ....:   def __init__(self, base_field, sub_field, length):
             ....:       super().__init__(base_field, sub_field, length,
             ....:                        "GeneratorMatrix", "NearestNeighbor")
             ....:       beta = base_field.gen()
             ....:       self._generator_matrix = matrix(base_field,
             ....:                                       [[beta^i for i in range(length)]])
             ....:   def generator_matrix(self):
             ....:       return self._generator_matrix
             ....:   def _repr_(self):
             ....:       return "[%d, %d] rank-metric repetition code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

        We now instantiate a member of our newly made code family::

            sage: C = RankRepetitionCode(GF(8), GF(2), 3)

        We can check its existence and parameters::

            sage: C
            [3, 1] rank-metric repetition code over GF(8)

        We can encode a vector::

            sage: word = vector(C.base_field(), [1])
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: codeword = E(word)
            sage: codeword
            (1, z3, z3^2)

        We can get the matrix representation of the codeword::

            sage: C.matrix_form_of_vector(codeword)
            [1 0 0]
            [0 1 0]
            [0 0 1]

        We can decode the vector representation of the codeword::

            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decode_to_code(codeword)
            (1, z3, z3^2)
            sage: D.decode_to_message(codeword)
            (1)

        We can check that it is truly a part of the framework category::

            sage: C.parent()
            <class '__main__.RankRepetitionCode_with_category'>
            sage: C.category()
            Category of facade finite dimensional vector spaces with basis
             over Finite Field in z3 of size 2^3

        And any method that works on rank metric linear codes works for our new dummy code::

            sage: C.minimum_distance()
            3
            sage: C.metric()
            'rank'

        TESTS:

        If ``sub_field`` is not a field, an error is raised::

            sage: C = RankRepetitionCode(GF(8), ZZ, 3)
            Traceback (most recent call last):
            ...
            ValueError: 'sub_field' must be a field (and Integer Ring is not one)

        If ``sub_field`` is not a subfield of ``base_field``, an error is raised::

            sage: C = RankRepetitionCode(GF(8), GF(3), 2)
            Traceback (most recent call last):
            ...
            ValueError: 'sub_field' has to be a subfield of 'base_field'
        """
        self._registered_decoders["NearestNeighbor"] = LinearRankMetricCodeNearestNeighborDecoder

        if not sub_field.is_field():
            raise ValueError("'sub_field' must be a field (and {} is not one)".format(sub_field))
        if not (sub_field.degree().divides(base_field.degree()) and (sub_field.prime_subfield() == base_field.prime_subfield())):
            raise ValueError("'sub_field' has to be a subfield of 'base_field'")
        m = base_field.degree() // sub_field.degree()
        self._extension_degree = m
        self._sub_field = sub_field

        self._generic_constructor = LinearRankMetricCode
        super().__init__(base_field, length, default_encoder_name,
                         default_decoder_name, "rank")

    def sub_field(self):
        r"""
        Return the sub field of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.sub_field()
            Finite Field in z2 of size 2^2
        """
        return self._sub_field

    def extension_degree(self):
        r"""
        Return `m`, the degree of the field extension of ``self``.

        Let ``base_field`` be `\GF{q^m}` and ``sub_field`` be `\GF{q}`. Then this
        function returns `m`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.extension_degree()
            3
        """

        return self._extension_degree

    def field_extension(self):
        r"""
        Return the field extension of ``self``.

        Let ``base_field`` be some field `\GF{q^m}` and ``sub_field`` `\GF{q}`.
        This function returns the vector space of dimension `m` over `\GF{q}`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.field_extension()
            Vector space of dimension 3 over Finite Field in z2 of size 2^2
        """
        return self.base_field().vector_space(self.sub_field(), map=False)

    def rank_distance_between_vectors(self, left, right):
        r"""
        Return the rank of the matrix of ``left`` - ``right``.

        INPUT:

        - ``left`` -- a vector over the ``base_field`` of ``self``

        - ``right`` -- a vector over the ``base_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: b = vector(GF(64), (1, 0, 0))
            sage: C.rank_distance_between_vectors(a, b)
            2
        """
        return rank_distance(left, right, self.sub_field())

    def minimum_distance(self):
        r"""
        Return the minimum distance of ``self``.

        This algorithm simply iterates over all the elements of the code and
        returns the minimum weight.

        EXAMPLES::

            sage: F.<a> = GF(8)
            sage: G = Matrix(F, [[1,a,a^2,0]])
            sage: C = codes.LinearRankMetricCode(G, GF(2))
            sage: C.minimum_distance()
            3
        """
        d = Infinity
        for c in self:
            if c == self.zero():
                continue
            d = min(self.rank_weight_of_vector(c), d)
        return d

    def rank_weight_of_vector(self, word):
        r"""
        Return the weight of the word, i.e. its rank.

        INPUT:

        - ``word`` -- a vector over the ``base_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.rank_weight_of_vector(a)
            2
        """
        return rank_weight(word, self.sub_field())

    def rank_support_of_vector(self, word, sub_field=None, basis=None):
        r"""
        Return the rank support of ``word`` over ``sub_field``, i.e. the  vector space over
        ``sub_field`` generated by its coefficients.

        If ``word`` is a vector over some field `\GF{q^m}`, and ``sub_field`` is a subfield of
        `\GF{q^m}`, the function converts it into a matrix over ``sub_field``, with
        respect to the basis ``basis``.

        INPUT:

        - ``word`` -- a vector over the ``base_field`` of ``self``.

        - ``sub_field`` -- (default: ``None``) a sub field of the
          ``base_field`` of ``self``; if not specified, it is the prime
          subfield of `\GF{p}` the ``base_field`` of ``self``.

        - ``basis`` -- (default: ``None``) a basis of ``base_field`` of
          ``self`` as a vector space over ``sub_field``. If not specified,
          the default basis is `1,\beta,\ldots,\beta^{m-1}`, where `\beta` is
          the generator of the multiplicative group of `\GF{q^m}` given by Sage.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: a = GF(64).gen()
            sage: c = vector([a^4 + a^3 + 1, a^4 + a^3 + 1, a^4 + a^3 + a^2 + 1])
            sage: c in C
            True
            sage: C.rank_support_of_vector(c)
            Vector space of degree 6 and dimension 2 over Finite Field of size 2
            Basis matrix:
            [1 0 0 1 1 0]
            [0 0 1 0 0 0]

        An example with a non canonical basis::

            sage: K.<a> = GF(2^3)
            sage: G = Matrix(K, [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G)
            sage: c = vector([a^2, a^2, 0])
            sage: basis = [a, a+1, a^2]
            sage: C.rank_support_of_vector(c, basis=basis)
            Vector space of degree 3 and dimension 1 over Finite Field of size 2
            Basis matrix:
            [0 0 1]

        TESTS::

            sage: C.rank_support_of_vector(a)
            Traceback (most recent call last):
            ...
            TypeError: input must be a vector

            sage: C.rank_support_of_vector(c, GF(2^4))
            Traceback (most recent call last):
            ...
            TypeError: the input subfield Finite Field in z4 of size 2^4 is not a subfield of Finite Field in a of size 2^3
        """
        if not isinstance(word, Vector):
            raise TypeError("input must be a vector")
        if sub_field is not None:
            if self.base_field().degree() % sub_field.degree() != 0:
                raise TypeError(f"the input subfield {sub_field} is not a subfield of {self.base_field()}")
        return to_matrix_representation(word, sub_field, basis).column_module()

    def matrix_form_of_vector(self, word):
        r"""
        Return the matrix representation of a word.

        INPUT:

        - ``word`` -- a vector over the ``base_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.matrix_form_of_vector(a)
            [1 1 1]
            [1 1 0]
            [0 0 0]
        """
        return to_matrix_representation(word, self.sub_field())

    def vector_form_of_matrix(self, word):
        r"""
        Return the vector representation of a word.

        INPUT:

        - ``word`` -- a matrix over the ``sub_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
            sage: C.vector_form_of_matrix(m)
            (z6 + 1, z6 + 1, 1)
        """
        return from_matrix_representation(word, self.base_field())


class LinearRankMetricCode(AbstractLinearRankMetricCode):
    r"""
    Linear rank metric codes over a finite field, represented using a
    generator matrix.

    This class should be used for arbitrary and unstructured linear rank metric
    codes. This means that basic operations on the code, such as the computation
    of the minimum distance, will use generic, slow algorithms.

    If you are looking for constructing a code from a more specific family, see
    if the family has been implemented by investigating ``codes.<tab>``. These
    more specific classes use properties particular to that family to allow
    faster algorithms, and could also have family-specific methods.

    EXAMPLES::

        sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
        sage: C = codes.LinearRankMetricCode(G, GF(4))
        sage: C
        [3, 2] linear rank metric code over GF(64)/GF(4)
        sage: C.base_field()
        Finite Field in z6 of size 2^6
        sage: C.sub_field()
        Finite Field in z2 of size 2^2
        sage: C.length()
        3
        sage: C.dimension()
        2
        sage: C[2]
        (z6, z6, 0)
        sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
        sage: word = vector(C.base_field(), [1, 0])
        sage: E(word)
        (1, 1, 0)
    """

    def __init__(self, generator, sub_field=None, basis=None):
        r"""
        See the docstring for :meth:`LinearRankMetricCode`.

        INPUT:

        - ``generator`` -- a generator matrix over the ``base_field`` with
          dimension `k \times n`, where `k` is the dimension of the code and
          `n` its length; or a code over a finite field

        - ``sub_field`` -- (default: ``None``) the sub field of ``self``, if not
          specified, it is the prime field of ``base_field``

        - ``basis`` -- (default: ``None``) a basis of `\GF{q^m}` as a vector space over
          ``sub_field``. If not specified, the default basis is
          `1,\beta,\ldots,\beta^{m-1}` where `\beta` is the generator of the
          multiplicative group of `\GF{q^m}` given by Sage.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4)) # indirect doctest
            sage: C
            [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        base_field = generator.base_ring()
        if not base_field.is_field():
            raise ValueError("'generator' must be defined on a field (not a ring)")

        if not sub_field:
            sub_field = base_field.prime_subfield()

        try:
            gen_basis = None
            if hasattr(generator, "nrows"):  # generator matrix case
                if generator.rank() < generator.nrows():
                    gen_basis = generator.row_space().basis()
            else:
                gen_basis = generator.basis()  # vector space etc. case
            if gen_basis is not None:
                from sage.matrix.constructor import matrix
                generator = matrix(base_field, gen_basis)
                if generator.nrows() == 0:
                    raise ValueError("this linear code contains no nonzero vector")
        except AttributeError:
            # Assume input is an AbstractLinearRankMetricCode, extract its generator matrix
            generator = generator.generator_matrix()

        self._generator_matrix = generator
        self._length = generator.ncols()
        super().__init__(base_field, sub_field, self._length,
                         "GeneratorMatrix", "NearestNeighbor", basis)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C
            [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        R = self.base_field()
        S = self.sub_field()
        if R and S in Fields():
            return "[%s, %s] linear rank metric code over GF(%s)/GF(%s)" % (self.length(), self.dimension(), R.cardinality(), S.cardinality())
        else:
            return "[%s, %s] linear rank metric code over %s/%s" % (self.length(), self.dimension(), R, S)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: latex(C)
            [3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}/\Bold{F}_{2^{2}}
        """
        return "[%s, %s]\\textnormal{ Linear rank metric code over }%s/%s"\
                % (self.length(), self.dimension(), self.base_field()._latex_(), self.sub_field()._latex_())

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Return a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. ``self._generator_matrix``
          will be returned if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.generator_matrix()
            [1 1 0]
            [0 0 1]
        """
        if encoder_name is None or encoder_name == 'GeneratorMatrix':
            g = self._generator_matrix
        else:
            g = super().generator_matrix(encoder_name, **kwargs)
        g.set_immutable()
        return g


####################### decoders ###############################
class LinearRankMetricCodeNearestNeighborDecoder(Decoder):
    r"""
    Construct a decoder for Linear Rank Metric Codes.

    This decoder will decode to the nearest codeword found.
    """

    def __init__(self, code):
        r"""

        INPUT:

        - ``code`` -- a code associated to this decoder

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        super().__init__(code, code.ambient_space(), code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Test equality between LinearRankMetricCodeNearestNeighborDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D1 = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D2 = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D1 == D2
            True
        """
        return isinstance(other, LinearRankMetricCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: latex(D)
            \textnormal{Nearest neighbor decoder for }[3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}/\Bold{F}_{2^{2}}
        """
        return "\\textnormal{Nearest neighbor decoder for }%s" % self.code()._latex_()

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT: a vector of ``self``'s message space

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: G = Matrix(F, [[1,1,0]])
            sage: C = codes.LinearRankMetricCode(G, GF(2))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decode_to_code(vector(F, [a, a, 1]))
            (a, a, 0)
        """
        C = self.code()
        c_min = C.zero()
        h_min = C.rank_weight_of_vector(r)
        for c in C:
            if C.rank_weight_of_vector(c-r) < h_min:
                h_min = C.rank_weight_of_vector(c-r)
                c_min = c
        c_min.set_immutable()
        return c_min

    def decoding_radius(self):
        r"""
        Return maximal number of errors ``self`` can decode.

        EXAMPLES::

            sage: F.<a> = GF(8)
            sage: G = Matrix(F, [[1,a,a^2,0]])
            sage: C = codes.LinearRankMetricCode(G, GF(2))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decoding_radius()
            1
        """
        return (self.code().minimum_distance() - 1) // 2


####################### registration ###############################

LinearRankMetricCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder
