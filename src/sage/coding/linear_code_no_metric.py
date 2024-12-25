# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Generic structures for linear codes of any metric

Class supporting methods available for linear codes over any metric (Hamming,
rank).
"""

from copy import copy

from sage.coding.abstract_code import AbstractCode
from sage.modules.module import Module
from sage.categories.modules import Modules
from sage.modules.free_module import VectorSpace
from sage.coding.encoder import Encoder
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.rings.integer_ring import ZZ


class AbstractLinearCodeNoMetric(AbstractCode, Module):
    r"""
    Abstract class for linear codes of any metric.

    This class contains all the methods that can be used on any linear code
    of any metric. Every abstract class of linear codes over some metric (e.g.
    abstract class for linear codes over the Hamming metric,
    :class:`sage.coding.linear_code.AbstractLinearCode`) should inherit from
    this class.

    To create a new class of linear codes over some metrics, you need to:

    - inherit from AbstractLinearCodeNoMetric

    - call AbstractCode ``__init__`` method in the subclass constructor.
      Example: ``super().__init__(length, "EncoderName",
      "DecoderName", "metric")``.

    - add the following two lines on the class level::

          _registered_encoders = {}
          _registered_decoders = {}


    - fill the dictionary of its encoders in ``sage.coding.__init__.py`` file.
      Example: I want to link the encoder ``MyEncoderClass`` to ``MyNewCodeClass``
      under the name ``MyEncoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_encoders["NameOfMyEncoder"] = MyEncoderClass``
      and all instances of ``MyNewCodeClass`` will be able to use instances of
      ``MyEncoderClass``.

    - fill the dictionary of its decoders in ``sage.coding.__init__`` file.
      Example: I want to link the encoder ``MyDecoderClass`` to ``MyNewCodeClass``
      under the name ``MyDecoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_decoders["NameOfMyDecoder"] = MyDecoderClass``
      and all instances of ``MyNewCodeClass`` will be able to use instances of
      ``MyDecoderClass``.

    - create a generic constructor representative of you abstract class. This
      generic constructor is a class for unstructured linear codes given by some
      generator and considered over the given metric. A good example of this is
      :class:`sage.coding.linear_code.LinearCode`, which is a generic constructor
      for :class:`sage.coding.linear_code.AbstractLinearCode`, an abstract class
      for linear codes over the Hamming metric.

    - set a private field in the ``__init__`` method specifying the generic
      constructor, (e.g. ``MyAbstractCode._generic_constructor = MyCode``)


    It is assumed that the subclass codes are linear over ``base_field``. To
    test this, it is recommended to add a test suite test to the generic
    constructor. To do this, create a representative of your code `C` and run
    ``TestSuite(C).run()``. A good example of this is in
    :class:`sage.coding.linear_code.LinearCode`.


    As the class :class:`AbstractLinearCodeNoMetric` is not designed to be
    instantiated, it does not have any representation methods. You should
    implement ``_repr_`` and ``_latex_`` methods in the subclass.

    .. WARNING::

        A lot of methods of the abstract class rely on the knowledge of a generator matrix.
        It is thus strongly recommended to set an encoder with a generator matrix implemented
        as a default encoder.

    TESTS:

    If the name of the default decoder is not known by the class, it will raise
    a exception::

        sage: from sage.coding.linear_code_no_metric import AbstractLinearCodeNoMetric
        sage: class MyCodeFamily(AbstractLinearCodeNoMetric):
        ....:   def __init__(self, field, length, dimension, generator_matrix):
        ....:       AbstractLinearCodeNoMetric.__init__(self, field, length, "Systematic", "Fail", "MyMetric")
        ....:       self._dimension = dimension
        ....:       self._generator_matrix = generator_matrix
        ....:   def generator_matrix(self):
        ....:       return self._generator_matrix
        ....:   def _repr_(self):
        ....:       return "[%d, %d] dummy code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

        sage: generator_matrix = matrix(GF(17), 5, 10,
        ....:                           {(i,i):1 for i in range(5)})
        sage: C = MyCodeFamily(GF(17), 10, 5, generator_matrix)
        Traceback (most recent call last):
        ...
        ValueError: You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders

    If the name of the default encoder is not known by the class, it will raise
    an exception::

        sage: class MyCodeFamily2(sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric):
        ....:   def __init__(self, field, length, dimension, generator_matrix):
        ....:       sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric.__init__(self, field, length, "Fail", "Syndrome", "MyMetric")
        ....:       self._dimension = dimension
        ....:       self._generator_matrix = generator_matrix
        ....:   def generator_matrix(self):
        ....:       return self._generator_matrix
        ....:   def _repr_(self):
        ....:       return "[%d, %d] dummy code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

        sage: C = MyCodeFamily2(GF(17), 10, 5, generator_matrix)
        Traceback (most recent call last):
        ...
        ValueError: You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders

    A ring instead of a field::

        sage: MyCodeFamily2(IntegerModRing(4), 4, 4, matrix.ones(4))
        Traceback (most recent call last):
        ...
        ValueError: 'base_field' must be a field (and Ring of integers modulo 4 is not one)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, default_encoder_name, default_decoder_name, metric='Hamming'):
        """
        Initialize mandatory parameters that any linear code shares.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear code. The class
        :class:`sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric`
        should never be directly instantiated.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer, must be > 0)

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``

        - ``metric`` -- (default: ``Hamming``) the metric of ``self``

        EXAMPLES:

            sage: from sage.coding.linear_code_no_metric import AbstractLinearCodeNoMetric
            sage: from sage.coding.linear_code import LinearCodeSyndromeDecoder
            sage: class MyLinearCode(AbstractLinearCodeNoMetric):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       self._registered_decoders['Syndrome'] = LinearCodeSyndromeDecoder
            ....:       AbstractLinearCodeNoMetric.__init__(self, field, length, "Systematic", "Syndrome")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())
            sage: C = MyLinearCode(GF(2), 1, 1, matrix(GF(2), [1]))
        """

        self._registered_encoders['Systematic'] = LinearCodeSystematicEncoder

        if not base_field.is_field():
            raise ValueError("'base_field' must be a field (and {} is not one)".format(base_field))
        if default_encoder_name not in self._registered_encoders:
            raise ValueError("You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders")
        if default_decoder_name not in self._registered_decoders:
            raise ValueError("You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders")

        # if not self.dimension() <= length:
        #     raise ValueError("The dimension of the code can be at most its length, {}".format(length))

        super().__init__(length, default_encoder_name, default_decoder_name, metric)
        cat = Modules(base_field).FiniteDimensional().WithBasis().Finite()
        facade_for = VectorSpace(base_field, self._length)
        self.Element = facade_for.Element
        Parent.__init__(self, base=base_field, facade=facade_for, category=cat)

    def base_field(self):
        r"""
        Return the base field of ``self``.

        EXAMPLES::

            sage: G  = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                     [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C.base_field()
            Finite Field of size 2
        """
        return self.base_ring()

    def ambient_space(self):
        r"""
        Return the ambient vector space of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.ambient_space()
            Vector space of dimension 7 over Finite Field of size 2
        """
        return VectorSpace(self.base_ring(),self.length())

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Return a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. The default encoder of ``self``
          will be used if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used

        EXAMPLES::

            sage: G = matrix(GF(3), 2, [1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix()
            [1 2 1]
            [2 1 1]
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.generator_matrix()

    def __eq__(self, other):
        r"""
        Test equality between two linear codes.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C1 = LinearCode(G)
            sage: C1 == 5
            False
            sage: C2 = LinearCode(G)
            sage: C1 == C2
            True
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,1,1]])
            sage: C2 = LinearCode(G)
            sage: C1 == C2
            False
            sage: G = Matrix(GF(3), [[1,2,1,0,0,0,0]])
            sage: C3 = LinearCode(G)
            sage: C1 == C3
            False
        """
        # Fail without computing the generator matrix if possible:
        if not (isinstance(other, AbstractLinearCodeNoMetric)
                and self.length() == other.length()
                and self.dimension() == other.dimension()
                and self.base_ring() == other.base_ring()):
            return False
        # Check that basis elements of `other` are all in `self.`
        # Since we're over a field and since the dimensions match, the codes
        # must be equal.
        # This implementation may avoid linear algebra altogether, if `self`
        # implements an efficient way to obtain a parity check matrix, and in
        # the worst case does only one system solving.
        return all(c in self for c in other.gens())

    def __ne__(self, other):
        r"""
        Test inequality of ``self`` and ``other``.

        This is a generic implementation, which returns the inverse of
        ``__eq__`` for ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C1 = LinearCode(G)
            sage: C2 = LinearCode(G)
            sage: C1 != C2
            False
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,1,1]])
            sage: C2 = LinearCode(G)
            sage: C1 != C2
            True
        """
        return not self == other

    def dimension(self):
        r"""
        Return the dimension of this code.

        EXAMPLES::

            sage: G = matrix(GF(2), [[1,0,0], [1,1,0]])
            sage: C = LinearCode(G)
            sage: C.dimension()
            2

        TESTS:

        Check that :issue:`21156` is fixed::

            sage: from sage.coding.linear_code import AbstractLinearCode
            sage: from sage.coding.encoder import Encoder
            sage: class MonkeyCode(AbstractLinearCode):
            ....:     _registered_encoders = {}
            ....:     _registered_decoders = {}
            ....:     def __init__(self):
            ....:         super().__init__(GF(5), 10, "Monkey", "Syndrome")
            ....:
            sage: class MonkeyEncoder(Encoder):
            ....:     def __init__(self, code):
            ....:         super().__init__(C)
            ....:     @cached_method
            ....:     def generator_matrix(self):
            ....:         G = identity_matrix(GF(5), 5).augment(matrix(GF(5), 5, 7))
            ....:         return G
            ....:
            sage: MonkeyCode._registered_encoders["Monkey"] = MonkeyEncoder
            sage: C = MonkeyCode()
            sage: C.dimension()
            5
        """
        try:
            return self._dimension
        except AttributeError:
            dimension = self.generator_matrix().nrows()
            self._dimension = dimension
            return self._dimension

    def cardinality(self):
        r"""
        Return the size of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.cardinality()
            16
            sage: len(C)
            16
        """
        return self.base_ring().order()**self.dimension()

    __len__ = cardinality

    def rate(self):
        r"""
        Return the ratio of the number of information symbols to
        the code length.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.rate()
            4/7
        """
        return self.dimension() / self.length()

    @cached_method
    def gens(self) -> list:
        r"""
        Return the generators of this code as a list of vectors.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.gens()
             [(1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1),
              (0, 0, 1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 1, 1)]
        """
        return self.generator_matrix().rows()

    def basis(self):
        r"""
        Return a basis of ``self``.

        OUTPUT: ``Sequence`` -- an immutable sequence whose universe is ambient space of ``self``

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.basis()
            [
            (1, 0, 0, 0, 0, 1, 1),
            (0, 1, 0, 0, 1, 0, 1),
            (0, 0, 1, 0, 1, 1, 0),
            (0, 0, 0, 1, 1, 1, 1)
            ]
            sage: C.basis().universe()
            Vector space of dimension 7 over Finite Field of size 2
        """
        gens = self.gens()
        from sage.structure.sequence import Sequence
        return Sequence(gens, universe=self.ambient_space(), check=False, immutable=True, cr=True)

    @cached_method
    def parity_check_matrix(self):
        r"""
        Return the parity check matrix of ``self``.

        The parity check matrix of a linear code `C` corresponds to the
        generator matrix of the dual code of `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Cperp = C.dual_code()
            sage: C; Cperp
            [7, 4] Hamming Code over GF(2)
            [7, 3] linear code over GF(2)
            sage: C.generator_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C.parity_check_matrix()
             [1 0 1 0 1 0 1]
             [0 1 1 0 0 1 1]
             [0 0 0 1 1 1 1]
            sage: Cperp.parity_check_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: Cperp.generator_matrix()
             [1 0 1 0 1 0 1]
             [0 1 1 0 0 1 1]
             [0 0 0 1 1 1 1]
        """
        G = self.generator_matrix()
        H = G.right_kernel()
        M = H.basis_matrix()
        M.set_immutable()
        return M

    def syndrome(self, r):
        r"""
        Return the syndrome of ``r``.

        The syndrome of ``r`` is the result of `H \times r` where `H` is
        the parity check matrix of ``self``. If ``r`` belongs to ``self``,
        its syndrome equals to the zero vector.

        INPUT:

        - ``r`` -- a vector of the same length as ``self``

        OUTPUT: a column vector

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: r = vector(GF(2), (1,0,1,0,1,0,1))
            sage: r in C
            True
            sage: C.syndrome(r)
            (0, 0, 0)

        If ``r`` is not a codeword, its syndrome is not equal to zero::

            sage: r = vector(GF(2), (1,0,1,0,1,1,1))
            sage: r in C
            False
            sage: C.syndrome(r)
            (0, 1, 1)

        Syndrome computation works fine on bigger fields::

            sage: C = codes.random_linear_code(GF(59), 12, 4)
            sage: c = C.random_element()
            sage: C.syndrome(c)
            (0, 0, 0, 0, 0, 0, 0, 0)
        """
        return self.parity_check_matrix()*r

    def __contains__(self, v):
        r"""
        Return ``True`` if `v` can be coerced into ``self``.
        Otherwise, return ``False``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: vector((1, 0, 0, 0, 0, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 2, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 0, 1/2, 1)) in C # indirect doctest
            False
        """
        if v not in self.ambient_space() or len(v) != self.length():
            return False
        return self.syndrome(v) == 0

    def systematic_generator_matrix(self, systematic_positions=None):
        """
        Return a systematic generator matrix of the code.

        A generator matrix of a code is called systematic if it contains
        a set of columns forming an identity matrix.

        INPUT:

        - ``systematic_positions`` -- (default: ``None``) if supplied, the set
          of systematic positions in the systematic generator matrix. See the
          documentation for :class:`LinearCodeSystematicEncoder` details.

        EXAMPLES::

            sage: G = matrix(GF(3), [[ 1, 2, 1, 0],\
                                     [ 2, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: C.generator_matrix()
            [1 2 1 0]
            [2 1 1 1]
            sage: C.systematic_generator_matrix()
            [1 2 0 1]
            [0 0 1 2]

        Specific systematic positions can also be requested::

            sage: C.systematic_generator_matrix(systematic_positions=[3,2])
            [1 2 0 1]
            [1 2 1 0]
        """
        systematic_positions = tuple(systematic_positions) if systematic_positions else None
        return self.encoder("Systematic", systematic_positions=systematic_positions).generator_matrix()

    def standard_form(self, return_permutation=True):
        r"""
        Return a linear code which is permutation-equivalent to ``self`` and
        admits a generator matrix in standard form.

        A generator matrix is in standard form if it is of the form `[I \vert
        A]`, where `I` is the `k \times k` identity matrix. Any code admits a
        generator matrix in systematic form, i.e. where a subset of the columns
        form the identity matrix, but one might need to permute columns to allow
        the identity matrix to be leading.

        INPUT:

        - ``return_permutation`` -- boolean (default: ``True``); if ``True``,
          the column permutation which brings ``self`` into the returned code
          is also returned

        OUTPUT: a :class:`LinearCode` whose :meth:`systematic_generator_matrix`
        is guaranteed to be of the form `[I \vert A]`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
            sage: Cs,p = C.standard_form()
            sage: p
            []
            sage: Cs is C
            True
            sage: C = LinearCode(matrix(GF(2), [[1,0,0,0,1,1,0],\
                                                [0,1,0,1,0,1,0],\
                                                [0,0,0,0,0,0,1]]))
            sage: Cs, p = C.standard_form()
            sage: p
            [1, 2, 7, 3, 4, 5, 6]
            sage: Cs.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 0 0 0]
        """
        E = self.encoder("Systematic")
        if E.systematic_positions() == tuple(range(self.dimension())):
            from sage.combinat.permutation import Permutation
            return self, Permutation([])
        else:
            perm = E.systematic_permutation()
            return self.permuted_code(perm), perm

    def redundancy_matrix(self):
        r"""
        Return the non-identity columns of a systematic generator matrix for
        ``self``.

        A systematic generator matrix is a generator matrix such that a subset
        of its columns forms the identity matrix. This method returns the
        remaining part of the matrix.

        For any given code, there can be many systematic generator matrices
        (depending on which positions should form the identity). This method
        will use the matrix returned by
        :meth:`AbstractLinearCode.systematic_generator_matrix`.

        OUTPUT:

        - An `k \times (n-k)` matrix.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.generator_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C.redundancy_matrix()
             [0 1 1]
             [1 0 1]
             [1 1 0]
             [1 1 1]
            sage: C = LinearCode(matrix(GF(3),2,[1,2,0,\
                                                 2,1,1]))
            sage: C.systematic_generator_matrix()
            [1 2 0]
            [0 0 1]
            sage: C.redundancy_matrix()
            [2]
            [0]
        """
        E = self.encoder("Systematic")
        G = E.generator_matrix()
        return G.delete_columns(E.systematic_positions())

    @cached_method
    def information_set(self):
        """
        Return an information set of the code.

        Return value of this method is cached.

        A set of column positions of a generator matrix of a code
        is called an information set if the corresponding columns
        form a square matrix of full rank.

        OUTPUT: information set of a systematic generator matrix of the code

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,2,0,
            ....:                     2,1,1])
            sage: code = LinearCode(G)
            sage: code.systematic_generator_matrix()
            [1 2 0]
            [0 0 1]
            sage: code.information_set()
            (0, 2)
        """
        return self.encoder("Systematic").systematic_positions()

    def is_information_set(self, positions):
        """
        Return whether the given positions form an information set.

        INPUT:

        - A list of positions, i.e. integers in the range 0 to `n-1` where `n`
          is the length of ``self``.

        OUTPUT: boolean indicating whether the positions form an information set

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,2,0,
            ....:                     2,1,1])
            sage: code = LinearCode(G)
            sage: code.is_information_set([0,1])
            False
            sage: code.is_information_set([0,2])
            True
        """
        try:
            self.encoder("Systematic", systematic_positions=tuple(positions))
            return True
        except ValueError:
            return False

    def __iter__(self):
        """
        Return an iterator over the elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: [list(c) for c in C if c.hamming_weight() < 4]
            [[0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 1],
             [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0],
             [1, 1, 1, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0, 1, 0], [0, 0, 1, 1, 0, 0, 1]]

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: L = list(C)
            sage: L[10].is_immutable()
            True
        """
        from sage.modules.finite_submodule_iter import \
                                                FiniteFieldsubspace_iterator
        return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)

    def __getitem__(self, i):
        r"""
        Return the `i`-th codeword of this code.

        The implementation of this depends on the implementation of the
        :meth:`.__iter__` method.

        The implementation is as follows. Suppose that:

        - the primitive element of the base_ring of this code is `a`,
        - the prime subfield is `p`,
        - the field has order `p^m`,
        - the code has dimension `k`,
        - and the generator matrix is `G`.

        Then the :meth:`.__iter__` method returns the elements in this order:

        1. first, the following ordered list is returned:
           ``[i*a^0 * G[0] for i in range(p)]``
        2. Next, the following ordered list is returned:
           ``[i*a^0 * G[0] + a^1*G[0] for i in range(p)]``
        3. This continues till we get
           ``[(i*a^0 +(p-1)*a^1 +...+ (p-1)*a^(m-1))*G[0] for i in range(p)]``
        4. Then, we move to G[1]:
           ``[i*a^0 * G[0] + a^0*G[1] for i in range(p)]``,

        and so on.
        Hence the `i`-th element can be obtained by the `p`-adic expansion
        of `i` as ``[i_0, i_1, ...,i_{m-1}, i_m, i_{m+1}, ..., i_{km-1}].``

        The element that is generated is:

        .. MATH::

             \begin{aligned}
             & (i_0 a^0 + i_1 a^1 + \cdots + i_{m-1} a^{m-1}) G[0] + \\
             & (i_m a^0 + i_{m+1} a^1 + \cdots + i_{2m-1} a^{m-1}) G[1] + \\
             & \vdots\\
             & (i_{(k-1)m} a^0 + \cdots + i_{km-1} a^{m-1}) G[k-1]
             \end{aligned}

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C[24]
            (2, 2, 0, 1, 2, 2, 0)
            sage: C[24] == C.list()[24]
            True

        TESTS::

            sage: C = random_matrix(GF(25,'a'), 2, 7).row_space()
            sage: C = LinearCode(C.basis_matrix())
            sage: Clist = C.list()
            sage: all(C[i] == Clist[i] for i in range(len(C)))
            True

        Check that only the indices less than the size of the code are
        allowed::

            sage: C[25**2]
            Traceback (most recent call last):
            ...
            IndexError: The value of the index 'i' (=625) must be between
            0 and 'q^k -1' (=624), inclusive, where 'q' is the size of the
            base field and 'k' is the dimension of the code.

        Check that codewords are immutable. See :issue:`16338`::

            sage: C[0].is_immutable()
            True
        """
        # IMPORTANT: If the __iter__() function implementation is changed
        # then the implementation here must also be changed so that
        # list(self)[i] and self[i] both return the same element.

        F = self.base_ring()
        maxindex = F.order()**self.dimension()-1
        if i < 0 or i > maxindex:
            raise IndexError("The value of the index 'i' (={}) must be between "
                             "0 and 'q^k -1' (={}), inclusive, where 'q' is "
                             "the size of the base field and 'k' is the "
                             "dimension of the code.".format(i, maxindex))

        a = F.primitive_element()
        m = F.degree()
        p = F.prime_subfield().order()
        A = [a ** k for k in range(m)]
        G = self.generator_matrix()
        N = self.dimension()*F.degree() # the total length of p-adic vector
        ivec = Integer(i).digits(p, padto=N)

        codeword = 0
        row = 0
        for g in G:
            codeword += sum(ivec[j+row*m]*A[j] for j in range(m)) * g
            row += 1

        # The codewords for a specific code can not change. So, we set them
        # to be immutable.
        codeword.set_immutable()
        return codeword

    def __hash__(self):
        r"""
        Return the hash value of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: hash(C) #random
            9015017528451745710

        If ``C1`` and ``C2`` are two codes which only differ by the
        coefficients of their generator matrices, their hashes are
        different (we check that the bug found in :issue:`18813` is
        fixed)::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C1 = LinearCode(G)
            sage: G = Matrix(GF(2), [[1,0,0,1,0,1,0],[0,1,0,0,1,0,0],[0,0,1,1,0,1,0],[0,0,0,0,0,0,1]])
            sage: C2 = LinearCode(G)
            sage: hash(C1) != hash(C2)
            True
        """
        Str = str(self)
        G = self.generator_matrix()
        return hash((Str, G)) ^ hash(Str) ^ hash(G)

    def is_subcode(self, other):
        """
        Return ``True`` if ``self`` is a subcode of ``other``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: G1 = C1.generator_matrix()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True
            sage: C1.is_subcode(C2)
            False
            sage: C3 = C1.extended_code()
            sage: C1.is_subcode(C3)
            False
            sage: C4 = C1.punctured([1])
            sage: C4.is_subcode(C1)
            False
            sage: C5 = C1.shortened([1])
            sage: C5.is_subcode(C1)
            False
            sage: C1 = codes.HammingCode(GF(9,"z"), 3)
            sage: G1 = C1.generator_matrix()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True
        """
        G = self.generator_matrix()
        return all(r in other for r in G.rows())

    def is_permutation_automorphism(self, g):
        r"""
        Return `1` if `g` is an element of `S_n` (`n` = length of ``self``) and
        if `g` is an automorphism of ``self``.

        EXAMPLES::

            sage: # needs sage.groups
            sage: C = codes.HammingCode(GF(3), 3)
            sage: g = SymmetricGroup(13).random_element()
            sage: C.is_permutation_automorphism(g)
            0
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0], [0,1,1,1,0,0,0,0],
            ....:          [0,0,0,0,0,0,0,1], [0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: S8 = SymmetricGroup(8)
            sage: g = S8("(2,3)")
            sage: C.is_permutation_automorphism(g)
            1
            sage: g = S8("(1,2,3,4)")
            sage: C.is_permutation_automorphism(g)
            0
        """
        basis = self.generator_matrix().rows()
        H = self.parity_check_matrix()
        V = H.column_space()
        HGm = H*g.matrix()
        for c in basis:
            if HGm*c != V(0):
                return False
        return True

    def permuted_code(self, p):
        r"""
        Return the permuted code, which is equivalent to ``self`` via the
        column permutation ``p``.

        EXAMPLES::

            sage: # needs sage.groups
            sage: C = codes.HammingCode(GF(2), 3)
            sage: G = C.permutation_automorphism_group(); G
            Permutation Group with generators
             [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: g = G("(2,3)(6,7)")
            sage: Cg = C.permuted_code(g)
            sage: Cg
            [7, 4] linear code over GF(2)
            sage: C.generator_matrix() == Cg.systematic_generator_matrix()
            True
        """
        if not hasattr(self, "_generic_constructor"):
            raise NotImplementedError("Generic constructor not set for the class of codes")
        G = copy(self.generator_matrix())
        G.permute_columns(p)
        return self._generic_constructor(G)

    def dual_code(self):
        r"""
        Return the dual code `C^{\perp}` of the code `C`,

        .. MATH::

            C^{\perp} = \{ v \in V\ |\ v\cdot c = 0,\ \forall c \in C \}.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.dual_code()
            [7, 3] linear code over GF(2)
            sage: C = codes.HammingCode(GF(4, 'a'), 3)
            sage: C.dual_code()
            [21, 3] linear code over GF(4)
        """
        if not hasattr(self, "_generic_constructor"):
            raise NotImplementedError("Generic constructor not set for the class of codes")
        return self._generic_constructor(self.parity_check_matrix())

    def is_self_dual(self) -> bool:
        """
        Return ``True`` if the code is self-dual (in the usual Hamming inner
        product) and ``False`` otherwise.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: C.is_self_dual()
            True
            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.is_self_dual()
            False
        """
        return self == self.dual_code()

    def is_self_orthogonal(self):
        """
        Return ``True`` if this code is self-orthogonal and ``False``
        otherwise.

        A code is self-orthogonal if it is a subcode of its dual.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: C.is_self_orthogonal()
            True
            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.is_self_orthogonal()
            False
            sage: C = codes.QuasiQuadraticResidueCode(11)  # optional - gap_package_guava
            sage: C.is_self_orthogonal()                   # optional - gap_package_guava
            True
        """
        return self.is_subcode(self.dual_code())

    @cached_method
    def zero(self):
        r"""
        Return the zero vector of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.zero()
            (0, 0, 0, 0, 0, 0, 0)
            sage: C.sum(()) # indirect doctest
            (0, 0, 0, 0, 0, 0, 0)
            sage: C.sum((C.gens())) # indirect doctest
            (1, 1, 1, 1, 1, 1, 1)
        """
        return self.ambient_space().zero()


####################### encoders ###############################


class LinearCodeSystematicEncoder(Encoder):
    r"""
    Encoder based on a generator matrix in systematic form for Linear codes.

    To encode an element of its message space, this encoder first builds a
    generator matrix in systematic form. What is called systematic form here
    is the reduced row echelon form of a matrix, which is not necessarily
    `[I \vert H]`, where `I` is the identity block and `H` the parity block.
    One can refer to :meth:`LinearCodeSystematicEncoder.generator_matrix`
    for a concrete example.
    Once such a matrix has been computed, it is used to encode any message
    into a codeword.

    This encoder can also serve as the default encoder of a code defined by a
    parity check matrix: if the :class:`LinearCodeSystematicEncoder` detects
    that it is the default encoder, it computes a generator matrix as the
    reduced row echelon form of the right kernel of the parity check matrix.

    INPUT:

    - ``code`` -- the associated code of this encoder

    - ``systematic_positions`` -- (default: ``None``) the positions in codewords that
      should correspond to the message symbols. A list of `k` distinct integers in
      the range 0 to `n-1` where `n` is the length of the code and `k` its
      dimension. The `0`-th symbol of a message will then be at position
      ``systematic_positions[0]``, the 1st index at position
      ``systematic_positions[1]``, etc. A :exc:`ValueError` is raised at
      construction time if the supplied indices do not form an information set.

    EXAMPLES:

    The following demonstrates the basic usage of :class:`LinearCodeSystematicEncoder`::

        sage: LCSE = codes.encoders.LinearCodeSystematicEncoder
        sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0,0],\
                                 [1,0,0,1,1,0,0,0],\
                                 [0,1,0,1,0,1,0,0],\
                                 [1,1,0,1,0,0,1,1]])
        sage: C = LinearCode(G)
        sage: E = LCSE(C)
        sage: E.generator_matrix()
        [1 0 0 0 0 1 1 1]
        [0 1 0 0 1 0 1 1]
        [0 0 1 0 1 1 0 0]
        [0 0 0 1 1 1 1 1]
        sage: E2 = LCSE(C, systematic_positions=[5,4,3,2])
        sage: E2.generator_matrix()
        [1 0 0 0 0 1 1 1]
        [0 1 0 0 1 0 1 1]
        [1 1 0 1 0 0 1 1]
        [1 1 1 0 0 0 0 0]

    An error is raised if one specifies systematic positions which do not form
    an information set::

        sage: E3 = LCSE(C, systematic_positions=[0,1,6,7])
        Traceback (most recent call last):
        ...
        ValueError: systematic_positions are not an information set


    We exemplify how to use :class:`LinearCodeSystematicEncoder` as the default
    encoder. The following class is the dual of the repetition code::

        sage: class DualRepetitionCode(sage.coding.linear_code.AbstractLinearCode):
        ....:   def __init__(self, field, length):
        ....:       super().__init__(field, length, "Systematic", "Syndrome")
        ....:
        ....:   def parity_check_matrix(self):
        ....:       return Matrix(self.base_field(), [1]*self.length())
        ....:
        ....:   def _repr_(self):
        ....:       return "Dual of the [%d, 1] Repetition Code over GF(%s)" % (self.length(), self.base_field().cardinality())
        ....:
        sage: DualRepetitionCode(GF(3), 5).generator_matrix()
        [1 0 0 0 2]
        [0 1 0 0 2]
        [0 0 1 0 2]
        [0 0 0 1 2]


    An exception is thrown if :class:`LinearCodeSystematicEncoder` is the default encoder but no
    parity check matrix has been specified for the code::

        sage: class BadCodeFamily(sage.coding.linear_code.AbstractLinearCode):
        ....:   def __init__(self, field, length):
        ....:       super().__init__(field, length, "Systematic", "Syndrome")
        ....:
        ....:   def _repr_(self):
        ....:       return "I am a badly defined code"
        ....:
        sage: BadCodeFamily(GF(3), 5).generator_matrix()
        Traceback (most recent call last):
        ...
        ValueError: a parity check matrix must be specified
        if LinearCodeSystematicEncoder is the default encoder
    """

    def __init__(self, code, systematic_positions=None):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: E
            Systematic encoder for [7, 4] linear code over GF(2)
        """
        super().__init__(code)
        self._systematic_positions = tuple(systematic_positions) if systematic_positions else None
        if systematic_positions:
            # Test that systematic_positions consists of integers in the right
            # range. We test that len(systematic_positions) = code.dimension()
            # in self.generator_matrix() to avoid possible infinite recursion.
            if (not all( e in ZZ and e >= 0 and e < code.length() for e in systematic_positions)) \
               or len(systematic_positions) != len(set(systematic_positions)):
                raise ValueError("systematic positions must be a tuple of distinct integers in the range 0 to n-1 where n is the length of the code")
            # Test that the systematic positions are an information set
            self.generator_matrix()

    def __eq__(self, other):
        r"""
        Test equality between LinearCodeSystematicEncoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: E1 = codes.encoders.LinearCodeSystematicEncoder(LinearCode(G))
            sage: E2 = codes.encoders.LinearCodeSystematicEncoder(LinearCode(G))
            sage: E1 == E2
            True
            sage: E1.systematic_positions()
            (0, 1, 2)
            sage: E3 = codes.encoders.LinearCodeSystematicEncoder(LinearCode(G), systematic_positions=(2,5,6))
            sage: E3.systematic_positions()
            (2, 5, 6)
            sage: E1 == E3
            False
        """
        return isinstance(other, LinearCodeSystematicEncoder)\
                and self.code() == other.code()\
                and self.systematic_positions() == other.systematic_positions()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: E
            Systematic encoder for [7, 4] linear code over GF(2)
        """
        return "Systematic encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: latex(E)
            \textnormal{Systematic encoder for }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Systematic encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Return a generator matrix in systematic form of the associated code of ``self``.

        Systematic form here means that a subsets of the columns of the matrix
        forms the identity matrix.

        .. NOTE::

            The matrix returned by this method will not necessarily be `[I \vert H]`, where `I`
            is the identity block and `H` the parity block. If one wants to know which columns
            create the identity block, one can call :meth:`systematic_positions`

        EXAMPLES::

            sage: LCSE = codes.encoders.LinearCodeSystematicEncoder
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],\
                                     [1,0,0,1,1,0,0],\
                                     [0,1,0,1,0,1,0],\
                                     [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = LCSE(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]

        We can ask for different systematic positions::

            sage: E2 = LCSE(C, systematic_positions=[5,4,3,2])
            sage: E2.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [1 1 0 1 0 0 1]
            [1 1 1 0 0 0 0]

        Another example where there is no generator matrix of the form `[I \vert H]`::

            sage: G = Matrix(GF(2), [[1,1,0,0,1,0,1],\
                                     [1,1,0,0,1,0,0],\
                                     [0,0,1,0,0,1,0],\
                                     [0,0,1,0,1,0,1]])
            sage: C = LinearCode(G)
            sage: E = LCSE(C)
            sage: E.generator_matrix()
            [1 1 0 0 0 1 0]
            [0 0 1 0 0 1 0]
            [0 0 0 0 1 1 0]
            [0 0 0 0 0 0 1]
        """
        C = self.code()
        # This if statement detects if this encoder is itself the default encoder.
        # In this case, attempt building the generator matrix from the parity
        # check matrix
        if hasattr(self, "_use_pc_matrix"):
            if self._use_pc_matrix == 1:
                self._use_pc_matrix = 2
                return C.parity_check_matrix().right_kernel_matrix()
            else:
                raise ValueError("a parity check matrix must be specified if LinearCodeSystematicEncoder is the default encoder")
        else:
            self._use_pc_matrix = 1
            M = copy(C.generator_matrix())
        if not self._systematic_positions:
            M.echelonize()
        else:
            k = M.nrows() # it is important that k is *not* computed as C.dimension() to avoid possible cyclic dependency
            if len(self._systematic_positions) != k:
                raise ValueError("systematic_positions must be a tuple of length equal to the dimension of the code")
            # Permute the columns of M and bring to reduced row echelon formb
            perm = self.systematic_permutation()
            M.permute_columns(perm)
            M.echelonize()
            if M[:,:k].is_singular():
                raise ValueError("systematic_positions are not an information set")
            M.permute_columns(perm.inverse())
        M.set_immutable()
        return M

    def systematic_permutation(self):
        r"""
        Return a permutation which would take the systematic positions into [0,..,k-1].

        EXAMPLES::

            sage: C = LinearCode(matrix(GF(2), [[1,0,0,0,1,1,0],\
                                                [0,1,0,1,0,1,0],\
                                                [0,0,0,0,0,0,1]]))
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: E.systematic_positions()
            (0, 1, 6)
            sage: E.systematic_permutation()
            [1, 2, 7, 3, 4, 5, 6]
        """
        n = self.code().length()
        systematic_positions = self.systematic_positions()
        k = len(systematic_positions)
        lp = [None] * n
        for (i, j) in zip(range(k), systematic_positions):
            lp[i] = j
        j = k
        set_sys_pos = set(systematic_positions)
        for i in range(n):
            if i not in set_sys_pos:
                lp[j] = i
                j += 1
        from sage.combinat.permutation import Permutation
        return Permutation([1 + e for e in lp])

    def systematic_positions(self):
        r"""
        Return a tuple containing the indices of the columns which form an
        identity matrix when the generator matrix is in systematic form.

        EXAMPLES::

            sage: LCSE = codes.encoders.LinearCodeSystematicEncoder
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],\
                                     [1,0,0,1,1,0,0],\
                                     [0,1,0,1,0,1,0],\
                                     [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = LCSE(C)
            sage: E.systematic_positions()
            (0, 1, 2, 3)

        We take another matrix with a less nice shape::

            sage: G = Matrix(GF(2), [[1,1,0,0,1,0,1],\
                                     [1,1,0,0,1,0,0],\
                                     [0,0,1,0,0,1,0],\
                                     [0,0,1,0,1,0,1]])
            sage: C = LinearCode(G)
            sage: E = LCSE(C)
            sage: E.systematic_positions()
            (0, 2, 4, 6)

        The systematic positions correspond to the positions which carry information in a codeword::

            sage: MS = E.message_space()
            sage: m = MS.random_element()
            sage: c = m * E.generator_matrix()
            sage: pos = E.systematic_positions()
            sage: info = MS([c[i] for i in pos])
            sage: m == info
            True

        When constructing a systematic encoder with specific systematic
        positions, then it is guaranteed that this method returns exactly those
        positions (even if another choice might also be systematic)::

            sage: G = Matrix(GF(2), [[1,0,0,0],\
                                     [0,1,0,0],\
                                     [0,0,1,1]])
            sage: C = LinearCode(G)
            sage: E = LCSE(C, systematic_positions=[0,1,3])
            sage: E.systematic_positions()
            (0, 1, 3)
        """
        return self._systematic_positions if self._systematic_positions else self.generator_matrix().pivots()
