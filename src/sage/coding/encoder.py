r"""
Encoder

Representation of a bijection between a message space and a code.
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modules.free_module_element import vector
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject

class Encoder(SageObject):
    r"""
    Abstract top-class for :class:`Encoder` objects.

    Every encoder class should inherit from this abstract class.

    To implement an encoder, you need to:

    - inherit from :class:`Encoder`,

    - call ``Encoder.__init__`` in the subclass constructor.
      Example: ``super(SubclassName, self).__init__(code)``.
      By doing that, your subclass will have its ``code`` parameter initialized.

    - Then, if the message space is a vector space, default implementations of :meth:`encode` and
      :meth:`unencode_nocheck` methods are provided. These implementations rely on :meth:`generator_matrix`
      which you need to override to use the default implementations.

    - If the message space is not of the form `F^k`, where `F` is a finite field,
      you cannot have a generator matrix.
      In that case, you need to override :meth:`encode`, :meth:`unencode_nocheck` and
      :meth:`message_space`.

    - By default, comparison of :class:`Encoder` (using methods ``__eq__`` and ``__ne__`` ) are
      by memory reference: if you build the same encoder twice, they will be different. If you
      need something more clever, override ``__eq__`` and ``__ne__`` in your subclass.

    - As :class:`Encoder` is not designed to be instantiated, it does not have any representation
      methods. You should implement ``_repr_`` and ``_latex_`` methods in the subclass.

    REFERENCES:

    .. [Nielsen] Johan S. R. Nielsen, (https://bitbucket.org/jsrn/codinglib/)
    """

    def __init__(self, code):
        r"""
        Initializes mandatory parameters for an :class:`Encoder` object.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear code. An abstract
        encoder object should never be created.

        INPUT:

        - ``code`` -- the associated code of ``self``

        EXAMPLES:

        We first create a new :class:`Encoder` subclass::

            sage: from sage.coding.encoder import Encoder
            sage: class EncoderExample(Encoder):
            ....:   def __init__(self, code):
            ....:       super(EncoderExample, self).__init__(code)

        We now create a member of our newly made class::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: E = EncoderExample(C)

        We can check its parameters::

            sage: E.code()
            Linear code of length 4, dimension 2 over Finite Field of size 2
        """
        self._code = code

    def __ne__(self, other):
        r"""
        Tests inequality of ``self`` and ``other``.

        This is a generic implementation, which returns the inverse of ``__eq__`` for self.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: E1 = LinearCode(G).encoder()
            sage: E2 = LinearCode(G).encoder()
            sage: E1 != E2
            False
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,1,1]])
            sage: E2 = LinearCode(G).encoder()
            sage: E1 != E2
            True
        """
        return not self == other

    def encode(self, word):
        r"""
        Transforms an element of the message space into a codeword.

        This is a default implementation which assumes that the message
        space of the encoder is `F^{k}`, where `F` is
        :meth:`sage.coding.linear_code.AbstractLinearCode.base_field`
        and `k` is :meth:`sage.coding.linear_code.AbstractLinearCode.dimension`.
        If this is not the case, this method should be overwritten by the subclass.

        .. NOTE::

            :meth:`encode` might be a partial function over ``self``'s :meth:`message_space`.
            One should use the exception :class:`EncodingError` to catch attempts
            to encode words that are outside of the message space.

        INPUT:

        - ``word`` -- a vector of the message space of the ``self``.

        OUTPUT:

        - a vector of :meth:`code`.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (0, 1, 1, 0))
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.encode(word)
            (1, 1, 0, 0, 1, 1, 0)

        If ``word`` is not in the message space of ``self``, it will return an exception::

            sage: word = random_vector(GF(7), 4)
            sage: E.encode(word)
            Traceback (most recent call last):
            ...
            ValueError: The value to encode must be in Vector space of dimension 4 over Finite Field of size 2
        """
        M = self.message_space()
        if word not in M:
            raise ValueError("The value to encode must be in %s" % M)
        return vector(word) * self.generator_matrix()

    def unencode(self, c, nocheck=False):
        r"""
        Returns the message corresponding to the codeword ``c``.

        This is the inverse of :meth:`encode`.

        INPUT:

        - ``c`` -- a codeword of :meth:`code`.

        - ``nocheck`` -- (default: ``False``) checks if ``c`` is in :meth:`code`. You might set
          this to ``True`` to disable the check for saving computation. Note that if ``c`` is
          not in :meth:`self` and ``nocheck = True``, then the output of :meth:`unencode` is
          not defined (except that it will be in the message space of ``self``).

        OUTPUT:

        - an element of the message space of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: c in C
            True
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.unencode(c)
            (0, 1, 1, 0)

        TESTS:

        If ``nocheck`` is set to ``False``, and one provides a word which is not in
        :meth:`code`, :meth:`unencode` will return an error::

            sage: c = vector(GF(2), (0, 1, 0, 0, 1, 1, 0))
            sage: c in C
            False
            sage: E.unencode(c, False)
            Traceback (most recent call last):
            ...
            EncodingError: Given word is not in the code

        Note that since ticket :trac: `21326`, codes cannot be of length zero::

            sage: G = Matrix(GF(17), [])
            sage: C = LinearCode(G)
            Traceback (most recent call last):
            ...
            ValueError: length must be a non-zero positive integer
       """
        if nocheck == False and c not in self.code():
            raise EncodingError("Given word is not in the code")
        return self.unencode_nocheck(c)

    @cached_method
    def _unencoder_matrix(self):
        r"""
        Finds an information set for the matrix ``G`` returned by :meth:`generator_matrix`,
        and returns the inverse of that submatrix of ``G``.

        AUTHORS:

            This function is taken from codinglib [Nielsen]_

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E._unencoder_matrix()
            (
            [0 0 1 1]
            [0 1 0 1]
            [1 1 1 0]
            [0 1 1 1], (0, 1, 2, 3)
            )
        """
        info_set = self.code().information_set()
        Gt = self.generator_matrix().matrix_from_columns(info_set)
        return (Gt.inverse(), info_set)

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.

        When ``c`` is not a codeword, the output is unspecified.

        AUTHORS:

            This function is taken from codinglib [Nielsen]_

        INPUT:


        - ``c`` -- a codeword of :meth:`code`.

        OUTPUT:

        - an element of the message space of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: c in C
            True
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.unencode_nocheck(c)
            (0, 1, 1, 0)

        Taking a vector that does not belong to ``C`` will not raise an error but
        probably just give a non-sensical result::

            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 1))
            sage: c in C
            False
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.unencode_nocheck(c)
            (0, 1, 1, 0)
            sage: m = vector(GF(2), (0, 1, 1, 0))
            sage: c1 = E.encode(m)
            sage: c == c1
            False
        """
        U, info_set = self._unencoder_matrix()
        cc = vector(self.code().base_ring(), [c[i] for i in info_set])
        return cc * U

    def code(self):
        r"""
        Returns the code for this :class:`Encoder`.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E.code() == C
            True
        """
        return self._code

    def message_space(self):
        r"""
        Returns the ambient space of allowed input to :meth:`encode`.
        Note that :meth:`encode` is possibly a partial function over
        the ambient space.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.code().base_field()**(self.code().dimension())

    @abstract_method(optional = True)
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        This is an abstract method and it should be implemented separately.
        Reimplementing this for each subclass of :class:`Encoder` is not mandatory
        (as a generator matrix only makes sense when the message space is of the `F^k`,
        where `F` is the base field of :meth:`code`.)

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E.generator_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
        """

class EncodingError(Exception):
    r"""
    Special exception class to indicate an error during encoding or unencoding.
    """
    pass
