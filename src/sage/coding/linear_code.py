# -*- coding: utf-8 -*-
r"""
Linear codes

Linear Codes
============

Let `F = \GF{q}` be a finite field. A rank `k` linear subspace of the vector
space `F^n` is called an `[n, k]`-linear code, `n` being the length of the code
and `k` its dimension. Elements of a code `C` are called codewords.

A linear map from `F^k` to an `[n,k]` code `C` is called an "encoding", and it
can be represented as a `k \times n` matrix, called a generator matrix.
Alternatively, `C` can be represented by its orthogonal complement in `F^n`,
i.e. the `(n-k)`-dimensional vector space `C^\perp` such that the inner product
of any element from `C` and any element from `C^\perp` is zero. `C^\perp` is
called the dual code of `C`, and any generator matrix for `C^\perp` is called a
parity check matrix for `C`.

We commonly endow `F^n` with the Hamming metric, i.e. the weight of a vector is
the number of non-zero elements in it. The central operation of a linear code
is then "decoding": given a linear code `C \subset F^n` and a "received word"
`r \in F^n` , retrieve the codeword `c \in C` such that the Hamming distance
between `r` and `c` is minimal.

Families or Generic codes
=========================

Linear codes are either studied as generic vector spaces without any known
structure, or as particular sub-families with special properties.

The class :class:`sage.coding.linear_code.LinearCode` is used to represent the
former.

For the latter, these will be represented by specialised classes; for instance,
the family of Hamming codes are represented by the class
:class:`sage.coding.hamming_code.HammingCode`. Type ``codes.<tab>`` for a list
of all code families known to Sage. Such code family classes should inherit from
the abstract base class :class:`sage.coding.linear_code.AbstractLinearCode`.

``AbstractLinearCode``
----------------------

This is a base class designed to contain methods, features and parameters
shared by every linear code. For instance, generic algorithms for computing the
minimum distance, the covering radius, etc. Many of these algorithms are slow,
e.g. exponential in the code length. For specific subfamilies, better algorithms
or even closed formulas might be known, in which case the respective method
should be overridden.

``AbstractLinearCode`` is an abstract class for linear codes, so any linear code
class should inherit from this class. Also ``AbstractLinearCode`` should never
itself be instantiated.

See :class:`sage.coding.linear_code.AbstractLinearCode` for details and
examples.

``LinearCode``
--------------

This class is used to represent arbitrary and unstructured linear codes.  It
mostly rely directly on generic methods provided by ``AbstractLinearCode``,
which means that basic operations on the code (e.g. computation of the minimum
distance) will use slow algorithms.

A ``LinearCode`` is instantiated by providing a generator matrix::

    sage: M = matrix(GF(2), [[1, 0, 0, 1, 0],\
                             [0, 1, 0, 1, 1],\
                             [0, 0, 1, 1, 1]])
    sage: C = codes.LinearCode(M)
    sage: C
    [5, 3] linear code over GF(2)
    sage: C.generator_matrix()
    [1 0 0 1 0]
    [0 1 0 1 1]
    [0 0 1 1 1]

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C = LinearCode(G)
    sage: C.basis()
    [
    (1, 1, 1, 0, 0, 0, 0),
    (1, 0, 0, 1, 1, 0, 0),
    (0, 1, 0, 1, 0, 1, 0),
    (1, 1, 0, 1, 0, 0, 1)
    ]
    sage: c = C.basis()[1]
    sage: c in C
    True
    sage: c.nonzero_positions()
    [0, 3, 4]
    sage: c.support()
    [0, 3, 4]
    sage: c.parent()
    Vector space of dimension 7 over Finite Field of size 2

Further references
------------------

If you want to get started on Sage's linear codes library, see
https://doc.sagemath.org/html/en/thematic_tutorials/coding_theory.html

If you want to learn more on the design of this library, see
https://doc.sagemath.org/html/en/thematic_tutorials/structures_in_coding_theory.html

REFERENCES:

- [HP2003]_

- [Gu]_

AUTHORS:

- David Joyner (2005-11-22, 2006-12-03): initial version

- William Stein (2006-01-23): Inclusion in Sage

- David Joyner (2006-01-30, 2006-04): small fixes

- David Joyner (2006-07): added documentation, group-theoretical methods,
  ToricCode

- David Joyner (2006-08): hopeful latex fixes to documentation, added list and
  __iter__ methods to LinearCode and examples, added hamming_weight function,
  fixed random method to return a vector, TrivialCode, fixed subtle bug in
  dual_code, added galois_closure method, fixed mysterious bug in
  permutation_automorphism_group (GAP was over-using "G" somehow?)

- David Joyner (2006-08): hopeful latex fixes to documentation, added
  CyclicCode, best_known_linear_code, bounds_minimum_distance,
  assmus_mattson_designs (implementing Assmus-Mattson Theorem).

- David Joyner (2006-09): modified decode syntax, fixed bug in
  is_galois_closed, added LinearCode_from_vectorspace, extended_code,
  zeta_function

- Nick Alexander (2006-12-10): factor GUAVA code to guava.py

- David Joyner (2007-05): added methods punctured, shortened, divisor,
  characteristic_polynomial, binomial_moment, support for
  LinearCode. Completely rewritten zeta_function (old version is now
  zeta_function2) and a new function, LinearCodeFromVectorSpace.

- David Joyner (2007-11): added zeta_polynomial, weight_enumerator,
  chinen_polynomial; improved best_known_code; made some pythonic revisions;
  added is_equivalent (for binary codes)

- David Joyner (2008-01): fixed bug in decode reported by Harald Schilly,
  (with Mike Hansen) added some doctests.

- David Joyner (2008-02): translated standard_form, dual_code to Python.

- David Joyner (2008-03): translated punctured, shortened, extended_code,
  random (and renamed random to random_element), deleted zeta_function2,
  zeta_function3, added wrapper automorphism_group_binary_code to Robert
  Miller's code), added direct_sum_code, is_subcode, is_self_dual,
  is_self_orthogonal, redundancy_matrix, did some alphabetical reorganizing
  to make the file more readable. Fixed a bug in permutation_automorphism_group
  which caused it to crash.

- David Joyner (2008-03): fixed bugs in spectrum and zeta_polynomial, which
  misbehaved over non-prime base rings.

- David Joyner (2008-10): use CJ Tjhal's MinimumWeight if char = 2 or 3 for
  min_dist; add is_permutation_equivalent and improve
  permutation_automorphism_group using an interface with Robert Miller's code;
  added interface with Leon's code for the spectrum method.

- David Joyner (2009-02): added native decoding methods (see module_decoder.py)

- David Joyner (2009-05): removed dependence on Guava, allowing it to be an
  option. Fixed errors in some docstrings.

- Kwankyu Lee (2010-01): added methods generator_matrix_systematic,
  information_set, and magma interface for linear codes.

- Niles Johnson (2010-08): :trac:`3893`: ``random_element()`` should pass on
  ``*args`` and ``**kwds``.

- Thomas Feulner (2012-11): :trac:`13723`: deprecation of ``hamming_weight()``

- Thomas Feulner (2013-10): added methods to compute a canonical representative
  and the automorphism group

TESTS::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C  = LinearCode(G)
    sage: C == loads(dumps(C))
    True
"""
#******************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#******************************************************************************
from __future__ import division, print_function, absolute_import

from six.moves import range
from six import iteritems

import inspect
from copy import copy

from sage.cpython.string import bytes_to_str
from sage.interfaces.all import gap
from sage.categories.modules import Modules
from sage.categories.cartesian_product import cartesian_product
from sage.categories.fields import Fields
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.modules.module import Module
from sage.modules.free_module_element import vector
from sage.arith.all import GCD, binomial
from sage.groups.all import SymmetricGroup
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.structure.parent import Parent
from sage.misc.all import prod
from sage.misc.functional import is_even
from sage.misc.cachefunc import cached_method
from sage.misc.sageinspect import sage_getargspec
from sage.misc.randstate import current_randstate
from sage.combinat.subset import Subsets
from sage.features.gap import GapPackage

from .encoder import Encoder
from .decoder import Decoder

#******************************************************************************
# coding theory functions
#******************************************************************************

def _dump_code_in_leon_format(C):
    r"""
    Writes a file in Sage's temp directory representing the code C, returning
    the absolute path to the file.

    This is the Sage translation of the GuavaToLeon command in Guava's
    codefun.gi file.

    INPUT:

    - ``C`` - a linear code (over GF(p), p < 11)

    OUTPUT:

    - Absolute path to the file written

    EXAMPLES::

        sage: C = codes.HammingCode(GF(2), 3); C
        [7, 4] Hamming Code over GF(2)
        sage: file_loc = sage.coding.linear_code._dump_code_in_leon_format(C)
        sage: f = open(file_loc); print(f.read())
        LIBRARY code;
        code=seq(2,4,7,seq(
        1,0,0,0,0,1,1,
        0,1,0,0,1,0,1,
        0,0,1,0,1,1,0,
        0,0,0,1,1,1,1
        ));
        FINISH;
        sage: f.close()

    """
    from sage.misc.temporary_file import tmp_filename
    F = C.base_ring()
    p = F.order()  # must be prime and <11
    s = "LIBRARY code;\n"+"code=seq(%s,%s,%s,seq(\n"%(p,C.dimension(),C.length())
    Gr = [str(r)[1:-1].replace(" ","") for r in C.generator_matrix().rows()]
    s += ",\n".join(Gr) + "\n));\nFINISH;"
    file_loc = tmp_filename()
    f = open(file_loc,"w")
    f.write(s)
    f.close()

    return file_loc

def _explain_constructor(cl):
    r"""
    Internal function for use error messages when constructing encoders and decoders.

    EXAMPLES::

        sage: from sage.coding.linear_code import _explain_constructor, LinearCodeSyndromeDecoder
        sage: cl = LinearCodeSyndromeDecoder
        sage: _explain_constructor(cl)
        "The constructor requires no arguments.\nIt takes the optional
        arguments ['maximum_error_weight'].\nSee the documentation of
        sage.coding.linear_code.LinearCodeSyndromeDecoder for more details."

        sage: from sage.coding.information_set_decoder import LinearCodeInformationSetDecoder
        sage: cl = LinearCodeInformationSetDecoder
        sage: _explain_constructor(cl)
        "The constructor requires the arguments ['number_errors'].\nIt takes the optional arguments ['algorithm'].\nIt accepts unspecified arguments as well.\nSee the documentation of sage.coding.information_set_decoder.LinearCodeInformationSetDecoder for more details."
    """
    if inspect.isclass(cl):
        argspec = sage_getargspec(cl.__init__)
        skip = 2 # skip the self and code arguments
    else:
        # Not a class, assume it's a factory function posing as a class
        argspec = sage_getargspec(cl)
        skip = 1 # skip code argument
    if argspec.defaults:
        args = argspec.args[skip:-len(argspec.defaults)]
        kwargs = argspec.args[-len(argspec.defaults):]
        opts = "It takes the optional arguments {}.".format(kwargs)
    else:
        args = argspec.args[skip:]
        opts = "It takes no optional arguments."
    if args:
        reqs = "The constructor requires the arguments {}.".format(args)
    else:
        reqs = "The constructor requires no arguments."
    if argspec.varargs or argspec.keywords:
        var = "It accepts unspecified arguments as well.\n"
    else:
        var = ""
    return("{}\n{}\n{}See the documentation of {}.{} for more details."\
            .format(reqs, opts, var, cl.__module__, cl.__name__))


class AbstractLinearCode(Module):
    """
    Abstract base class for linear codes.

    This class contains all methods that can be used on Linear Codes and on
    Linear Codes families.  So, every Linear Code-related class should inherit
    from this abstract class.

    To implement a linear code, you need to:

    - inherit from AbstractLinearCode

    - call AbstractLinearCode ``__init__`` method in the subclass constructor. Example:
      ``super(SubclassName, self).__init__(base_field, length, "EncoderName", "DecoderName")``.
      By doing that, your subclass will have its ``length`` parameter
      initialized and will be properly set as a member of the category framework.
      You need of course to complete the constructor by adding any additional parameter
      needed to describe properly the code defined in the subclass.

    - Add the following two lines on the class level::

          _registered_encoders = {}
          _registered_decoders = {}


    - fill the dictionary of its encoders in ``sage.coding.__init__.py`` file. Example:
      I want to link the encoder ``MyEncoderClass`` to ``MyNewCodeClass``
      under the name ``MyEncoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_encoders["NameOfMyEncoder"] = MyEncoderClass`` and all instances of
      ``MyNewCodeClass`` will be able to use instances of ``MyEncoderClass``.

    - fill the dictionary of its decoders in ``sage.coding.__init__`` file. Example:
      I want to link the encoder ``MyDecoderClass`` to ``MyNewCodeClass``
      under the name ``MyDecoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_decoders["NameOfMyDecoder"] = MyDecoderClass`` and all instances of
      ``MyNewCodeClass`` will be able to use instances of ``MyDecoderClass``.


    As AbstractLinearCode is not designed to be implemented, it does not have any representation
    methods. You should implement ``_repr_`` and ``_latex_`` methods in the subclass.

    .. NOTE::

        :class:`AbstractLinearCode` has a generic implementation of the
        method ``__eq__`` which uses the generator matrix and is quite
        slow. In subclasses you are encouraged to override ``__eq__``
        and ``__hash__``.

    .. WARNING::

        The default encoder should always have `F^{k}` as message space, with `k` the dimension
        of the code and `F` is the base ring of the code.

        A lot of methods of the abstract class rely on the knowledge of a generator matrix.
        It is thus strongly recommended to set an encoder with a generator matrix implemented
        as a default encoder.

        TESTS::

        This class uses the following experimental feature:
        :class:`sage.coding.relative_finite_field_extension.RelativeFiniteFieldExtension`.
        This test block is here only to trigger the experimental warning so it does not
        interferes with doctests::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: RelativeFiniteFieldExtension(Fqm, Fq)
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See http://trac.sagemath.org/20284 for details.
            Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, default_encoder_name, default_decoder_name):
        """
        Initializes mandatory parameters that any linear code shares.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear code. The class
        :class:`sage.coding.linear_code.AbstractLinearCode` should never be
        directly instantiated.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer, must be > 0)

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``

        EXAMPLES:

        The following example demonstrates how to subclass `AbstractLinearCode`
        for representing a new family of codes. The example family is non-sensical::

            sage: class MyCodeFamily(sage.coding.linear_code.AbstractLinearCode):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "GeneratorMatrix", "Syndrome")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

        We now instantiate a member of our newly made code family::

            sage: generator_matrix = matrix(GF(17), 5, 10,
            ....:                           {(i,i):1 for i in range(5)})
            sage: C = MyCodeFamily(GF(17), 10, 5, generator_matrix)

        We can check its existence and parameters::

            sage: C
            [10, 5] dummy code over GF(17)

        We can check that it is truly a part of the framework category::

            sage: C.parent()
            <class '__main__.MyCodeFamily_with_category'>
            sage: C.category()
            Category of facade finite dimensional vector spaces with basis over Finite Field of size 17

        And any method that works on linear codes works for our new dummy code::

            sage: C.minimum_distance()
            1
            sage: C.is_self_orthogonal()
            False
            sage: print(C.divisor()) #long time
            1

        TESTS:

        If the length field is neither a Python int nor a Sage Integer, it will
        raise a exception::

            sage: C = MyCodeFamily(GF(17), 10.0, 5, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: length must be a positive integer

        If the length of the code is not a non-zero positive integer
        (See :trac:`21326`), it will raise an exception::

            sage: empty_generator_matrix = Matrix(GF(17),0,1)
            sage: C = MyCodeFamily(GF(17), 0, 1, empty_generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: length must be a positive integer

        If the name of the default decoder is not known by the class, it will raise
        a exception::

            sage: class MyCodeFamily2(sage.coding.linear_code.AbstractLinearCode):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "GeneratorMatrix", "Fail")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

            sage: C = MyCodeFamily2(GF(17), 10, 5, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders

        If the name of the default encoder is not known by the class, it will raise
        an exception::

            sage: class MyCodeFamily3(sage.coding.linear_code.AbstractLinearCode):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "Fail", "Syndrome")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

            sage: C = MyCodeFamily3(GF(17), 10, 5, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders

        A ring instead of a field::

            sage: codes.LinearCode(IntegerModRing(4),matrix.ones(4))
            Traceback (most recent call last):
            ...
            ValueError: 'generator' must be defined on a field (not a ring)
        """
        from sage.coding.information_set_decoder import LinearCodeInformationSetDecoder

        # Add here any generic encoder or decoder. This allows any class which
        # inherits from AbstractLinearCode to use generic decoders/encoders
        self._registered_encoders['Systematic'] = LinearCodeSystematicEncoder
        self._registered_decoders['Syndrome'] = LinearCodeSyndromeDecoder
        self._registered_decoders['NearestNeighbor'] = LinearCodeNearestNeighborDecoder
        self._registered_decoders['InformationSet'] = LinearCodeInformationSetDecoder

        if not isinstance(length, (int, Integer)) or length <= 0:
            raise ValueError("length must be a positive integer")
        if not base_field.is_field():
            raise ValueError("{} is not a field".format(base_field))
        if not default_encoder_name in self._registered_encoders:
            raise ValueError("You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders")
        if not default_decoder_name in self._registered_decoders:
            raise ValueError("You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders")

        self._length = Integer(length)
        self._default_decoder_name = default_decoder_name
        self._default_encoder_name = default_encoder_name
        cat = Modules(base_field).FiniteDimensional().WithBasis().Finite()
        facade_for = VectorSpace(base_field, self._length)
        self.Element = type(facade_for.an_element()) #for when we made this a non-facade parent
        Parent.__init__(self, base=base_field, facade=facade_for, category=cat)

    def __getstate__(self):
        """
        Used for pickling codes.

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: '_registered_encoders' in C.__getstate__()
            True
        """
        d = super(AbstractLinearCode, self).__getstate__()
        d['_registered_encoders'] = self._registered_encoders
        d['_registered_decoders'] = self._registered_decoders
        return d

    def _repr_(self):
        r"""
        Return an error message requiring to override ``_repr_`` in ``self``.

        As one has to implement specific representation methods (`_repr_` and `_latex_`)
        when writing a new code class which inherits from :class:`AbstractLinearCode`,
        the generic call to `_repr_` has to fail.

        EXAMPLES:

        This was taken from :trac:`20899` (and thus ensures this method fixes what was
        described in this ticket).

        We create a new code class, its dedicated encoder
        and set appropriate parameters::

            sage: from sage.coding.linear_code import AbstractLinearCode
            sage: from sage.coding.encoder import Encoder
            sage: class MyCode(AbstractLinearCode):
            ....:    _registered_encoders = {}
            ....:    _registered_decoders = {}
            ....:    def __init__(self):
            ....:        super(MyCode, self).__init__(GF(5), 10, "Monkey", "Syndrome")
            ....:        self._dimension = 2

            sage: class MonkeyEncoder(Encoder):
            ....:    def __init__(self, C):
            ....:        super(MonkeyEncoder, self).__init__(C)
            ....:    @cached_method
            ....:    def generator_matrix(self):
            ....:        return matrix(GF(5), 2, 10, [ [1]*5 + [0]*5, [0]*5 + [1]*5 ])
            sage: MyCode._registered_encoders["Monkey"] = MonkeyEncoder
            sage: MyCode._registered_decoders["Syndrome"] = codes.decoders.LinearCodeSyndromeDecoder

        We check we get a sensible error message while asking for a string
        representation of an instance of our new class:

            sage: C = MyCode()
            sage: C #random
            Traceback (most recent call last):
            ...
            RuntimeError: Please override _repr_ in the implementation of <class '__main__.MyCode_with_category'>
        """
        raise RuntimeError("Please override _repr_ in the implementation of {}".format(self.parent()))

    def _latex_(self):
        r"""
        Return an error message requiring to override ``_latex_`` in ``self``.

        As one has to implement specific representation methods (`_repr_` and `_latex_`)
        when writing a new code class which inherits from :class:`AbstractLinearCode`,
        the generic call to `_latex_` has to fail.

        EXAMPLES:

        This was taken from :trac:`20899` (and thus ensures this method fixes what was
        described in this ticket).

        We create a new code class, its dedicated encoder
        and set appropriate parameters::

            sage: from sage.coding.linear_code import AbstractLinearCode
            sage: from sage.coding.encoder import Encoder
            sage: class MyCode(AbstractLinearCode):
            ....:    _registered_encoders = {}
            ....:    _registered_decoders = {}
            ....:    def __init__(self):
            ....:        super(MyCode, self).__init__(GF(5), 10, "Monkey", "Syndrome")
            ....:        self._dimension = 2

            sage: class MonkeyEncoder(Encoder):
            ....:    def __init__(self, C):
            ....:        super(MonkeyEncoder, self).__init__(C)
            ....:    @cached_method
            ....:    def generator_matrix(self):
            ....:        return matrix(GF(5), 2, 10, [ [1]*5 + [0]*5, [0]*5 + [1]*5 ])
            sage: MyCode._registered_encoders["Monkey"] = MonkeyEncoder
            sage: MyCode._registered_decoders["Syndrome"] = codes.decoders.LinearCodeSyndromeDecoder

        We check we get a sensible error message while asking for a string
        representation of an instance of our new class:

            sage: C = MyCode()
            sage: latex(C)
            Traceback (most recent call last):
            ...
            RuntimeError: Please override _latex_ in the implementation of <class '__main__.MyCode_with_category'>
        """
        raise RuntimeError("Please override _latex_ in the implementation of {}".format(self.parent()))

    def _an_element_(self):
        r"""
        Return an element of the linear code. Currently, it simply returns
        the first row of the generator matrix.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.an_element()
            (1, 0, 0, 0, 0, 1, 1)
            sage: C2 = C.cartesian_product(C)
            sage: C2.an_element()
            ((1, 0, 0, 0, 0, 1, 1), (1, 0, 0, 0, 0, 1, 1))
        """
        return self.gens()[0]

    def add_decoder(self, name, decoder):
        r"""
        Adds an decoder to the list of registered decoders of ``self``.

        .. NOTE::

            This method only adds ``decoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.decoder.Decoder`, please refer
            to the documentation of :class:`AbstractLinearCode`.

        INPUT:

        - ``name`` -- the string name for the decoder

        - ``decoder`` -- the class name of the decoder

        EXAMPLES:

        First of all, we create a (very basic) new decoder::

            sage: class MyDecoder(sage.coding.decoder.Decoder):
            ....:   def __init__(self, code):
            ....:       super(MyDecoder, self).__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyDecoder decoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new decoder to the list of available decoders of C::

            sage: C.add_decoder("MyDecoder", MyDecoder)
            sage: sorted(C.decoders_available())
            ['InformationSet', 'MyDecoder', 'NearestNeighbor', 'Syndrome']

        We can verify that any new code will not know MyDecoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: sorted(C2.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']

        TESTS:

        It is impossible to use a name which is in the dictionary of available decoders::

            sage: C.add_decoder("Syndrome", MyDecoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered decoder with this name
        """
        if self._registered_decoders == self.__class__._registered_decoders:
            self._registered_decoders = copy(self._registered_decoders)
            reg_dec = self._registered_decoders
            if name in reg_dec:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder
        else:
            if name in self._registered_decoders:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder

    def add_encoder(self, name, encoder):
        r"""
        Adds an encoder to the list of registered encoders of ``self``.

        .. NOTE::

            This method only adds ``encoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.encoder.Encoder`, please refer
            to the documentation of :class:`AbstractLinearCode`.

        INPUT:

        - ``name`` -- the string name for the encoder

        - ``encoder`` -- the class name of the encoder

        EXAMPLES:

        First of all, we create a (very basic) new encoder::

            sage: class MyEncoder(sage.coding.encoder.Encoder):
            ....:   def __init__(self, code):
            ....:       super(MyEncoder, self).__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyEncoder encoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new encoder to the list of available encoders of C::

            sage: C.add_encoder("MyEncoder", MyEncoder)
            sage: sorted(C.encoders_available())
            ['MyEncoder', 'Systematic']

        We can verify that any new code will not know MyEncoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: sorted(C2.encoders_available())
            ['Systematic']

        TESTS:

        It is impossible to use a name which is in the dictionary of available encoders::

            sage: C.add_encoder("Systematic", MyEncoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered encoder with this name
        """
        if self._registered_encoders == self.__class__._registered_encoders:
            self._registered_encoders = copy(self._registered_encoders)
            reg_enc = self._registered_encoders
            if name in reg_enc:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder
        else:
            if name in self._registered_encoders:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder

    def automorphism_group_gens(self, equivalence="semilinear"):
        r"""
        Return generators of the automorphism group of ``self``.

        INPUT:

        - ``equivalence`` (optional) -- which defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        OUTPUT:

        - generators of the automorphism group of ``self``
        - the order of the automorphism group of ``self``

        EXAMPLES:

        Note, this result can depend on the PRNG state in libgap in a way that
        depends on which packages are loaded, so we must re-seed GAP to ensure
        a consistent result for this example::

            sage: libgap.set_seed(0)
            0
            sage: C = codes.HammingCode(GF(4, 'z'), 3)
            sage: C.automorphism_group_gens()
            ([((1, z, z + 1, z, z, 1, 1, z, z + 1, z, z, 1, z, z + 1, z, z, 1, z, z + 1, z, z); (1,5,18,7,11,8)(2,12,21)(3,20,14,10,19,15)(4,9)(13,17,16), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z + 1),
              ((1, 1, z, z + 1, z, z, z + 1, z + 1, z, 1, 1, z, z, z + 1, z + 1, 1, z, z, 1, z, z + 1); (2,11)(3,13)(4,14)(5,20)(6,17)(8,15)(16,19), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z + 1),
              ((z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z); (), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z)],
             362880)
            sage: C.automorphism_group_gens(equivalence="linear")
            ([((z, 1, 1, z, z + 1, z, z, z + 1, z + 1, z + 1, 1, z + 1, z, z, 1, 1, 1, z, z, z + 1, z); (1,6)(2,20,9,16)(3,10,8,11)(4,15,21,5)(12,17)(13,14,19,18), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z),
              ((1, z, z + 1, z, z, z, z + 1, z + 1, 1, z, z, z, 1, z, 1, z + 1, z, z + 1, z, z + 1, 1); (1,15,20,5,8,6,12,14,13,7,16,11,19,3,21,4,9,10,18,17,2), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z),
              ((z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1); (), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z)],
             181440)
            sage: C.automorphism_group_gens(equivalence="permutational")
            ([((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (1,19)(3,17)(4,21)(5,20)(7,14)(9,12)(10,16)(11,15), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z),
              ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (1,11)(3,10)(4,9)(5,7)(12,21)(14,20)(15,19)(16,17), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z),
              ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (1,17)(2,8)(3,14)(4,10)(7,12)(9,19)(13,18)(15,20), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z),
              ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (2,13)(3,14)(4,20)(5,11)(8,18)(9,19)(10,15)(16,21), Ring endomorphism of Finite Field in z of size 2^2
                Defn: z |--> z)],
            64)
        """
        aut_group_can_label = self._canonize(equivalence)
        return aut_group_can_label.get_autom_gens(), \
               aut_group_can_label.get_autom_order()

    def ambient_space(self):
        r"""
        Returns the ambient vector space of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.ambient_space()
            Vector space of dimension 7 over Finite Field of size 2
        """
        return VectorSpace(self.base_ring(),self.length())

    def assmus_mattson_designs(self, t, mode=None):
        r"""
        Assmus and Mattson Theorem (section 8.4, page 303 of [HP2003]_): Let
        `A_0, A_1, ..., A_n` be the weights of the codewords in a binary
        linear `[n , k, d]` code `C`, and let `A_0^*, A_1^*, ..., A_n^*` be
        the weights of the codewords in its dual `[n, n-k, d^*]` code `C^*`.
        Fix a `t`, `0<t<d`, and let

        .. MATH::

           s = |\{ i\ |\ A_i^* \not= 0, 0< i \leq n-t\}|.

        Assume `s\leq d-t`.

        1. If `A_i\not= 0` and `d\leq i\leq n`
           then `C_i = \{ c \in C\ |\ wt(c) = i\}` holds a simple t-design.

        2. If `A_i^*\not= 0` and `d*\leq i\leq n-t` then
           `C_i^* = \{ c \in C^*\ |\ wt(c) = i\}` holds a simple t-design.

        A block design is a pair `(X,B)`, where `X` is a non-empty finite set
        of `v>0` elements called points, and `B` is a non-empty finite
        multiset of size b whose elements are called blocks, such that each
        block is a non-empty finite multiset of `k` points. `A` design without
        repeated blocks is called a simple block design. If every subset of
        points of size `t` is contained in exactly `\lambda` blocks the block
        design is called a `t-(v,k,\lambda)` design (or simply a `t`-design
        when the parameters are not specified). When `\lambda=1` then the
        block design is called a `S(t,k,v)` Steiner system.

        In the Assmus and Mattson Theorem (1), `X` is the set `\{1,2,...,n\}`
        of coordinate locations and `B = \{supp(c)\ |\ c \in C_i\}` is the set
        of supports of the codewords of `C` of weight `i`. Therefore, the
        parameters of the `t`-design for `C_i` are

        ::

            t =       given
            v =       n
            k =       i   (k not to be confused with dim(C))
            b =       Ai
            lambda = b*binomial(k,t)/binomial(v,t) (by Theorem 8.1.6,
                                                       p 294, in [HP2003]_)

        Setting the ``mode="verbose"`` option prints out the values of the
        parameters.

        The first example below means that the binary [24,12,8]-code C has
        the property that the (support of the) codewords of weight 8 (resp.,
        12, 16) form a 5-design. Similarly for its dual code `C^*` (of course
        `C=C^*` in this case, so this info is extraneous). The test fails to
        produce 6-designs (ie, the hypotheses of the theorem fail to hold,
        not that the 6-designs definitely don't exist). The command
        assmus_mattson_designs(C,5,mode="verbose") returns the same value
        but prints out more detailed information.

        The second example below illustrates the blocks of the 5-(24, 8, 1)
        design (i.e., the S(5,8,24) Steiner system).

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))             #  example 1
            sage: C.assmus_mattson_designs(5)
            ['weights from C: ',
            [8, 12, 16, 24],
            'designs from C: ',
            [[5, (24, 8, 1)], [5, (24, 12, 48)], [5, (24, 16, 78)], [5, (24, 24, 1)]],
            'weights from C*: ',
            [8, 12, 16],
            'designs from C*: ',
            [[5, (24, 8, 1)], [5, (24, 12, 48)], [5, (24, 16, 78)]]]
            sage: C.assmus_mattson_designs(6)
            0
            sage: X = range(24)                           #  example 2
            sage: blocks = [c.support() for c in C if c.hamming_weight()==8]; len(blocks)  # long time computation
            759
        """
        C = self
        ans = []
        G = C.generator_matrix()
        n = len(G.columns())
        Cp = C.dual_code()
        wts = C.weight_distribution()
        d = min([i for i in range(1,len(wts)) if wts[i]!=0])
        if t>=d:
            return 0
        nonzerowts = [i for i in range(len(wts)) if wts[i]!=0 and i<=n and i>=d]
        if mode=="verbose":
            for w in nonzerowts:
                print("The weight w={} codewords of C* form a t-(v,k,lambda) design, where\n \
                        t={}, v={}, k={}, lambda={}. \nThere are {} block of this design.".format(\
                        w,t,n,w,wts[w]*binomial(w,t)//binomial(n,t),wts[w]))
        wtsp = Cp.weight_distribution()
        dp = min([i for i in range(1,len(wtsp)) if wtsp[i]!=0])
        nonzerowtsp = [i for i in range(len(wtsp)) if wtsp[i]!=0 and i<=n-t and i>=dp]
        s = len([i for i in range(1,n) if wtsp[i]!=0 and i<=n-t and i>0])
        if mode=="verbose":
            for w in nonzerowtsp:
                print("The weight w={} codewords of C* form a t-(v,k,lambda) design, where\n \
                        t={}, v={}, k={}, lambda={}. \nThere are {} block of this design.".format(\
                        w,t,n,w,wts[w]*binomial(w,t)//binomial(n,t),wts[w]))
        if s<=d-t:
            des = [[t,(n,w,wts[w]*binomial(w,t)//binomial(n,t))] for w in nonzerowts]
            ans = ans + ["weights from C: ",nonzerowts,"designs from C: ",des]
            desp = [[t,(n,w,wtsp[w]*binomial(w,t)//binomial(n,t))] for w in nonzerowtsp]
            ans = ans + ["weights from C*: ",nonzerowtsp,"designs from C*: ",desp]
            return ans
        return 0

    def base_field(self):
        r"""
        Return the base field of ``self``.

        EXAMPLES::

            sage: G  = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C.base_field()
            Finite Field of size 2
        """
        return self.base_ring()

    def basis(self):
        r"""
        Returns a basis of ``self``.

        OUTPUT:

        -  ``Sequence`` - an immutable sequence whose universe is ambient space of ``self``.

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
        return Sequence(gens, universe=self.ambient_space(), check = False, immutable=True, cr=True)

    # S. Pancratz, 19 Jan 2010:  In the doctests below, I removed the example
    # ``C.binomial_moment(3)``, which was also marked as ``#long``.  This way,
    # we shorten the doctests time while still maintaining a zero and a
    # non-zero example.
    def binomial_moment(self, i):
        r"""
        Returns the i-th binomial moment of the `[n,k,d]_q`-code `C`:

        .. MATH::

            B_i(C) = \sum_{S, |S|=i} \frac{q^{k_S}-1}{q-1}

        where `k_S` is the dimension of the shortened code `C_{J-S}`,
        `J=[1,2,...,n]`. (The normalized binomial moment is
        `b_i(C) = \binom(n,d+i)^{-1}B_{d+i}(C)`.) In other words, `C_{J-S}`
        is isomorphic to the subcode of C of codewords supported on S.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.binomial_moment(2)
            0
            sage: C.binomial_moment(4)    # long time
            35

        .. warning::

            This is slow.

        REFERENCE:

        - [Du2004]_
        """
        n = self.length()
        k = self.dimension()
        d = self.minimum_distance()
        F = self.base_ring()
        q = F.order()
        J = range(1,n+1)
        Cp = self.dual_code()
        dp = Cp.minimum_distance()
        if i<d:
            return 0
        if i>n-dp and i<=n:
            return binomial(n,i)*(q**(i+k-n) -1)//(q-1)
        from sage.combinat.set_partition import SetPartitions
        P = SetPartitions(J,2).list()
        b = QQ(0)
        for p in P:
            p = list(p)
            S = p[0]
            if len(S)==n-i:
                C_S = self.shortened(S)
                k_S = C_S.dimension()
                b = b + (q**(k_S) -1)//(q-1)
        return b

    @cached_method
    def _canonize(self, equivalence):
        r"""
        Compute a canonical representative and the automorphism group
        under the action of the semimonomial transformation group.

        INPUT:

        - ``equivalence`` -- which defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        EXAMPLES:

        Note, this result can depend on the PRNG state in libgap in a way that
        depends on which packages are loaded, so we must re-seed GAP to ensure
        a consistent result for this example::

            sage: libgap.set_seed(0)
            0
            sage: C = codes.HammingCode(GF(4, 'z'), 3)
            sage: aut_group_can_label = C._canonize("semilinear")
            sage: C_iso = LinearCode(aut_group_can_label.get_transporter()*C.generator_matrix())
            sage: C_iso == aut_group_can_label.get_canonical_form()
            True
            sage: aut_group_can_label.get_autom_gens()
            [((1, z, z + 1, z, z, 1, 1, z, z + 1, z, z, 1, z, z + 1, z, z, 1, z, z + 1, z, z); (1,5,18,7,11,8)(2,12,21)(3,20,14,10,19,15)(4,9)(13,17,16), Ring endomorphism of Finite Field in z of size 2^2
               Defn: z |--> z + 1),
             ((1, 1, z, z + 1, z, z, z + 1, z + 1, z, 1, 1, z, z, z + 1, z + 1, 1, z, z, 1, z, z + 1); (2,11)(3,13)(4,14)(5,20)(6,17)(8,15)(16,19), Ring endomorphism of Finite Field in z of size 2^2
               Defn: z |--> z + 1),
             ((z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z); (), Ring endomorphism of Finite Field in z of size 2^2
               Defn: z |--> z)]
        """
        from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
        return LinearCodeAutGroupCanLabel(self, algorithm_type=equivalence)

    def canonical_representative(self, equivalence="semilinear"):
        r"""
        Compute a canonical orbit representative under the action of the
        semimonomial transformation group.

        See :mod:`sage.coding.codecan.autgroup_can_label`
        for more details, for example if you would like to compute
        a canonical form under some more restrictive notion of equivalence,
        i.e. if you would like to restrict the permutation group
        to a Young subgroup.

        INPUT:

        - ``equivalence`` (optional) -- which defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        OUTPUT:

        - a canonical representative of ``self``
        - a semimonomial transformation mapping ``self`` onto its representative

        EXAMPLES::

            sage: F.<z> = GF(4)
            sage: C = codes.HammingCode(F, 3)
            sage: CanRep, transp = C.canonical_representative()

        Check that the transporter element is correct::

            sage: LinearCode(transp*C.generator_matrix()) == CanRep
            True

        Check if an equivalent code has the same canonical representative::

            sage: f = F.hom([z**2])
            sage: C_iso = LinearCode(C.generator_matrix().apply_map(f))
            sage: CanRep_iso, _ = C_iso.canonical_representative()
            sage: CanRep_iso == CanRep
            True

        Since applying the Frobenius automorphism could be extended to an
        automorphism of `C`, the following must also yield ``True``::

            sage: CanRep1, _ = C.canonical_representative("linear")
            sage: CanRep2, _ = C_iso.canonical_representative("linear")
            sage: CanRep2 == CanRep1
            True

        TESTS:

        Check that interrupting this does not segfault
        (see :trac:`21651`)::

            sage: C = LinearCode(random_matrix(GF(47), 25, 35))
            sage: alarm(0.5); C.canonical_representative()
            Traceback (most recent call last):
            ...
            AlarmInterrupt
        """
        aut_group_can_label = self._canonize(equivalence)
        return aut_group_can_label.get_canonical_form(), \
               aut_group_can_label.get_transporter()

    def __contains__(self, v):
        r"""
        Returns True if `v` can be coerced into ``self``. Otherwise, returns False.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: vector((1, 0, 0, 0, 0, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 2, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 0, 1/2, 1)) in C # indirect doctest
            False
        """
        if not v in self.ambient_space() or len(v) != self.length():
            return False
        return self.syndrome(v) == 0

    def characteristic(self):
        r"""
        Returns the characteristic of the base ring of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.characteristic()
            2
        """
        return (self.base_ring()).characteristic()

    def characteristic_polynomial(self):
        r"""
        Returns the characteristic polynomial of a linear code, as defined in
        [Lin1999]_.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: C.characteristic_polynomial()
            -4/3*x^3 + 64*x^2 - 2816/3*x + 4096
        """
        R = PolynomialRing(QQ,"x")
        x = R.gen()
        C = self
        Cd = C.dual_code()
        Sd = Cd.support()
        k = C.dimension()
        n = C.length()
        q = (C.base_ring()).order()
        return q**(n-k)*prod([1-x/j for j in Sd if j>0])

    def chinen_polynomial(self):
        """
        Returns the Chinen zeta polynomial of the code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.chinen_polynomial()       # long time
            1/5*(2*sqrt(2)*t^3 + 2*sqrt(2)*t^2 + 2*t^2 + sqrt(2)*t + 2*t + 1)/(sqrt(2) + 1)
            sage: C = codes.GolayCode(GF(3), False)
            sage: C.chinen_polynomial()       # long time
            1/7*(3*sqrt(3)*t^3 + 3*sqrt(3)*t^2 + 3*t^2 + sqrt(3)*t + 3*t + 1)/(sqrt(3) + 1)

        This last output agrees with the corresponding example given in
        Chinen's paper below.

        REFERENCES:

        - Chinen, K. "An abundance of invariant polynomials satisfying the
          Riemann hypothesis", April 2007 preprint.
        """
        from sage.functions.all import sqrt
        C = self
        n = C.length()
        RT = PolynomialRing(QQ,2,"Ts")
        T,s = RT.fraction_field().gens()
        t = PolynomialRing(QQ,"t").gen()
        Cd = C.dual_code()
        k = C.dimension()
        q = (C.base_ring()).characteristic()
        d = C.minimum_distance()
        dperp = Cd.minimum_distance()
        if dperp > d:
            P = RT(C.zeta_polynomial())
            # Sage does not find dealing with sqrt(int) *as an algebraic object*
            # an easy thing to do. Some tricky gymnastics are used to
            # make Sage deal with objects over QQ(sqrt(q)) nicely.
            if is_even(n):
                Pd = q**(k-n//2) * RT(Cd.zeta_polynomial()) * T**(dperp - d)
            else:
                Pd = s * q**(k-(n+1)//2) * RT(Cd.zeta_polynomial()) * T**(dperp - d)
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))
        if dperp < d:
            P = RT(C.zeta_polynomial())*T**(d - dperp)
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))
        if dperp == d:
            P = RT(C.zeta_polynomial())
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

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

    @cached_method
    def covering_radius(self):
        r"""
        Return the minimal integer `r` such that any element in the ambient space of ``self`` has distance at most `r` to a codeword of ``self``.

        This method requires the optional GAP package Guava.

        If the covering radius a code equals its minimum distance, then the code is called perfect.

        .. NOTE::

            This method is currently not implemented on codes over base fields
            of cardinality greater than 256 due to limitations in the underlying
            algorithm of GAP.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 5)
            sage: C.covering_radius()  # optional - gap_packages (Guava package)
            1

            sage: C = codes.random_linear_code(GF(263), 5, 1)
            sage: C.covering_radius()  # optional - gap_packages (Guava package)
            Traceback (most recent call last):
            ...
            NotImplementedError: the GAP algorithm that Sage is using is limited to computing with fields of size at most 256
        """
        GapPackage("guava", spkg="gap_packages").require()
        gap.load_package("guava")
        F = self.base_ring()
        if F.cardinality() > 256:
            raise NotImplementedError("the GAP algorithm that Sage is using "
                                      "is limited to computing with fields "
                                      "of size at most 256")
        gapG = gap(self.generator_matrix())
        C = gapG.GeneratorMatCode(gap(F))
        r = C.CoveringRadius()
        try:
            return ZZ(r)
        except TypeError:
            raise RuntimeError("the covering radius of this code cannot be computed by Guava")

    def decode_to_code(self, word, decoder_name=None, *args, **kwargs):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``word`` -- a vector of the same length as ``self`` over
          the base field of ``self``

        - ``decoder_name`` -- (default: ``None``) Name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT:

        - A vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: C.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decode_to_code(w_err, 'NearestNeighbor')
            (1, 1, 0, 0, 1, 1, 0)
        """
        D = self.decoder(decoder_name, *args, **kwargs)
        return D.decode_to_code(word)

    def decode_to_message(self, word, decoder_name=None, *args, **kwargs):
        r"""
        Correct the errors in word and decodes it to the message space.

        INPUT:

        - ``word`` -- a vector of the same length as ``self`` over the
          base field of ``self``

        - ``decoder_name`` -- (default: ``None``) Name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT:

        - A vector of the message space of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.decode_to_message(word)
            (0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decode_to_message(word, 'NearestNeighbor')
            (0, 1, 1, 0)
        """
        return self.unencode(self.decode_to_code(word, decoder_name, *args, **kwargs), **kwargs)

    @cached_method
    def decoder(self, decoder_name=None, *args, **kwargs):
        r"""
        Return a decoder of ``self``.

        INPUT:

        - ``decoder_name`` -- (default: ``None``) name of the decoder which will be
          returned. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments will be forwarded to the constructor of the decoder
          that will be returned by this method

        OUTPUT:

        - a decoder object

        Besides creating the decoder and returning it, this method also stores
        the decoder in a cache. With this behaviour, each decoder will be created
        at most one time for ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoder()
            Syndrome decoder for [7, 4] linear code over GF(2) handling errors of weight up to 1


        If the name of a decoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decoder('Try')
            Traceback (most recent call last):
            ...
            ValueError: There is no Decoder named 'Try'. The known Decoders are: ['InformationSet', 'NearestNeighbor', 'Syndrome']

        Some decoders take extra arguments. If the user forgets to supply these,
        the error message attempts to be helpful::

            sage: C.decoder('InformationSet')
            Traceback (most recent call last):
            ...
            ValueError: Constructing the InformationSet decoder failed, possibly due to missing or incorrect parameters.
            The constructor requires the arguments ['number_errors'].
            It takes the optional arguments ['algorithm'].
            It accepts unspecified arguments as well.
            See the documentation of sage.coding.information_set_decoder.LinearCodeInformationSetDecoder for more details.

        """
        if decoder_name is None:
            decoder_name = self._default_decoder_name
        if decoder_name in self._registered_decoders:
            decClass = self._registered_decoders[decoder_name]
            try:
                return decClass(self, *args, **kwargs)
            except TypeError:
                raise ValueError(
                        "Constructing the {0} decoder failed, possibly due "
                        "to missing or incorrect parameters.\n{1}".format(
                            decoder_name, _explain_constructor(decClass)))
        else:
            raise ValueError(
                    "There is no Decoder named '{0}'. The known Decoders are: "
                    "{1}".format(decoder_name, self.decoders_available()))

    def decoders_available(self, classes=False):
        r"""
        Returns a list of the available decoders' names for ``self``.

        INPUT:

        - ``classes`` -- (default: ``False``) if ``classes`` is set to ``True``,
          return instead a ``dict`` mapping available decoder name to the
          associated decoder class.

        OUTPUT: a list of strings, or a `dict` mapping strings to classes.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoders_available()
            ['InformationSet', 'NearestNeighbor', 'Syndrome']

            sage: dictionary = C.decoders_available(True)
            sage: sorted(dictionary.keys())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: dictionary['NearestNeighbor']
            <class 'sage.coding.linear_code.LinearCodeNearestNeighborDecoder'>
        """
        if classes:
            return copy(self._registered_decoders)

        return sorted(self._registered_decoders)

    def divisor(self):
        r"""
        Returns the greatest common divisor of the weights of the nonzero codewords.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: C.divisor()   # Type II self-dual
            4
            sage: C = codes.QuadraticResidueCodeEvenPair(17,GF(2))[0]
            sage: C.divisor()
            2
        """
        C = self
        A = C.weight_distribution()
        n = C.length()
        V = VectorSpace(QQ,n+1)
        S = V(A).nonzero_positions()
        S0 = [S[i] for i in range(1,len(S))]
        if len(S)>1: return GCD(S0)
        return 1

    def is_projective(self):
        r"""
        Test  whether the code is projective.

        A linear code `C` over a field is called *projective* when its dual `Cd`
        has minimum weight `\geq 3`, i.e. when no two coordinate positions of
        `C` are linearly independent (cf. definition 3 from [BS2011]_ or 9.8.1 from
        [BH12]_).

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2), False)
            sage: C.is_projective()
            True
            sage: C.dual_code().minimum_distance()
            8

        A non-projective code::

            sage: C = codes.LinearCode(matrix(GF(2),[[1,0,1],[1,1,1]]))
            sage: C.is_projective()
            False
        """
        M = self.generator_matrix().transpose()

        def projectivize(row):
            if not row.is_zero():
                for i in range(len(row)):
                    if row[i]:
                        break
                row = ~(row[i]) * row
            row.set_immutable()
            return row

        rows = set()
        for row in M.rows():
            row = projectivize(row)
            if row in rows:
                return False
            rows.add(row)

        return True

    def dual_code(self):
        r"""
        Returns the dual code `C^{\perp}` of the code `C`,

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
        return LinearCode(self.parity_check_matrix())

    def dimension(self):
        r"""
        Returns the dimension of this code.

        EXAMPLES::

            sage: G = matrix(GF(2),[[1,0,0],[1,1,0]])
            sage: C = LinearCode(G)
            sage: C.dimension()
            2

        TESTS:

        Check that :trac:`21156` is fixed::

            sage: from sage.coding.linear_code import AbstractLinearCode
            sage: from sage.coding.encoder import Encoder
            sage: class MonkeyCode(AbstractLinearCode):
            ....:     _registered_encoders = {}
            ....:     _registered_decoders = {}
            ....:     def __init__(self):
            ....:         super(MonkeyCode, self).__init__(GF(5), 10, "Monkey", "Syndrome")
            ....:
            sage: class MonkeyEncoder(Encoder):
            ....:     def __init__(self, code):
            ....:         super(MonkeyEncoder, self).__init__(C)
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

    def direct_sum(self, other):
        """
        Direct sum of the codes ``self`` and ``other``

        Returns the code given by the direct sum of the codes ``self`` and
        ``other``, which must be linear codes defined over the same base ring.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: C2 = C1.direct_sum(C1); C2
            [14, 8] linear code over GF(2)
            sage: C3 = C1.direct_sum(C2); C3
            [21, 12] linear code over GF(2)
        """
        C1 = self; C2 = other
        G1 = C1.generator_matrix()
        G2 = C2.generator_matrix()
        F = C1.base_ring()
        n1 = len(G1.columns())
        k1 = len(G1.rows())
        n2 = len(G2.columns())
        k2 = len(G2.rows())
        MS1 = MatrixSpace(F,k2,n1)
        MS2 = MatrixSpace(F,k1,n2)
        Z1 = MS1(0)
        Z2 = MS2(0)
        top = G1.augment(Z2)
        bottom = Z1.augment(G2)
        G = top.stack(bottom)
        return LinearCode(G)

    def juxtapose(self, other):
        """
        Juxtaposition of ``self`` and ``other``

        The two codes must have equal dimension.

        EXAMPLES::

           sage: C1 = codes.HammingCode(GF(2), 3)
           sage: C2 = C1.juxtapose(C1)
           sage: C2
           [14, 4] linear code over GF(2)
        """
        G1 = self.generator_matrix()
        G2 = other.generator_matrix()
        G = G1.augment(G2)
        return LinearCode(G)

    def u_u_plus_v_code(self, other):
        r"""
        The `(u|u+v)`-construction with ``self=u`` and ``other=v``

        Returns the code obtained through `(u|u+v)`-construction with ``self`` as `u`
        and ``other`` as `v`. Note that `u` and `v` must have equal lengths.
        For `u` a `[n, k_1, d_1]`-code and `v` a `[n, k_2, d_2]`-code this returns
        a `[2n, k_1+k_2, d]`-code, where `d=\min(2d_1,d_2)`.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: D = C1.u_u_plus_v_code(C2)
            sage: D
            [14, 8] linear code over GF(2)
        """
        F = self.base_ring()
        G1 = self.generator_matrix()
        G2 = other.generator_matrix()
        k2 = len(G2.rows())
        n2 = len(G2.columns())
        MS = MatrixSpace(F,k2,n2)
        Z = MS(0)
        top = G1.augment(G1)
        bot = Z.augment(G2)
        G = top.stack(bot)
        return LinearCode(G)

    def product_code(self, other):
        """
        Combines ``self`` with ``other`` to give the tensor product code.

        If ``self`` is a `[n_1, k_1, d_1]`-code and ``other`` is
        a `[n_2, k_2, d_2]`-code, the product is a `[n_1n_2, k_1k_2, d_1d_2]`-code.

        Note that the two codes have to be over the same field.

            EXAMPLES::

                sage: C = codes.HammingCode(GF(2), 3)
                sage: C
                [7, 4] Hamming Code over GF(2)
                sage: D = codes.ReedMullerCode(GF(2), 2, 2)
                sage: D
                Binary Reed-Muller Code of order 2 and number of variables 2
                sage: A = C.product_code(D)
                sage: A
                [28, 16] linear code over GF(2)
                sage: A.length() == C.length()*D.length()
                True
                sage: A.dimension() == C.dimension()*D.dimension()
                True
                sage: A.minimum_distance() == C.minimum_distance()*D.minimum_distance()
                True

        """
        G1 = self.generator_matrix()
        G2 = other.generator_matrix()
        G = G1.tensor_product(G2)
        return LinearCode(G)

    def construction_x(self, other, aux):
        r"""
        Construction X applied to ``self=C_1``, ``other=C_2`` and ``aux=C_a``.

        ``other`` must be a subcode of ``self``.

        If `C_1` is a `[n, k_1, d_1]` linear code and `C_2` is
        a `[n, k_2, d_2]` linear code, then `k_1 > k_2` and `d_1 < d_2`. `C_a` must
        be a `[n_a, k_a, d_a]` linear code, such that `k_a + k_2 = k_1`
        and `d_a + d_1 \leq d_2`.

        The method will then return a `[n+n_a, k_1, d_a+d_1]` linear code.

            EXAMPLES::

                sage: C = codes.BCHCode(GF(2),15,7)
                sage: C
                [15, 5] BCH Code over GF(2) with designed distance 7
                sage: D = codes.BCHCode(GF(2),15,5)
                sage: D
                [15, 7] BCH Code over GF(2) with designed distance 5
                sage: C.is_subcode(D)
                True
                sage: C.minimum_distance()
                7
                sage: D.minimum_distance()
                5
                sage: aux = codes.HammingCode(GF(2),2)
                sage: aux = aux.dual_code()
                sage: aux.minimum_distance()
                2
                sage: Cx = D.construction_x(C,aux)
                sage: Cx
                [18, 7] linear code over GF(2)
                sage: Cx.minimum_distance()
                7
        """
        if other.is_subcode(self) == False:
            raise ValueError("%s is not a subcode of %s"%(self,other))

        G2 = self.generator_matrix()
        left = other.generator_matrix()  # G1
        k = self.dimension()

        for r in G2.rows():
            if r not in left.row_space():
                left = left.stack(r)

        Ga = aux.generator_matrix()
        na = aux.length()
        ka = aux.dimension()

        F = self.base_field()
        MS = MatrixSpace(F,k-ka,na)
        Z = MS(0)
        right = Z.stack(Ga)
        G = left.augment(right)
        return LinearCode(G)

    def __eq__(self, right):
        """
        Checks if ``self`` is equal to ``right``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: C1 == C2
            True

        TESTS:

        We check that :trac:`16644` is fixed::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C == ZZ
            False
        """
        if not (isinstance(right, LinearCode)\
                and self.length() == right.length()\
                and self.dimension() == right.dimension()\
                and self.base_ring() == right.base_ring()):
            return False
        Ks = self.parity_check_matrix().right_kernel()
        rbas = right.gens()
        if not all(c in Ks for c in rbas):
            return False
        Kr = right.parity_check_matrix().right_kernel()
        sbas = self.gens()
        if not all(c in Kr for c in sbas):
            return False
        return True

    def __ne__(self, other):
        r"""
        Tests inequality of ``self`` and ``other``.

        This is a generic implementation, which returns the inverse of ``__eq__`` for self.

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

    def encode(self, word, encoder_name=None, *args, **kwargs):
        r"""
        Transforms an element of a message space into a codeword.

        INPUT:

        - ``word`` -- a vector of a message space of the code.

        - ``encoder_name`` -- (default: ``None``) Name of the encoder which will be used
          to encode ``word``. The default encoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        One can use the following shortcut to encode a word ::

            C(word)

        OUTPUT:

        - a vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word)
            (1, 1, 0, 0, 1, 1, 0)
            sage: C(word)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the encoder amongst the list of the available ones::

            sage: sorted(C.encoders_available())
            ['GeneratorMatrix', 'Systematic']
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word, 'GeneratorMatrix')
            (1, 1, 0, 0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, *args, **kwargs)
        return E.encode(word)

    def __call__(self, m):
        r"""
        Returns either ``m`` if it is a codeword or ``self.encode(m)``
        if it is an element of the message space of the encoder used by
        ``encode``.

        INPUT:

        - ``m`` -- a vector whose length equals to code's length or an element
          of the message space used by ``encode``

        - ``**kwargs`` -- extra arguments are forwarded to ``encode``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0))
            sage: C(word)
            (1, 1, 0, 0, 1, 1, 0)

            sage: c = C.random_element()
            sage: C(c) == c
            True

        TESTS:

        If one passes a vector which belongs to the ambient space, it has to be a codeword.
        Otherwise, an exception is raised::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0, 0, 1, 0))
            sage: C(word)
            Traceback (most recent call last):
            ...
            ValueError: If the input is a vector which belongs to the ambient space, it has to be a codeword
        """
        if m in self.ambient_space():
            if m in self:
                return m
            else:
                raise ValueError("If the input is a vector which belongs to the ambient space, it has to be a codeword")
        else:
            return self.encode(m)

    @cached_method
    def encoder(self, encoder_name=None, *args, **kwargs):
        r"""
        Returns an encoder of ``self``.

        The returned encoder provided by this method is cached.

        This methods creates a new instance of the encoder subclass designated by ``encoder_name``.
        While it is also possible to do the same by directly calling the subclass' constructor,
        it is strongly advised to use this method to take advantage of the caching mechanism.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          returned. The default encoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the constructor of the encoder
          this method will return.

        OUTPUT:

        - an Encoder object.

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoder()
            Generator matrix-based encoder for [7, 4] linear code over GF(2)

        We check that the returned encoder is cached::

            sage: C.encoder.is_in_cache()
            True

        If the name of an encoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: sorted(C.encoders_available())
            ['GeneratorMatrix', 'Systematic']
            sage: C.encoder('NonExistingEncoder')
            Traceback (most recent call last):
            ...
            ValueError: There is no Encoder named 'NonExistingEncoder'. The known Encoders are: ['GeneratorMatrix', 'Systematic']

        Some encoders take extra arguments. If the user incorrectly supplies
        these, the error message attempts to be helpful::

            sage: C.encoder('Systematic', strange_parameter=True)
            Traceback (most recent call last):
            ...
            ValueError: Constructing the Systematic encoder failed, possibly due to missing or incorrect parameters.
            The constructor requires no arguments.
            It takes the optional arguments ['systematic_positions'].
            See the documentation of sage.coding.linear_code.LinearCodeSystematicEncoder for more details.
        """
        if encoder_name is None:
            encoder_name = self._default_encoder_name
        if encoder_name in self._registered_encoders:
            encClass = self._registered_encoders[encoder_name]
            try:
                return encClass(self, *args, **kwargs)
            except TypeError:
                raise ValueError(
                        "Constructing the {0} encoder failed, possibly due "
                        "to missing or incorrect parameters.\n{1}".format(
                            encoder_name, _explain_constructor(encClass)))
        else:
            raise ValueError(
                    "There is no Encoder named '{0}'. The known Encoders are: "
                    "{1}".format(encoder_name, self.encoders_available()))

    def encoders_available(self, classes=False):
        r"""
        Returns a list of the available encoders' names for ``self``.

        INPUT:

        - ``classes`` -- (default: ``False``) if ``classes`` is set to ``True``,
          return instead a ``dict`` mapping available encoder name to the
          associated encoder class.

        OUTPUT: a list of strings, or a `dict` mapping strings to classes.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoders_available()
            ['GeneratorMatrix', 'Systematic']
            sage: dictionary = C.encoders_available(True)
            sage: sorted(dictionary.items())
            [('GeneratorMatrix', <class 'sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder'>),
             ('Systematic', <class 'sage.coding.linear_code.LinearCodeSystematicEncoder'>)]
        """
        if classes:
            return copy(self._registered_encoders)

        return sorted(self._registered_encoders)

    def extended_code(self):
        r"""
        Returns `self` as an extended code.

        See documentation of :class:`sage.coding.extended_code.ExtendedCode`
        for details.
        EXAMPLES::


            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: C
            [21, 18] Hamming Code over GF(4)
            sage: Cx = C.extended_code()
            sage: Cx
            Extension of [21, 18] Hamming Code over GF(4)
        """
        from .extended_code import ExtendedCode
        return ExtendedCode(self)

    def galois_closure(self, F0):
        r"""
        If ``self`` is a linear code defined over `F` and `F_0` is a subfield
        with Galois group `G = Gal(F/F_0)` then this returns the `G`-module
        `C^-` containing `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: Cc = C.galois_closure(GF(2))
            sage: C; Cc
            [21, 18] Hamming Code over GF(4)
            [21, 20] linear code over GF(4)
            sage: c = C.basis()[2]
            sage: V = VectorSpace(GF(4,'a'),21)
            sage: c2 = V([x^2 for x in c.list()])
            sage: c2 in C
            False
            sage: c2 in Cc
            True
        """
        G = self.generator_matrix()
        F = self.base_ring()
        q = F.order()
        q0 = F0.order()
        a = q.log(q0)  # test if F/F0 is a field extension
        if not isinstance(a, Integer):
            raise ValueError("Base field must be an extension of given field %s"%F0)
        n = len(G.columns())
        k = len(G.rows())
        G0 = [[x**q0 for x in g.list()] for g in G.rows()]
        G1 = [[x for x in g.list()] for g in G.rows()]
        G2 = G0+G1
        MS = MatrixSpace(F,2*k,n)
        G3 = MS(G2)
        r = G3.rank()
        MS = MatrixSpace(F,r,n)
        Grref = G3.echelon_form()
        G = MS([Grref.row(i) for i in range(r)])
        return LinearCode(G)

    def __getitem__(self, i):
        r"""
        Returns the `i`-th codeword of this code.

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
         Hence the `i`-th element can be obtained by the p-adic expansion
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

        Check that codewords are immutable. See :trac:`16338`::

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

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Returns a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. The default encoder of ``self``
          will be used if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix()
            [1 2 1]
            [2 1 1]
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.generator_matrix()

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

        Specific systematic positions can also be requested:

            sage: C.systematic_generator_matrix(systematic_positions=[3,2])
            [1 2 0 1]
            [1 2 1 0]
        """
        systematic_positions = tuple(systematic_positions) if systematic_positions else None
        return self.encoder("Systematic", systematic_positions=systematic_positions).generator_matrix()

    @cached_method
    def gens(self):
        r"""
        Returns the generators of this code as a list of vectors.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.gens()
             [(1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1), (0, 0, 1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 1, 1)]
        """
        return self.generator_matrix().rows()

    def genus(self):
        r"""
        Returns the "Duursma genus" of the code, `\gamma_C = n+1-k-d`.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3); C1
            [7, 4] Hamming Code over GF(2)
            sage: C1.genus()
            1
            sage: C2 = codes.HammingCode(GF(4,"a"), 2); C2
            [5, 3] Hamming Code over GF(4)
            sage: C2.genus()
            0

        Since all Hamming codes have minimum distance 3, these computations
        agree with the definition, `n+1-k-d`.
        """
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        gammaC = n+1-k-d
        return gammaC

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

    @cached_method
    def information_set(self):
        """
        Return an information set of the code.

        Return value of this method is cached.

        A set of column positions of a generator matrix of a code
        is called an information set if the corresponding columns
        form a square matrix of full rank.

        OUTPUT:

        - Information set of a systematic generator matrix of the code.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,2,0,\
                                      2,1,1])
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
          is the length of `self`.

        OUTPUT:

        - A boolean indicating whether the positions form an information set.


        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,2,0,\
                                      2,1,1])
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

    def is_permutation_automorphism(self,g):
        r"""
        Returns `1` if `g` is an element of `S_n` (`n` = length of self) and
        if `g` is an automorphism of self.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(3), 3)
            sage: g = SymmetricGroup(13).random_element()
            sage: C.is_permutation_automorphism(g)
            0
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
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

    def is_permutation_equivalent(self,other,algorithm=None):
        """
        Returns ``True`` if ``self`` and ``other`` are permutation equivalent
        codes and ``False`` otherwise.

        The ``algorithm="verbose"`` option also returns a permutation (if
        ``True``) sending ``self`` to ``other``.

        Uses Robert Miller's double coset partition refinement work.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(2),"x")
            sage: g = x^3+x+1
            sage: C1 = codes.CyclicCode(length = 7, generator_pol = g); C1
            [7, 4] Cyclic Code over GF(2)
            sage: C2 = codes.HammingCode(GF(2), 3); C2
            [7, 4] Hamming Code over GF(2)
            sage: C1.is_permutation_equivalent(C2)
            True
            sage: C1.is_permutation_equivalent(C2,algorithm="verbose")
            (True, (3,4)(5,7,6))
            sage: C1 = codes.random_linear_code(GF(2), 10, 5)
            sage: C2 = codes.random_linear_code(GF(3), 10, 5)
            sage: C1.is_permutation_equivalent(C2)
            False
        """
        from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct
        F = self.base_ring()
        F_o = other.base_ring()
        q = F.order()
        G = self.generator_matrix()
        n = self.length()
        n_o = other.length()
        if F != F_o or n != n_o:
            return False
        k = len(G.rows())
        MS = MatrixSpace(F,q**k,n)
        CW1 = MS(self.list())
        CW2 = MS(other.list())
        B1 = NonlinearBinaryCodeStruct(CW1)
        B2 = NonlinearBinaryCodeStruct(CW2)
        ans = B1.is_isomorphic(B2)
        if ans is not False:
            if algorithm=="verbose":
                Sn = SymmetricGroup(n)
                return True, Sn([i+1 for i in ans])**(-1)
            return True
        return False

    def is_self_dual(self):
        """
        Returns ``True`` if the code is self-dual (in the usual Hamming inner
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
        Returns ``True`` if this code is self-orthogonal and ``False``
        otherwise.

        A code is self-orthogonal if it is a subcode of its dual.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: C.is_self_orthogonal()
            True
            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.is_self_orthogonal()
            False
            sage: C = codes.QuasiQuadraticResidueCode(11)  # optional - gap_packages (Guava package)
            sage: C.is_self_orthogonal()             # optional - gap_packages (Guava package)
            True
        """
        return self.is_subcode(self.dual_code())

    def is_galois_closed(self):
        r"""
        Checks if ``self`` is equal to its Galois closure.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,"a"), 3)
            sage: C.is_galois_closed()
            False
        """
        F = self.base_ring()
        p = F.characteristic()
        return self == self.galois_closure(GF(p))

    def is_subcode(self, other):
        """
        Returns ``True`` if ``self`` is a subcode of ``other``.

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
        for r in G.rows():
            if not(r in other):
                return False
        return True

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

    def length(self):
        r"""
        Returns the length of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.length()
            7
        """
        return self._length

    def list(self):
        r"""
        Return a list of all elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 1, 0, 1)
            True
        """
        return [x for x in self]

    def _magma_init_(self, magma):
        r"""
        Retun a string representation in Magma of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Cm = magma(C)                 # optional - magma, indirect doctest
            sage: Cm.MinimumWeight()            # optional - magma
            3

        """
        G = magma(self.generator_matrix())._ref()
        s = 'LinearCode(%s)' % G
        return s

    @cached_method
    def minimum_distance(self, algorithm=None):
        r"""
        Returns the minimum distance of ``self``.

        .. NOTE::

            When using GAP, this raises a ``NotImplementedError`` if
            the base field of the code has size greater than 256 due
            to limitations in GAP.

        INPUT:

        -  ``algorithm`` -- (default: ``None``) the name of the algorithm to use
           to perform minimum distance computation. If set to ``None``,
           GAP methods will be used. ``algorithm`` can be:
           - ``"Guava"``, which will use optional GAP package Guava

        OUTPUT:

        - Integer, minimum distance of this code

        EXAMPLES::

            sage: MS = MatrixSpace(GF(3),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            3

        If ``algorithm`` is provided, then the minimum distance will be
        recomputed even if there is a stored value from a previous run.::

            sage: C.minimum_distance(algorithm="gap")
            3
            sage: C.minimum_distance(algorithm="guava")  # optional - gap_packages (Guava package)
            3

        TESTS::

            sage: C = codes.random_linear_code(GF(4,"a"), 5, 2)
            sage: C.minimum_distance(algorithm='something')
            Traceback (most recent call last):
            ...
            ValueError: The algorithm argument must be one of None, 'gap' or 'guava'; got 'something'

        The field must be size at most 256::

            sage: C = codes.random_linear_code(GF(257,"a"), 5, 2)
            sage: C.minimum_distance()
            Traceback (most recent call last):
            ...
            NotImplementedError: the GAP algorithm that Sage is using
             is limited to computing with fields of size at most 256
        """
        if algorithm == "guava":
            GapPackage("guava", spkg="gap_packages").require()

        # If the minimum distance has already been computed or provided by
        # the user then simply return the stored value.
        # This is done only if algorithm is None.
        if algorithm not in (None, "gap", "guava"):
            raise ValueError("The algorithm argument must be one of None, "
                        "'gap' or 'guava'; got '{0}'".format(algorithm))

        F = self.base_ring()
        q = F.order()
        if q > 256:
            raise NotImplementedError("the GAP algorithm that Sage is using "
                                      "is limited to computing with fields "
                                      "of size at most 256")

        G = self.generator_matrix()
        if (q == 2 or q == 3) and algorithm=="guava":
            gap.load_package("guava")
            C = gap(G).GeneratorMatCode(gap(F))
            d = C.MinimumWeight()
            return ZZ(d)
        return self._minimum_weight_codeword(algorithm).hamming_weight()

    def _minimum_weight_codeword(self, algorithm = None):
        r"""
        Returns a minimum weight codeword of ``self``.

        INPUT:

        -  ``algorithm`` -- (default: ``None``) the name of the algorithm to use
           to perform minimum weight codeword search. If set to ``None``,
           a search using GAP methods will be done. ``algorithm`` can be:
           - ``"Guava"``, which will use optional GAP package Guava

        REMARKS:

        - The code in the default case allows one (for free) to also compute the
          message vector `m` such that `m\*G = v`, and the (minimum) distance, as
          a triple.  however, this output is not implemented.
        - The binary case can presumably be done much faster using Robert Miller's
          code (see the docstring for the spectrum method). This is also not (yet)
          implemented.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C._minimum_weight_codeword()
            (0, 1, 0, 1, 0, 1, 0)

        TESTS:

        We check that :trac:`18480` is fixed::

            sage: codes.HammingCode(GF(2), 2).minimum_distance()
            3

        AUTHORS:

        - David Joyner (11-2005)
        """
        n, k = self.length(), self.dimension()
        F = self.base_field()
        Gmat = self.generator_matrix()._gap_init_()

        current_randstate().set_seed_gap()

        if algorithm=="guava":
            GapPackage("guava", spkg="gap_packages").require()
            gap.load_package("guava")
            from sage.interfaces.gap import gfq_gap_to_sage
            gap.eval("G:="+Gmat)
            C = gap(Gmat).GeneratorMatCode(F)
            cg = C.MinimumDistanceCodeword()
            c = [gfq_gap_to_sage(cg[j],F) for j in range(1,n+1)]
            return vector(F, c)

        q = F.order()
        ans = None
        dist_min = n + 1
        gap.eval('Gmat:='+Gmat)
        gap.eval('K:=GF({})'.format(q))
        gap.eval('v:=Z({})*{}'.format(q,[0]*n))
        for i in range(1,k+1):
            gap.eval("P:=AClosestVectorCombinationsMatFFEVecFFECoords(Gmat,K,v,{},1)".format(i))
            gap.eval("d:=WeightVecFFE(P[1])")
            v = gap("P[1]")
            dist = gap("d")
            if dist and dist < dist_min:
                dist_min = dist
                ans = list(v)

        if ans is None:
            raise RuntimeError("Computation failed due to some GAP error")

        # return the result as a vector (and not a 1xn matrix)
        return vector(F, ans)

    def module_composition_factors(self, gp):
        r"""
        Prints the GAP record of the Meataxe composition factors module in
        Meataxe notation. This uses GAP but not Guava.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: gp = C.permutation_automorphism_group()

        Now type "C.module_composition_factors(gp)" to get the record printed.
        """
        F = self.base_ring()
        q = F.order()
        gens = gp.gens()
        G = self.generator_matrix()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        mats = [] # initializing list of mats by which the gens act on self
        Fn = VectorSpace(F, n)
        W = Fn.subspace_with_basis(G.rows()) # this is self
        for g in gens:
            p = MS(g.matrix())
            m = [W.coordinate_vector(r*p) for r in G.rows()]
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("M:=GModuleByMats("+mats_str+", GF("+str(q)+"))")
        print(gap("MTX.CompositionFactors( M )"))

    def permutation_automorphism_group(self, algorithm="partition"):
        r"""
        If `C` is an `[n,k,d]` code over `F`, this function computes the
        subgroup `Aut(C) \subset S_n` of all permutation automorphisms of `C`.
        The binary case always uses the (default) partition refinement
        algorithm of Robert Miller.

        Note that if the base ring of `C` is `GF(2)` then this is the full
        automorphism group. Otherwise, you could use
        :meth:`~sage.coding.linear_code.LinearCode.automorphism_group_gens`
        to compute generators of the full automorphism group.

        INPUT:

        - ``algorithm`` - If ``"gap"`` then GAP's MatrixAutomorphism function
          (written by Thomas Breuer) is used. The implementation combines an
          idea of mine with an improvement suggested by Cary Huffman. If
          ``"gap+verbose"`` then code-theoretic data is printed out at
          several stages of the computation. If ``"partition"`` then the
          (default) partition refinement algorithm of Robert Miller is used.
          Finally, if ``"codecan"`` then the partition refinement algorithm
          of Thomas Feulner is used, which also computes a canonical
          representative of ``self`` (call
          :meth:`~sage.coding.linear_code.LinearCode.canonical_representative`
          to access it).

        OUTPUT:

        - Permutation automorphism group

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: C
            [8, 4] linear code over GF(2)
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            144
            sage: GG = C.permutation_automorphism_group("codecan")
            sage: GG == G
            True

        A less easy example involves showing that the permutation
        automorphism group of the extended ternary Golay code is the
        Mathieu group `M_{11}`.

        ::

            sage: C = codes.GolayCode(GF(3))
            sage: M11 = MathieuGroup(11)
            sage: M11.order()
            7920
            sage: G = C.permutation_automorphism_group()  # long time (6s on sage.math, 2011)
            sage: G.is_isomorphic(M11)                    # long time
            True
            sage: GG = C.permutation_automorphism_group("codecan") # long time
            sage: GG == G # long time
            True

        Other examples::

            sage: C = codes.GolayCode(GF(2))
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            244823040
            sage: C = codes.HammingCode(GF(2), 5)
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            9999360
            sage: C = codes.HammingCode(GF(3), 2); C
            [4, 2] Hamming Code over GF(3)
            sage: C.permutation_automorphism_group(algorithm="partition")
            Permutation Group with generators [(1,3,4)]
            sage: C = codes.HammingCode(GF(4,"z"), 2); C
            [5, 3] Hamming Code over GF(4)
            sage: G = C.permutation_automorphism_group(algorithm="partition"); G
            Permutation Group with generators [(1,3)(4,5), (1,4)(3,5)]
            sage: GG = C.permutation_automorphism_group(algorithm="codecan") # long time
            sage: GG == G # long time
            True
            sage: C.permutation_automorphism_group(algorithm="gap")  # optional - gap_packages (Guava package)
            Permutation Group with generators [(1,3)(4,5), (1,4)(3,5)]
            sage: C = codes.GolayCode(GF(3), True)
            sage: C.permutation_automorphism_group(algorithm="gap")  # optional - gap_packages (Guava package)
            Permutation Group with generators [(5,7)(6,11)(8,9)(10,12), (4,6,11)(5,8,12)(7,10,9), (3,4)(6,8)(9,11)(10,12), (2,3)(6,11)(8,12)(9,10), (1,2)(5,10)(7,12)(8,9)]

        However, the option ``algorithm="gap+verbose"``, will print out::

            Minimum distance: 5 Weight distribution: [1, 0, 0, 0, 0, 132, 132,
            0, 330, 110, 0, 24]

            Using the 132 codewords of weight 5 Supergroup size: 39916800

        in addition to the output of
        ``C.permutation_automorphism_group(algorithm="gap")``.
        """
        F = self.base_ring()
        q = F.order()
        G = self.generator_matrix() if 2*self.dimension() <= self.length() else self.dual_code().generator_matrix()
        n = len(G.columns())
        if "gap" in algorithm:
            GapPackage("guava", spkg="gap_packages").require()
            gap.load_package('guava')
            wts = self.weight_distribution()                          # bottleneck 1
            nonzerowts = [i for i in range(len(wts)) if wts[i]!=0]
            Sn = SymmetricGroup(n)
            Gp = gap("SymmetricGroup(%s)"%n)               # initializing G in gap
            Gstr = str(gap(G))
            gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
            gap.eval("eltsC:=Elements(C)")
            if algorithm=="gap+verbose":
                print("\n Minimum distance: %s \n Weight distribution: \n %s" % (nonzerowts[1], wts))
            stop = 0                                          # only stop if all gens are autos
            for i in range(1,len(nonzerowts)):
                if stop == 1:
                    break
                wt = nonzerowts[i]
                if algorithm=="gap+verbose":
                    size = Gp.Size()
                    print("\n Using the %s codewords of weight %s \n Supergroup size: \n %s\n " % (wts[wt], wt, size))
                gap.eval("Cwt:=Filtered(eltsC,c->WeightCodeword(c)=%s)"%wt)   # bottleneck 2 (repeated
                gap.eval("matCwt:=List(Cwt,c->VectorCodeword(c))")            # for each i until stop = 1)
                if gap("Length(matCwt)") > 0:
                    A = gap("MatrixAutomorphisms(matCwt)")
                    G2 = gap("Intersection2(%s,%s)"%(str(A).replace("\n",""),str(Gp).replace("\n",""))) #  bottleneck 3
                    Gp = G2
                    if Gp.Size()==1:
                        return PermutationGroup([()])
                    autgp_gens = Gp.GeneratorsOfGroup()
                    gens = [Sn(str(x).replace("\n","")) for x in autgp_gens]
                    stop = 1                    # get ready to stop
                    for x in gens:              # if one of these gens is not an auto then don't stop
                        if not(self.is_permutation_automorphism(x)):
                            stop = 0
                            break
            G = PermutationGroup(gens)
            return G
        if algorithm=="partition":
            if q == 2:
                from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct
                B = LinearBinaryCodeStruct(G)
                autgp = B.automorphism_group()
                L = [[j+1 for j in gen] for gen in autgp[0]]
                AutGp = PermutationGroup(L)
            else:
                from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
                from sage.matrix.constructor import matrix
                weights = {}
                for c in self:
                    wt = c.hamming_weight()
                    if wt not in weights:
                        weights[wt] = [c]
                    else:
                        weights[wt].append(c)
                weights.pop(0)
                AutGps = []
                for wt, words in iteritems(weights):
                    M = MatrixStruct(matrix(words))
                    autgp = M.automorphism_group()
                    L = [[j+1 for j in gen] for gen in autgp[0]]
                    G = PermutationGroup(L)
                    AutGps.append(G)
                if len(AutGps) > 0:
                    AutGp = AutGps[0]
                    for G in AutGps[1:]:
                        AutGp = AutGp.intersection(G)
                else:
                    return PermutationGroup([])
            return AutGp
        if algorithm=="codecan":
            gens, _ = self.automorphism_group_gens("permutational")
            return PermutationGroup([x.get_perm() for x in gens])
        raise NotImplementedError("The only algorithms implemented currently are 'gap', 'gap+verbose', and 'partition'.")

    def permuted_code(self, p):
        r"""
        Returns the permuted code, which is equivalent to ``self`` via the
        column permutation ``p``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: G = C.permutation_automorphism_group(); G
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: g = G("(2,3)(6,7)")
            sage: Cg = C.permuted_code(g)
            sage: Cg
            [7, 4] linear code over GF(2)
            sage: C.generator_matrix() == Cg.systematic_generator_matrix()
            True
        """
        G = copy(self.generator_matrix())
        G.permute_columns(p)
        return LinearCode(G)

    def punctured(self, L):
        r"""
        Returns a :class:`sage.coding.punctured_code` object from ``L``.

        INPUT:

        - ``L`` - List of positions to puncture

        OUTPUT:

        - an instance of :class:`sage.coding.punctured_code`

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.punctured([1,2])
            Puncturing of [7, 4] Hamming Code over GF(2) on position(s) [1, 2]
        """
        from .punctured_code import PuncturedCode
        return PuncturedCode(self, set(L))

    def _punctured_form(self, points):
        r"""
        Returns a representation of self as a :class:`LinearCode` punctured in ``points``.

        INPUT:

        - ``points`` -- a set of positions where to puncture ``self``

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 4)
            sage: C._punctured_form({3})
            [10, 4] linear code over GF(7)
        """
        if not isinstance(points, (Integer, int, set)):
            raise TypeError("points must be either a Sage Integer, a Python int, or a set")
        M = self.generator_matrix()
        G = M.delete_columns(list(points))
        G = G.echelon_form()
        k = G.rank()
        return LinearCode(G[:k])

    def random_element(self, *args, **kwds):
        """
        Returns a random codeword; passes other positional and keyword
        arguments to ``random_element()`` method of vector space.

        OUTPUT:

        - Random element of the vector space of this code

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: C.random_element() # random test
            (1, 0, 0, a + 1, 1, a, a, a + 1, a + 1, 1, 1, 0, a + 1, a, 0, a, a, 0, a, a, 1)

        Passes extra positional or keyword arguments through::

            sage: C.random_element(prob=.5, distribution='1/n') # random test
            (1, 0, a, 0, 0, 0, 0, a + 1, 0, 0, 0, 0, 0, 0, 0, 0, a + 1, a + 1, 1, 0, 0)

        TESTS:

        Test that the codeword returned is immutable (see :trac:`16469`)::

            sage: c = C.random_element()
            sage: c.is_immutable()
            True

        Test that codeword returned has the same parent as any non-random codeword
        (see :trac:`19653`)::

            sage: C = codes.random_linear_code(GF(16, 'a'), 10, 4)
            sage: c1 = C.random_element()
            sage: c2 = C[1]
            sage: c1.parent() == c2.parent()
            True
        """
        E = self.encoder()
        M = E.message_space()
        m = M.random_element(*args, **kwds)
        c = E.encode(m)
        c.set_immutable()
        return c

    def relative_distance(self):
        r"""
        Return the ratio of the minimum distance to the code length.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2),3)
            sage: C.relative_distance()
            3/7
        """
        return self.minimum_distance() / self.length()

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

    def redundancy_matrix(self):
        r"""
        Returns the non-identity columns of a systematic generator matrix for
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

    def shortened(self, L):
        r"""
        Returns the code shortened at the positions ``L``, where
        `L \subset \{1,2,...,n\}`.

        Consider the subcode `C(L)` consisting of all codewords `c\in C` which
        satisfy `c_i=0` for all `i\in L`. The punctured code `C(L)^L` is
        called the shortened code on `L` and is denoted `C_L`. The code
        constructed is actually only isomorphic to the shortened code defined
        in this way.

        By Theorem 1.5.7 in [HP2003]_, `C_L` is `((C^\perp)^L)^\perp`. This is used
        in the construction below.

        INPUT:

        - ``L`` - Subset of `\{1,...,n\}`, where `n` is the length of this code

        OUTPUT:

        - Linear code, the shortened code described above

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.shortened([1,2])
            [5, 2] linear code over GF(2)
        """
        Cd = self.dual_code()
        Cdp = Cd.punctured(set(L))
        return Cdp.dual_code()

    @cached_method
    def weight_distribution(self, algorithm=None):
        r"""
        Returns the weight distribution, or spectrum, of ``self`` as a list.

        The weight distribution a code of length `n` is the sequence `A_0,
        A_1,..., A_n` where `A_i` is the number of codewords of weight `i`.

        INPUT:

        - ``algorithm`` - (optional, default: ``None``) If set to ``"gap"``,
          call GAP. If set to `"leon"`, call the option GAP package GUAVA and
          call a function therein by Jeffrey Leon (see warning below). If set to
          ``"binary"``, use an algorithm optimized for binary codes. The default
          is to use ``"binary"`` for binary codes and ``"gap"`` otherwise.

        OUTPUT:

        - A list of non-negative integers: the weight distribution.

        .. WARNING::

            Specifying ``algorithm = "leon"`` sometimes prints a traceback
            related to a stack smashing error in the C library. The result
            appears to be computed correctly, however. It appears to run much
            faster than the GAP algorithm in small examples and much slower than
            the GAP algorithm in larger examples.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G = MS([[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.weight_distribution()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: F.<z> = GF(2^2,"z")
            sage: C = codes.HammingCode(F, 2); C
            [5, 3] Hamming Code over GF(4)
            sage: C.weight_distribution()
            [1, 0, 0, 30, 15, 18]
            sage: C = codes.HammingCode(GF(2), 3); C
            [7, 4] Hamming Code over GF(2)
            sage: C.weight_distribution(algorithm="leon")   # optional - gap_packages (Guava package)
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.weight_distribution(algorithm="gap")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.weight_distribution(algorithm="binary")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C = codes.HammingCode(GF(3), 3); C
            [13, 10] Hamming Code over GF(3)
            sage: C.weight_distribution() == C.weight_distribution(algorithm="leon")   # optional - gap_packages (Guava package)
            True
            sage: C = codes.HammingCode(GF(5), 2); C
            [6, 4] Hamming Code over GF(5)
            sage: C.weight_distribution() == C.weight_distribution(algorithm="leon")   # optional - gap_packages (Guava package)
            True
            sage: C = codes.HammingCode(GF(7), 2); C
            [8, 6] Hamming Code over GF(7)
            sage: C.weight_distribution() == C.weight_distribution(algorithm="leon")   # optional - gap_packages (Guava package)
            True

        """
        if algorithm is None:
            if self.base_ring().order() == 2:
                algorithm = "binary"
            else:
                algorithm = "gap"
        F = self.base_ring()
        n = self.length()
        if algorithm=="gap":
            Gmat = self.generator_matrix()._gap_init_()
            q = self.base_ring().order()
            z = 'Z(%s)*%s'%(q, [0]*self.length())     # GAP zero vector as a string
            _ = gap.eval("w:=DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
            v = [eval(gap.eval("w["+str(i)+"]")) for i in range(1,self.length()+2)] # because GAP returns vectors in compressed form
            return v
        elif algorithm=="binary":
            from sage.coding.binary_code import weight_dist
            return weight_dist(self.generator_matrix())
        elif algorithm=="leon":
            if not(F.order() in [2,3,5,7]):
                raise NotImplementedError("The algorithm 'leon' is only implemented for q = 2,3,5,7.")
            # The GAP command DirectoriesPackageLibrary tells the location of the latest
            # version of the Guava libraries, so gives us the location of the Guava binaries too.
            guava_bin_dir = gap.eval('DirectoriesPackagePrograms("guava")[1]')
            guava_bin_dir = guava_bin_dir[guava_bin_dir.index('"') + 1:guava_bin_dir.rindex('"')]
            input = _dump_code_in_leon_format(self) + "::code"
            import os
            import subprocess
            lines = subprocess.check_output([os.path.join(guava_bin_dir, 'wtdist'), input])
            from six import StringIO
            # to use the already present output parser
            wts = [0] * (n + 1)
            for L in StringIO(bytes_to_str(lines)).readlines():
                L = L.strip()
                if L:
                    o = ord(L[0])
                    if o >= 48 and o <= 57:
                        wt, num = L.split()
                        wts[eval(wt)] = eval(num)
            return wts
        else:
            raise NotImplementedError("The only algorithms implemented currently are 'gap', 'leon' and 'binary'.")

    spectrum = weight_distribution

    def standard_form(self, return_permutation=True):
        r"""
        Returns a linear code which is permutation-equivalent to ``self`` and
        admits a generator matrix in standard form.

        A generator matrix is in standard form if it is of the form `[I \vert
        A]`, where `I` is the `k \times k` identity matrix. Any code admits a
        generator matrix in systematic form, i.e. where a subset of the columns
        form the identity matrix, but one might need to permute columns to allow
        the identity matrix to be leading.

        INPUT:

        - ``return_permutation`` -- (default: ``True``) if ``True``, the column
          permutation which brings ``self`` into the returned code is also
          returned.

        OUTPUT:

        - A :class:`LinearCode` whose :meth:`systematic_generator_matrix` is
          guaranteed to be of the form `[I \vert A]`.

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

    def support(self):
        r"""
        Returns the set of indices `j` where `A_j` is nonzero, where
        `A_j` is the number of codewords in `self` of Hamming weight `j`.

        OUTPUT:

        - List of integers

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.weight_distribution()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.support()
            [0, 3, 4, 7]
        """
        n = self.length()
        F = self.base_ring()
        V = VectorSpace(F,n+1)
        return V(self.weight_distribution()).support()

    def syndrome(self, r):
        r"""
        Returns the syndrome of ``r``.

        The syndrome of ``r`` is the result of `H \times r` where `H` is
        the parity check matrix of ``self``. If ``r`` belongs to ``self``,
        its syndrome equals to the zero vector.

        INPUT:

        - ``r`` -- a vector of the same length as ``self``

        OUTPUT:

        - a column vector

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

    def unencode(self, c, encoder_name=None, nocheck=False, **kwargs):
        r"""
        Returns the message corresponding to ``c``.

        This is the inverse of :meth:`encode`.

        INPUT:

        - ``c`` -- a codeword of ``self``.

        - ``encoder_name`` -- (default: ``None``) name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``nocheck`` -- (default: ``False``) checks if ``c`` is in ``self``. You might set
          this to ``True`` to disable the check for saving computation. Note that if ``c`` is
          not in ``self`` and ``nocheck = True``, then the output of :meth:`unencode` is
          not defined (except that it will be in the message space of ``self``).

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        OUTPUT:

        - an element of the message space of ``encoder_name`` of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.unencode(c)
            (0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.unencode(c, nocheck)

    def weight_enumerator(self, names=None, bivariate=True):
        r"""
        Return the weight enumerator polynomial of ``self``.

        This is the bivariate, homogeneous polynomial in `x` and `y` whose
        coefficient to `x^i y^{n-i}` is the number of codewords of `self` of
        Hamming weight `i`. Here, `n` is the length of `self`.

        INPUT:

        - ``names`` - (default: ``"xy"``) The names of the variables in the
          homogeneous polynomial. Can be given as a single string of length 2,
          or a single string with a comma, or as a tuple or list of two strings.

        - ``bivariate`` - (default: `True`) Whether to return a bivariate,
          homogeneous polynomial or just a univariate polynomial. If set to
          ``False``, then ``names`` will be interpreted as a single variable
          name and default to ``"x"``.

        OUTPUT:

        - The weight enumerator polynomial over `\ZZ`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.weight_enumerator()
            x^7 + 7*x^4*y^3 + 7*x^3*y^4 + y^7
            sage: C.weight_enumerator(names="st")
            s^7 + 7*s^4*t^3 + 7*s^3*t^4 + t^7
            sage: C.weight_enumerator(names="var1, var2")
            var1^7 + 7*var1^4*var2^3 + 7*var1^3*var2^4 + var2^7
            sage: C.weight_enumerator(names=('var1', 'var2'))
            var1^7 + 7*var1^4*var2^3 + 7*var1^3*var2^4 + var2^7
            sage: C.weight_enumerator(bivariate=False)
            x^7 + 7*x^4 + 7*x^3 + 1

        An example of a code with a non-symmetrical weight enumerator::

            sage: C = codes.GolayCode(GF(3), extended=False)
            sage: C.weight_enumerator()
            24*x^11 + 110*x^9*y^2 + 330*x^8*y^3 + 132*x^6*y^5 + 132*x^5*y^6 + y^11
        """
        if names is None:
            if bivariate:
                names = "xy"
            else:
                names = "x"
        spec = self.weight_distribution()
        n = self.length()
        if bivariate:
            R = PolynomialRing(ZZ,2,names)
            x,y = R.gens()
            return sum(spec[i]*x**i*y**(n-i) for i in range(n+1))
        else:
            R = PolynomialRing(ZZ,names)
            x, = R.gens()
            return sum(spec[i]*x**i for i in range(n+1))

    @cached_method
    def zero(self):
        r"""
        Returns the zero vector of ``self``.

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

    def zeta_polynomial(self, name="T"):
        r"""
        Returns the Duursma zeta polynomial of this code.

        Assumes that the minimum distances of this code and its dual are
        greater than 1.  Prints a warning to ``stdout`` otherwise.

        INPUT:

        - ``name`` - String, variable name (default: ``"T"``)

        OUTPUT:

        - Polynomial over `\QQ`

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C = codes.databases.best_linear_code_in_guava(6,3,GF(2))  # optional - gap_packages (Guava package)
            sage: C.minimum_distance()                   # optional - gap_packages (Guava package)
            3
            sage: C.zeta_polynomial()                    # optional - gap_packages (Guava package)
            2/5*T^2 + 2/5*T + 1/5
            sage: C = codes.HammingCode(GF(2), 4)
            sage: C.zeta_polynomial()
            16/429*T^6 + 16/143*T^5 + 80/429*T^4 + 32/143*T^3 + 30/143*T^2 + 2/13*T + 1/13
            sage: F.<z> = GF(4,"z")
            sage: MS = MatrixSpace(F, 3, 6)
            sage: G = MS([[1,0,0,1,z,z],[0,1,0,z,1,z],[0,0,1,z,z,1]])
            sage: C = LinearCode(G)  # the "hexacode"
            sage: C.zeta_polynomial()
            1

        REFERENCES:

        - [Du2001]_
        """
        n = self.length()
        q = (self.base_ring()).order()
        d = self.minimum_distance()
        dperp = (self.dual_code()).minimum_distance()
        if d == 1 or dperp == 1:
            print("\n WARNING: There is no guarantee this function works when the minimum distance")
            print("            of the code or of the dual code equals 1.\n")
        RT = PolynomialRing(QQ,"%s"%name)
        R = PolynomialRing(QQ,3,"xy%s"%name)
        x,y,T = R.gens()
        we = self.weight_enumerator()
        A = R(we)
        #B = A(x+y,y,T)-(x+y)**n
        B = A(x,x+y,T)-(x+y)**n
        Bs = B.coefficients()
        Bs.reverse()
        b = [Bs[i]/binomial(n,i+d) for i in range(len(Bs))]
        r = n-d-dperp+2
        P_coeffs = []
        for i in range(len(b)):
           if i == 0:
              P_coeffs.append(b[0])
           if i == 1:
              P_coeffs.append(b[1] - (q+1)*b[0])
           if i>1:
              P_coeffs.append(b[i] - (q+1)*b[i-1] + q*b[i-2])
        P = sum([P_coeffs[i]*T**i for i in range(r+1)])
        return RT(P)/RT(P)(1)

    def zeta_function(self, name="T"):
        r"""
        Returns the Duursma zeta function of the code.

        INPUT:

        - ``name`` - String, variable name (default: ``"T"``)

        OUTPUT:

        - Element of `\QQ(T)`

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.zeta_function()
            (1/5*T^2 + 1/5*T + 1/10)/(T^2 - 3/2*T + 1/2)
        """
        P =  self.zeta_polynomial()
        q = (self.base_ring()).characteristic()
        RT = PolynomialRing(QQ,"%s"%name)
        T = RT.gen()
        return P/((1-T)*(1-q*T))


############################ linear codes python class ########################

class LinearCode(AbstractLinearCode):
    r"""
    Linear codes over a finite field or finite ring, represented using a
    generator matrix.

    This class should be used for arbitrary and unstructured linear codes. This
    means that basic operations on the code, such as the computation of the
    minimum distance, will use generic, slow algorithms.

    If you are looking for constructing a code from a more specific family, see
    if the family has been implemented by investigating `codes.<tab>`. These
    more specific classes use properties particular to that family to allow
    faster algorithms, and could also have family-specific methods.

    See :wikipedia:`Linear_code` for more information on unstructured linear codes.

    INPUT:

    - ``generator`` -- a generator matrix over a finite field (``G`` can be
      defined over a finite ring but the matrices over that ring must have
      certain attributes, such as ``rank``); or a code over a finite field

    - ``d`` -- (optional, default: ``None``) the minimum distance of the code

    .. NOTE::

        The veracity of the minimum distance ``d``, if provided, is not
        checked.

    EXAMPLES::

        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        [7, 4] linear code over GF(2)
        sage: C.base_ring()
        Finite Field of size 2
        sage: C.dimension()
        4
        sage: C.length()
        7
        sage: C.minimum_distance()
        3
        sage: C.spectrum()
        [1, 0, 0, 7, 7, 0, 0, 1]
        sage: C.weight_distribution()
        [1, 0, 0, 7, 7, 0, 0, 1]

    The minimum distance of the code, if known, can be provided as an
    optional parameter.::

        sage: C  = LinearCode(G, d=3)
        sage: C.minimum_distance()
        3

    Another example.::

        sage: MS = MatrixSpace(GF(5),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        [7, 4] linear code over GF(5)

    Providing a code as the parameter in order to "forget" its structure (see
    :trac:`20198`)::

        sage: C = codes.GeneralizedReedSolomonCode(GF(23).list(), 12)
        sage: LinearCode(C)
        [23, 12] linear code over GF(23)

    Another example::

        sage: C = codes.HammingCode(GF(7), 3)
        sage: C
        [57, 54] Hamming Code over GF(7)
        sage: LinearCode(C)
        [57, 54] linear code over GF(7)

    AUTHORS:

    - David Joyner (11-2005)
    - Charles Prior (03-2016): :trac:`20198`, LinearCode from a code
    """
    def __init__(self, generator, d=None):
        r"""
        See the docstring for :meth:`LinearCode`.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)    # indirect doctest
            sage: C
            [7, 4] linear code over GF(2)

        The minimum distance of the code, if known, can be provided as an
        optional parameter.::

            sage: C  = LinearCode(G, d=3)
            sage: C.minimum_distance()
            3

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: TestSuite(C).run()

        Check that it works even with input matrix with non full rank (see
        :trac:`17452`)::

            sage: K.<a> = GF(4)
            sage: G = matrix([[a, a + 1, 1, a + 1, 1, 0, 0],
            ....:             [0, a, a + 1, 1, a + 1, 1, 0],
            ....:             [0, 0, a, a + 1, 1, a + 1, 1],
            ....:             [a + 1, 0, 1, 0, a + 1, 1, a + 1],
            ....:             [a, a + 1, a + 1, 0, 0, a + 1, 1],
            ....:             [a + 1, a, a, 1, 0, 0, a + 1],
            ....:             [a, a + 1, 1, a + 1, 1, 0, 0]])
            sage: C = LinearCode(G)
            sage: C.basis()
            [
            (1, 0, 0, a + 1, 0, 1, 0),
            (0, 1, 0, 0, a + 1, 0, 1),
            (0, 0, 1, a, a + 1, a, a + 1)
            ]
            sage: C.minimum_distance()
            3

        We can construct a linear code directly from a vector space
            sage: VS = matrix(GF(2), [[1,0,1],\
                                      [1,0,1]]).row_space()
            sage: C = LinearCode(VS); C
            [3, 1] linear code over GF(2)

        Forbid the zero vector space (see :trac:`17452` and :trac:`6486`)::

            sage: G = matrix(GF(2), [[0,0,0]])
            sage: C = LinearCode(G)
            Traceback (most recent call last):
            ...
            ValueError: this linear code contains no non-zero vector
        """

        base_ring = generator.base_ring()
        if not base_ring.is_field():
            raise ValueError("'generator' must be defined on a field (not a ring)")

        try:
            basis = None
            if hasattr(generator,"nrows"): # generator matrix case
                if generator.rank() < generator.nrows():
                    basis = generator.row_space().basis()
            else:
                basis = generator.basis() # vector space etc. case
            if not basis is None:
                from sage.matrix.constructor import matrix
                generator = matrix(base_ring, basis)
                if generator.nrows() == 0:
                    raise ValueError("this linear code contains no non-zero vector")
        except AttributeError:
            # Assume input is an AbstractLinearCode, extract its generator matrix
            generator = generator.generator_matrix()

        super(LinearCode, self).__init__(base_ring, generator.ncols(), "GeneratorMatrix", "Syndrome")
        self._generator_matrix = generator
        self._dimension = generator.rank()
        self._minimum_distance = d

    def _repr_(self):
        r"""
        See the docstring for :meth:`LinearCode`.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C                     # indirect doctest
            [7, 4] linear code over GF(2)
        """
        R = self.base_ring()
        if R in Fields():
            return "[%s, %s] linear code over GF(%s)"%(self.length(), self.dimension(), R.cardinality())
        else:
            return "[%s, %s] linear code over %s"%(self.length(), self.dimension(), R)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: latex(C)
            [7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "[%s, %s]\\textnormal{ Linear code over }%s"\
                % (self.length(), self.dimension(), self.base_ring()._latex_())

    def __hash__(self):
        r"""
        Returns the hash value of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: hash(C) #random
            9015017528451745710

        If ``C1`` and ``C2`` are two codes which only differ by the
        coefficients of their generator matrices, their hashes are
        different (we check that the bug found in :trac:`18813` is
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

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Returns a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. ``self._generator_matrix``
          will be returned if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix()
            [1 2 1]
            [2 1 1]
        """
        if encoder_name is None or encoder_name == 'GeneratorMatrix':
            g = self._generator_matrix
        else:
            g = super(LinearCode, self).generator_matrix(encoder_name, **kwargs)
        g.set_immutable()
        return g


####################### encoders ###############################

class LinearCodeGeneratorMatrixEncoder(Encoder):
    r"""
    Encoder based on generator_matrix for Linear codes.

    This is the default encoder of a generic linear code, and should never be used for other codes than
    :class:`LinearCode`.

    INPUT:

    - ``code`` -- The associated :class:`LinearCode` of this encoder.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-based encoder for [7, 4] linear code over GF(2)
        """
        super(LinearCodeGeneratorMatrixEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeGeneratorMatrixEncoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: E1 = LinearCode(G).encoder()
            sage: E2 = LinearCode(G).encoder()
            sage: E1 == E2
            True
        """
        return isinstance(other, LinearCodeGeneratorMatrixEncoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-based encoder for [7, 4] linear code over GF(2)
        """
        return "Generator matrix-based encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: latex(E)
            \textnormal{Generator matrix-based encoder for }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Generator matrix-based encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.generator_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
        """
        g = self.code().generator_matrix()
        g.set_immutable()
        return g


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

    - ``code`` -- The associated code of this encoder.

    - ``systematic_positions`` -- (default: ``None``) the positions in codewords that
      should correspond to the message symbols. A list of `k` distinct integers in
      the range 0 to `n-1` where `n` is the length of the code and `k` its
      dimension. The 0th symbol of a message will then be at position
      ``systematic_positions[0]``, the 1st index at position
      ``systematic_positions[1]``, etc. A ``ValueError`` is raised at
      construction time if the supplied indices do not form an information set.

    EXAMPLES:

    The following demonstrates the basic usage of :class:`LinearCodeSystematicEncoder`::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0,0],\
                                     [1,0,0,1,1,0,0,0],\
                                     [0,1,0,1,0,1,0,0],\
                                     [1,1,0,1,0,0,1,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1 1]
            [0 1 0 0 1 0 1 1]
            [0 0 1 0 1 1 0 0]
            [0 0 0 1 1 1 1 1]
            sage: E2 = codes.encoders.LinearCodeSystematicEncoder(C, systematic_positions=[5,4,3,2])
            sage: E2.generator_matrix()
            [1 0 0 0 0 1 1 1]
            [0 1 0 0 1 0 1 1]
            [1 1 0 1 0 0 1 1]
            [1 1 1 0 0 0 0 0]

    An error is raised if one specifies systematic positions which do not form
    an information set::

            sage: E3 = codes.encoders.LinearCodeSystematicEncoder(C, systematic_positions=[0,1,6,7])
            Traceback (most recent call last):
            ...
            ValueError: systematic_positions are not an information set


    We exemplify how to use :class:`LinearCodeSystematicEncoder` as the default
    encoder. The following class is the dual of the repetition code::

        sage: class DualRepetitionCode(sage.coding.linear_code.AbstractLinearCode):
        ....:   def __init__(self, field, length):
        ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "Systematic", "Syndrome")
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
        ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self, field, length, "Systematic", "Syndrome")
        ....:
        ....:   def _repr_(self):
        ....:       return "I am a badly defined code"
        ....:
        sage: BadCodeFamily(GF(3), 5).generator_matrix()
        Traceback (most recent call last):
        ...
        ValueError: a parity check matrix must be specified if LinearCodeSystematicEncoder is the default encoder
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
        super(LinearCodeSystematicEncoder, self).__init__(code)
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
        Tests equality between LinearCodeSystematicEncoder objects.

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
        Returns a generator matrix in systematic form of the associated code of ``self``.

        Systematic form here means that a subsets of the columns of the matrix
        forms the identity matrix.

        .. NOTE::

            The matrix returned by this method will not necessarily be `[I \vert H]`, where `I`
            is the identity block and `H` the parity block. If one wants to know which columns
            create the identity block, one can call :meth:`systematic_positions`

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],\
                                     [1,0,0,1,1,0,0],\
                                     [0,1,0,1,0,1,0],\
                                     [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]

        We can ask for different systematic positions::

            sage: E2 = codes.encoders.LinearCodeSystematicEncoder(C, systematic_positions=[5,4,3,2])
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
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
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
        Returns a permutation which would take the systematic positions into [0,..,k-1]

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
        lp = [ None ]*n
        for (i,j) in zip(range(k), systematic_positions):
            lp[i] = j
        j = k
        set_sys_pos = set(systematic_positions)
        for i in range(n):
            if not i in set_sys_pos:
                lp[j] = i
                j += 1
        from sage.combinat.permutation import Permutation
        return Permutation([1 + e for e in lp])

    def systematic_positions(self):
        r"""
        Returns a tuple containing the indices of the columns which form an
        identity matrix when the generator matrix is in systematic form.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],\
                                     [1,0,0,1,1,0,0],\
                                     [0,1,0,1,0,1,0],\
                                     [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: E.systematic_positions()
            (0, 1, 2, 3)

        We take another matrix with a less nice shape::

            sage: G = Matrix(GF(2), [[1,1,0,0,1,0,1],\
                                     [1,1,0,0,1,0,0],\
                                     [0,0,1,0,0,1,0],\
                                     [0,0,1,0,1,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
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
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C, systematic_positions=[0,1,3])
            sage: E.systematic_positions()
            (0, 1, 3)
        """
        return self._systematic_positions if self._systematic_positions else self.generator_matrix().pivots()


####################### decoders ###############################

class LinearCodeSyndromeDecoder(Decoder):
    r"""
    Constructs a decoder for Linear Codes based on syndrome lookup table.

    The decoding algorithm works as follows:

    - First, a lookup table is built by computing the syndrome of every error
      pattern of weight up to ``maximum_error_weight``.
    - Then, whenever one tries to decode a word ``r``, the syndrome of ``r`` is
      computed. The corresponding error pattern is recovered from the
      pre-computed lookup table.
    - Finally, the recovered error pattern is subtracted from ``r`` to recover
      the original word.

    ``maximum_error_weight`` need never exceed the covering radius of the code,
    since there are then always lower-weight errors with the same syndrome. If
    one sets ``maximum_error_weight`` to a value greater than the covering
    radius, then the covering radius will be determined while building the
    lookup-table. This lower value is then returned if you query
    ``decoding_radius`` after construction.

    If ``maximum_error_weight`` is left unspecified or set to a number at least
    the covering radius of the code, this decoder is complete, i.e. it decodes
    every vector in the ambient space.

    .. NOTE::

        Constructing the lookup table takes time exponential in the length of the
        code and the size of the code's base field. Afterwards, the individual
        decodings are fast.

    INPUT:

    - ``code`` -- A code associated to this decoder

    - ``maximum_error_weight`` -- (default: ``None``) the maximum number of
      errors to look for when building the table. An error is raised if it is
      set greater than `n-k`, since this is an upper bound on the covering
      radius on any linear code. If ``maximum_error_weight`` is kept
      unspecified, it will be set to `n - k`, where `n` is the length of
      ``code`` and `k` its dimension.

    EXAMPLES::

        sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
        sage: C = LinearCode(G)
        sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
        sage: D
        Syndrome decoder for [9, 3] linear code over GF(3) handling errors of weight up to 4

    If one wants to correct up to a lower number of errors, one can do as follows::

        sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, maximum_error_weight=2)
        sage: D
        Syndrome decoder for [9, 3] linear code over GF(3) handling errors of weight up to 2

    If one checks the list of types of this decoder before constructing it,
    one will notice it contains the keyword ``dynamic``.
    Indeed, the behaviour of the syndrome decoder depends on the maximum
    error weight one wants to handle, and how it compares to the minimum
    distance and the covering radius of ``code``.
    In the following examples, we illustrate this property by computing
    different instances of syndrome decoder for the same code.

    We choose the following linear code, whose covering radius equals to 4
    and minimum distance to 5 (half the minimum distance is 2)::

        sage: G = matrix(GF(5), [[1, 0, 0, 0, 0, 4, 3, 0, 3, 1, 0],
        ....:                    [0, 1, 0, 0, 0, 3, 2, 2, 3, 2, 1],
        ....:                    [0, 0, 1, 0, 0, 1, 3, 0, 1, 4, 1],
        ....:                    [0, 0, 0, 1, 0, 3, 4, 2, 2, 3, 3],
        ....:                    [0, 0, 0, 0, 1, 4, 2, 3, 2, 2, 1]])
        sage: C = LinearCode(G)

    In the following examples, we illustrate how the choice of
    ``maximum_error_weight`` influences the types of the instance of
    syndrome decoder, alongside with its decoding radius.

    We build a first syndrome decoder, and pick a ``maximum_error_weight``
    smaller than both the covering radius and half the minimum distance::

        sage: D = C.decoder("Syndrome", maximum_error_weight = 1)
        sage: D.decoder_type()
        {'always-succeed', 'bounded_distance', 'hard-decision'}
        sage: D.decoding_radius()
        1

    In that case, we are sure the decoder will always succeed. It is also
    a bounded distance decoder.

    We now build another syndrome decoder, and this time,
    ``maximum_error_weight`` is chosen to be bigger than half the minimum distance,
    but lower than the covering radius::

        sage: D = C.decoder("Syndrome", maximum_error_weight = 3)
        sage: D.decoder_type()
        {'bounded_distance', 'hard-decision', 'might-error'}
        sage: D.decoding_radius()
        3

    Here, we still get a bounded distance decoder.
    But because we have a maximum error weight bigger than half the
    minimum distance, we know it might return a codeword which was not
    the original codeword.

    And now, we build a third syndrome decoder, whose ``maximum_error_weight``
    is bigger than both the covering radius and half the minimum distance::

        sage: D = C.decoder("Syndrome", maximum_error_weight = 5)
        sage: D.decoder_type()
        {'complete', 'hard-decision', 'might-error'}
        sage: D.decoding_radius()
        4

    In that case, the decoder might still return an unexpected codeword, but
    it is now complete. Note the decoding radius is equal to 4: it was
    determined while building the syndrome lookup table that any error with
    weight more than 4 will be decoded incorrectly. That is because the covering
    radius for the code is 4.

    The minimum distance and the covering radius are both determined while
    computing the syndrome lookup table. They user did not explicitly ask to
    compute these on the code ``C``. The dynamic typing of the syndrome decoder
    might therefore seem slightly surprising, but in the end is quite
    informative.
    """

    def __init__(self, code, maximum_error_weight=None):
        r"""
        TESTS:

        If ``maximum_error_weight`` is greater or equal than `n-k`, where `n`
        is ``code``'s length, and `k` is ``code``'s dimension,
        an error is raised::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, 42)
            Traceback (most recent call last):
            ...
            ValueError: maximum_error_weight has to be less than code's length minus its dimension

        The Syndrome Decoder of a Hamming code should have types
        ``minimum-distance`` and ``always-succeed`` (see :trac:`20898`)::

            sage: C = codes.HammingCode(GF(5), 3)
            sage: D = C.decoder("Syndrome")
            sage: C.minimum_distance()
            3
            sage: D.maximum_error_weight()
            1
            sage: D.decoder_type()
            {'always-succeed', 'complete', 'hard-decision', 'minimum-distance'}
        """
        n_minus_k = code.length() - code.dimension()
        if maximum_error_weight is None:
            self._maximum_error_weight = n_minus_k
        elif not isinstance(maximum_error_weight, (Integer, int)):
            raise ValueError("maximum_error_weight has to be a Sage integer or a Python int")
        elif maximum_error_weight > n_minus_k:
            raise ValueError("maximum_error_weight has to be less than code's length minus its dimension")
        else:
            self._maximum_error_weight = maximum_error_weight
        super(LinearCodeSyndromeDecoder, self).__init__(code, code.ambient_space(),\
                code._default_encoder_name)
        self._lookup_table = self._build_lookup_table()

    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeSyndromeDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: D1 = codes.decoders.LinearCodeSyndromeDecoder(LinearCode(G))
            sage: D2 = codes.decoders.LinearCodeSyndromeDecoder(LinearCode(G))
            sage: D1 == D2
            True
        """
        return (isinstance(other, LinearCodeSyndromeDecoder) and
                self.code() == other.code() and
                self.maximum_error_weight() == other.maximum_error_weight())

    def __hash__(self):
        """
        Return the hash of self.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: D1 = codes.decoders.LinearCodeSyndromeDecoder(LinearCode(G))
            sage: D2 = codes.decoders.LinearCodeSyndromeDecoder(LinearCode(G))
            sage: hash(D1) == hash(D2)
            True
        """
        return hash((self.code(), self.maximum_error_weight()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D
            Syndrome decoder for [9, 3] linear code over GF(3) handling errors of weight up to 4
        """
        return "Syndrome decoder for %s handling errors of weight up to %s" % (self.code(), self.maximum_error_weight())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: latex(D)
            \textnormal{Syndrome decoder for [9, 3]\textnormal{ Linear code over }\Bold{F}_{3} handling errors of weight up to 4}
        """
        return "\\textnormal{Syndrome decoder for %s handling errors of weight up to %s}" % (self.code()._latex_(), self.maximum_error_weight())

    @cached_method
    def _build_lookup_table(self):
        r"""
        Builds lookup table for all possible error patterns of weight up to :meth:`maximum_error_weight`.

        EXAMPLES::

            sage: G = Matrix(GF(3),[
            ....:   [1, 0, 0, 0, 2, 2, 1, 1],
            ....:   [0, 1, 0, 0, 0, 0, 1, 1],
            ....:   [0, 0, 1, 0, 2, 0, 0, 2],
            ....:   [0, 0, 0, 1, 0, 2, 0, 1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, maximum_error_weight = 1)
            sage: D._build_lookup_table()
            {(0, 0, 0, 0): (0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 1): (0, 0, 0, 0, 1, 0, 0, 0),
             (0, 0, 0, 2): (0, 0, 0, 0, 2, 0, 0, 0),
             (0, 0, 1, 0): (0, 0, 1, 0, 0, 0, 0, 0),
             (0, 0, 1, 2): (0, 0, 0, 0, 0, 0, 0, 1),
             (0, 0, 2, 0): (0, 0, 2, 0, 0, 0, 0, 0),
             (0, 0, 2, 1): (0, 0, 0, 0, 0, 0, 0, 2),
             (0, 1, 0, 0): (0, 1, 0, 0, 0, 0, 0, 0),
             (0, 1, 1, 2): (0, 0, 0, 0, 0, 0, 2, 0),
             (0, 2, 0, 0): (0, 2, 0, 0, 0, 0, 0, 0),
             (0, 2, 2, 1): (0, 0, 0, 0, 0, 0, 1, 0),
             (1, 0, 0, 0): (1, 0, 0, 0, 0, 0, 0, 0),
             (1, 2, 0, 2): (0, 0, 0, 0, 0, 1, 0, 0),
             (1, 2, 2, 0): (0, 0, 0, 1, 0, 0, 0, 0),
             (2, 0, 0, 0): (2, 0, 0, 0, 0, 0, 0, 0),
             (2, 1, 0, 1): (0, 0, 0, 0, 0, 2, 0, 0),
             (2, 1, 1, 0): (0, 0, 0, 2, 0, 0, 0, 0)}

        TESTS:

        Check that :trac:`24114` is fixed::

            sage: R.<x> = PolynomialRing(GF(3))
            sage: f = x^2 + x + 2
            sage: K.<a> = f.root_field()
            sage: H = Matrix(K,[[1,2,1],[2*a+1,a,1]])
            sage: C = codes.from_parity_check_matrix(H)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.syndrome_table()
             {(0, 0): (0, 0, 0),
              (0, 1): (0, 1, 0),
              (0, 2): (0, 2, 0),
              (0, a): (0, a, 0),
             ...
              (2*a + 2, 2*a): (0, 0, 2),
              (2*a + 2, 2*a + 1): (2*a + 2, 2*a + 1, 0),
              (2*a + 2, 2*a + 2): (2*a + 2, 2*a + 2, 0)}
        """
        t = self._maximum_error_weight
        self._code_covering_radius = None
        self._code_minimum_distance = None
        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        C = self.code()
        n = C.length()
        k = C.dimension()
        H = C.parity_check_matrix()
        F = C.base_ring()
        l = list(F)
        zero = F.zero()
        #Builds a list of generators of all error positions for all
        #possible error weights
        if zero in l:
            l.remove(zero)
        # Remember to include the no-error-vector to handle codes of minimum
        # distance 1 gracefully
        zero_syndrome = vector(F,[F.zero()]*(n-k))
        zero_syndrome.set_immutable()
        lookup = { zero_syndrome : vector(F,[F.zero()]*n) }
        error_position_tables = [cartesian_product([l]*i) for i in range(1, t+1)]
        first_collision = True
        #Filling the lookup table
        for i in range(1, t+1):
            stop = True
            patterns = Subsets(range(n), i)
            basic = vector(F, n)
            for p in patterns:
                for error in error_position_tables[i-1]:
                    ind = 0
                    e = copy(basic)
                    for pos in p:
                        e[pos] = error[ind]
                        ind += 1
                    s = H * e
                    s.set_immutable()
                    try:
                        e_cur = lookup[s]
                        #if this is the first time we see a collision
                        #we learn the minimum distance of the code
                        if first_collision:
                            self._code_minimum_distance = e.hamming_weight() + e_cur.hamming_weight()
                            first_collision = False
                    except KeyError:
                        stop = False
                        lookup[s] = copy(e)
            #if we reached the early termination condition
            #we learn the covering radius of the code
            if stop:
                self._code_covering_radius = i - 1
                self._maximum_error_weight = self._code_covering_radius
                break
        # Update decoder types depending on whether we are decoding up to covering radius
        if self._code_covering_radius:
            self._decoder_type.add("complete")
        else:
            self._decoder_type.add("bounded_distance")
        # Update decoder types depending on whether we are decoding beyond d/2
        if self._code_minimum_distance:
            if self._maximum_error_weight == (self._code_minimum_distance-1)//2:
                self._decoder_type.update({"minimum-distance","always-succeed"})
            else:
                # then t > (d-1)/2
                self._decoder_type.add("might-error")
        else:
            self._decoder_type.add("always-succeed")
        return lookup


    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self``'s message space

        EXAMPLES::

            sage: G = Matrix(GF(3),[
            ....:   [1, 0, 0, 0, 2, 2, 1, 1],
            ....:   [0, 1, 0, 0, 0, 0, 1, 1],
            ....:   [0, 0, 1, 0, 2, 0, 0, 2],
            ....:   [0, 0, 0, 1, 0, 2, 0, 1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, maximum_error_weight = 2)
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: c = C.random_element()
            sage: r = Chan(c)
            sage: c == D.decode_to_code(r)
            True
        """
        lookup_table = self.syndrome_table()
        s = self.code().parity_check_matrix() * r
        s.set_immutable()
        if s.is_zero():
            return r
        err = lookup_table[s]
        r_corr = copy(r)
        for i in range(self.code().length()):
            r_corr[i] = r[i] - err[i]
        return r_corr

    def maximum_error_weight(self):
        r"""
        Returns the maximal number of errors a received word can have
        and for which ``self`` is guaranteed to return a most likely codeword.

        Same as ``self.decoding_radius``.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.maximum_error_weight()
            4
        """
        return self._maximum_error_weight

    def decoding_radius(self):
        r"""
        Returns the maximal number of errors a received word can have
        and for which ``self`` is guaranteed to return a most likely codeword.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.decoding_radius()
            4
        """
        return self._maximum_error_weight

    def syndrome_table(self):
        r"""
        Return the syndrome lookup table of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.syndrome_table()
            {(0, 0, 0): (0, 0, 0, 0, 0, 0, 0),
             (0, 0, 1): (0, 0, 0, 1, 0, 0, 0),
             (0, 1, 0): (0, 1, 0, 0, 0, 0, 0),
             (0, 1, 1): (0, 0, 0, 0, 0, 1, 0),
             (1, 0, 0): (1, 0, 0, 0, 0, 0, 0),
             (1, 0, 1): (0, 0, 0, 0, 1, 0, 0),
             (1, 1, 0): (0, 0, 1, 0, 0, 0, 0),
             (1, 1, 1): (0, 0, 0, 0, 0, 0, 1)}
        """
        return self._lookup_table


class LinearCodeNearestNeighborDecoder(Decoder):
    r"""
    Construct a decoder for Linear Codes. This decoder will decode to the
    nearest codeword found.

    INPUT:

    - ``code`` -- A code associated to this decoder
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [7, 4] linear code over GF(2)
        """
        super(LinearCodeNearestNeighborDecoder, self).__init__(code, code.ambient_space(), \
                code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeNearestNeighborDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: D1 = codes.decoders.LinearCodeNearestNeighborDecoder(LinearCode(G))
            sage: D2 = codes.decoders.LinearCodeNearestNeighborDecoder(LinearCode(G))
            sage: D1 == D2
            True
        """
        return isinstance(other, LinearCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [7, 4] linear code over GF(2)
        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: latex(D)
            \textnormal{Nearest neighbor decoder for }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Nearest neighbor decoder for }%s" % self.code()._latex_()

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self``'s message space

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: D.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)
        """
        c_min = self.code().zero()
        h_min = r.hamming_weight()
        for c in self.code():
            if (c-r).hamming_weight() < h_min:
                h_min = (c-r).hamming_weight()
                c_min = c
        c_min.set_immutable()
        return c_min

    def decoding_radius(self):
        r"""
        Return maximal number of errors ``self`` can decode.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: D.decoding_radius()
            1
        """
        return (self.code().minimum_distance()-1) // 2


####################### registration ###############################

LinearCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder

LinearCodeSyndromeDecoder._decoder_type = {"hard-decision", "dynamic"}
LinearCodeNearestNeighborDecoder._decoder_type = {"hard-decision", "always-succeed", "complete"}
