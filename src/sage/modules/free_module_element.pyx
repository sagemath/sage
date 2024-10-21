# sage_setup: distribution = sagemath-modules
r"""
Elements of free modules

AUTHORS:

- William Stein

- Josh Kantor

- Thomas Feulner (2012-11): Added :meth:`FreeModuleElement.hamming_weight` and
  :meth:`FreeModuleElement_generic_sparse.hamming_weight`

- Jeroen Demeyer (2015-02-24): Implement fast Cython methods
  ``get_unsafe`` and ``set_unsafe`` similar to other places in Sage
  (:issue:`17562`)

EXAMPLES: We create a vector space over `\QQ` and a
subspace of this space.

::

    sage: V = QQ^5
    sage: W = V.span([V.1, V.2])

Arithmetic operations always return something in the ambient space,
since there is a canonical map from `W` to `V` but
not from `V` to `W`.

::

    sage: parent(W.0 + V.1)
    Vector space of dimension 5 over Rational Field
    sage: parent(V.1 + W.0)
    Vector space of dimension 5 over Rational Field
    sage: W.0 + V.1
    (0, 2, 0, 0, 0)
    sage: W.0 - V.0
    (-1, 1, 0, 0, 0)

Next we define modules over `\ZZ` and a finite
field.

::

    sage: K = ZZ^5
    sage: M = GF(7)^5

Arithmetic between the `\QQ` and
`\ZZ` modules is defined, and the result is always
over `\QQ`, since there is a canonical coercion map
to `\QQ`.

::

    sage: K.0 + V.1
    (1, 1, 0, 0, 0)
    sage: parent(K.0 + V.1)
    Vector space of dimension 5 over Rational Field

Since there is no canonical coercion map to the finite field from
`\QQ` the following arithmetic is not defined::

    sage: V.0 + M.0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +:
     'Vector space of dimension 5 over Rational Field' and
     'Vector space of dimension 5 over Finite Field of size 7'

However, there is a map from `\ZZ` to the finite
field, so the following is defined, and the result is in the finite
field.

::

    sage: w = K.0 + M.0; w
    (2, 0, 0, 0, 0)
    sage: parent(w)
    Vector space of dimension 5 over Finite Field of size 7
    sage: parent(M.0 + K.0)
    Vector space of dimension 5 over Finite Field of size 7

Matrix vector multiply::

    sage: MS = MatrixSpace(QQ, 3)
    sage: A = MS([0,1,0,1,0,0,0,0,1])
    sage: V = QQ^3
    sage: v = V([1,2,3])
    sage: v * A
    (2, 1, 3)

TESTS::

    sage: D = 46341
    sage: u = 7
    sage: R = Integers(D)
    sage: p = matrix(R,[[84, 97, 55, 58, 51]])
    sage: 2*p.row(0)                                                                    # needs sage.libs.pari
    (168, 194, 110, 116, 102)

This is a test from :issue:`20211`::

    sage: MatrixSpace(ZZ, 1, 1)(vector([1]))
    [1]
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport cython
from cpython.slice cimport PySlice_GetIndicesEx

from sage.categories.rings import Rings
from sage.structure.sequence import Sequence
from sage.structure.element cimport Element, RingElement, Vector
from sage.structure.element import canonical_coercion
from sage.structure.richcmp cimport richcmp_not_equal, richcmp, rich_to_bool

import sage.rings.abc
from sage.rings.infinity import Infinity, AnInfinity
from sage.rings.integer_ring import ZZ
from sage.rings.abc import RealDoubleField, ComplexDoubleField

from sage.rings.integer cimport Integer, smallInteger
from sage.arith.numerical_approx cimport digits_to_bits

# For the norm function, we cache Sage integers 1 and 2
__one__ = smallInteger(1)
__two__ = smallInteger(2)


def is_FreeModuleElement(x):
    """
    EXAMPLES::

        sage: sage.modules.free_module_element.is_FreeModuleElement(0)
        doctest:warning...
        DeprecationWarning: The function is_FreeModuleElement is deprecated;
        use 'isinstance(..., FreeModuleElement)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        False
        sage: sage.modules.free_module_element.is_FreeModuleElement(vector([1,2,3]))
        True
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38184,
                       "The function is_FreeModuleElement is deprecated; "
                       "use 'isinstance(..., FreeModuleElement)' instead.")
    return isinstance(x, FreeModuleElement)


def vector(arg0, arg1=None, arg2=None, sparse=None, immutable=False):
    r"""
    Return a vector or free module element with specified entries.

    CALL FORMATS:

    This constructor can be called in several different ways.
    In each case, ``sparse=True`` or ``sparse=False`` as well
    as ``immutable=True`` or ``immutable=False`` can be
    supplied as an option.  ``free_module_element()`` is an
    alias for ``vector()``.

        1. vector(object)

        2. vector(ring, object)

        3. vector(object, ring)

        4. vector(ring, degree, object)

        5. vector(ring, degree)

    INPUT:

    - ``object`` -- list, dictionary, or other
      iterable containing the entries of the vector, including
      any object that is palatable to the ``Sequence`` constructor

    - ``ring`` -- a base ring (or field) for the vector space or free module,
      which contains all of the elements

    - ``degree`` -- integer specifying the number of
      entries in the vector or free module element

    - ``sparse`` -- boolean, whether the result should be a sparse vector

    - ``immutable`` -- boolean (default: ``False``); whether the result
      should be an immutable vector

    In call format 4, an error is raised if the ``degree`` does not match
    the length of ``object`` so this call can provide some safeguards.
    Note however that using this format when ``object`` is a dictionary
    is unlikely to work properly.

    OUTPUT:

    An element of the ambient vector space or free module with the
    given base ring and implied or specified dimension or rank,
    containing the specified entries and with correct degree.

    In call format 5, no entries are specified, so the element is
    populated with all zeros.

    If the ``sparse`` option is not supplied, the output will
    generally have a dense representation.  The exception is if
    ``object`` is a dictionary, then the representation will be sparse.

    EXAMPLES::

        sage: v = vector([1,2,3]); v
        (1, 2, 3)
        sage: v.parent()
        Ambient free module of rank 3 over the principal ideal domain Integer Ring
        sage: v = vector([1,2,3/5]); v
        (1, 2, 3/5)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field

    All entries must *canonically* coerce to some common ring::

        sage: v = vector([17, GF(11)(5), 19/3]); v
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements

    ::

        sage: v = vector([17, GF(11)(5), 19]); v
        (6, 5, 8)
        sage: v.parent()
        Vector space of dimension 3 over Finite Field of size 11
        sage: v = vector([17, GF(11)(5), 19], QQ); v
        (17, 5, 19)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field
        sage: v = vector((1,2,3), QQ); v
        (1, 2, 3)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field
        sage: v = vector(QQ, (1,2,3)); v
        (1, 2, 3)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field
        sage: v = vector(vector([1,2,3])); v
        (1, 2, 3)
        sage: v.parent()
        Ambient free module of rank 3 over the principal ideal domain Integer Ring

    You can also use ``free_module_element``, which is
    the same as ``vector``. ::

        sage: free_module_element([1/3, -4/5])
        (1/3, -4/5)

    We make a vector mod 3 out of a vector over `\ZZ`. ::

        sage: vector(vector([1,2,3]), GF(3))
        (1, 2, 0)

    The degree of a vector may be specified::

        sage: vector(QQ, 4, [1,1/2,1/3,1/4])
        (1, 1/2, 1/3, 1/4)

    But it is an error if the degree and size of the list of entries
    are mismatched::

        sage: vector(QQ, 5, [1,1/2,1/3,1/4])
        Traceback (most recent call last):
        ...
        ValueError: incompatible degrees in vector constructor

    Providing no entries populates the vector with zeros, but of course,
    you must specify the degree since it is not implied.  Here we use a
    finite field as the base ring. ::

        sage: w = vector(FiniteField(7), 4); w
        (0, 0, 0, 0)
        sage: w.parent()
        Vector space of dimension 4 over Finite Field of size 7

    The fastest method to construct a zero vector is to call the
    :meth:`~sage.modules.free_module.FreeModule_generic.zero_vector`
    method directly on a free module or vector space, since
    vector(...)  must do a small amount of type checking.  Almost as
    fast as the ``zero_vector()`` method is the
    :func:`~sage.modules.free_module_element.zero_vector` constructor,
    which defaults to the integers.  ::

        sage: vector(ZZ, 5)          # works fine
        (0, 0, 0, 0, 0)
        sage: (ZZ^5).zero_vector()   # very tiny bit faster
        (0, 0, 0, 0, 0)
        sage: zero_vector(ZZ, 5)     # similar speed to vector(...)
        (0, 0, 0, 0, 0)
        sage: z = zero_vector(5); z
        (0, 0, 0, 0, 0)
        sage: z.parent()
        Ambient free module of rank 5 over
        the principal ideal domain Integer Ring

    Here we illustrate the creation of sparse vectors by using a
    dictionary::

        sage: vector({1:1.1, 3:3.14})
        (0.000000000000000, 1.10000000000000, 0.000000000000000, 3.14000000000000)

    With no degree given, a dictionary of entries implicitly declares a
    degree by the largest index (key) present.  So you can provide a
    terminal element (perhaps a zero?) to set the degree.  But it is probably safer
    to just include a degree in your construction.  ::

        sage: v = vector(QQ, {0:1/2, 4:-6, 7:0}); v
        (1/2, 0, 0, 0, -6, 0, 0, 0)
        sage: v.degree()
        8
        sage: v.is_sparse()
        True
        sage: w = vector(QQ, 8, {0:1/2, 4:-6})
        sage: w == v
        True

    It is an error to specify a negative degree. ::

        sage: vector(RR, -4, [1.0, 2.0, 3.0, 4.0])
        Traceback (most recent call last):
        ...
        ValueError: cannot specify the degree of a vector as a negative integer (-4)

    It is an error to create a zero vector but not provide
    a ring as the first argument.  ::

        sage: vector('junk', 20)
        Traceback (most recent call last):
        ...
        TypeError: first argument must be base ring of zero vector, not junk

    And it is an error to specify an index in a dictionary
    that is greater than or equal to a requested degree. ::

        sage: vector(ZZ, 10, {3:4, 7:-2, 10:637})
        Traceback (most recent call last):
        ...
        ValueError: dictionary of entries has a key (index) exceeding the requested degree

    A 1-dimensional numpy array of type float or complex may be
    passed to vector. Unless an explicit ring is given, the result will
    be a vector in the appropriate dimensional vector space over the
    real double field or the complex double field. The data in the array
    must be contiguous, so column-wise slices of numpy matrices will
    raise an exception. ::

        sage: # needs numpy
        sage: import numpy
        sage: x = numpy.random.randn(10)
        sage: y = vector(x)
        sage: parent(y)
        Vector space of dimension 10 over Real Double Field
        sage: parent(vector(RDF, x))
        Vector space of dimension 10 over Real Double Field
        sage: parent(vector(CDF, x))
        Vector space of dimension 10 over Complex Double Field
        sage: parent(vector(RR, x))
        Vector space of dimension 10 over Real Field with 53 bits of precision
        sage: v = numpy.random.randn(10) * complex(0,1)
        sage: w = vector(v)
        sage: parent(w)
        Vector space of dimension 10 over Complex Double Field

    Multi-dimensional arrays are not supported::

        sage: # needs numpy
        sage: import numpy as np
        sage: a = np.array([[1, 2, 3], [4, 5, 6]], np.float64)
        sage: vector(a)
        Traceback (most recent call last):
        ...
        TypeError: cannot convert 2-dimensional array to a vector

    If any of the arguments to vector have Python type int, real,
    or complex, they will first be coerced to the appropriate Sage
    objects. This fixes :issue:`3847`. ::

        sage: v = vector([int(0)]); v
        (0)
        sage: v[0].parent()
        Integer Ring
        sage: v = vector(range(10)); v
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        sage: v[3].parent()
        Integer Ring
        sage: v = vector([float(23.4), int(2), complex(2+7*I), 1]); v                   # needs sage.symbolic
        (23.4, 2.0, 2.0 + 7.0*I, 1.0)
        sage: v[1].parent()                                                             # needs sage.symbolic
        Complex Double Field

    If the argument is a vector, it doesn't change the base ring. This
    fixes :issue:`6643`::

        sage: # needs sage.rings.number_field
        sage: K.<sqrt3> = QuadraticField(3)
        sage: u = vector(K, (1/2, sqrt3/2))
        sage: vector(u).base_ring()
        Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?
        sage: v = vector(K, (0, 1))
        sage: vector(v).base_ring()
        Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?

    Constructing a vector from a numpy array behaves as expected::

        sage: # needs numpy
        sage: import numpy
        sage: a = numpy.array([1,2,3])
        sage: v = vector(a); v
        (1, 2, 3)
        sage: parent(v)
        Ambient free module of rank 3 over the principal ideal domain Integer Ring

    Complex numbers can be converted naturally to a sequence of length 2.  And
    then to a vector.  ::

        sage: c = CDF(2 + 3*I)                                                          # needs sage.rings.complex_double sage.symbolic
        sage: v = vector(c); v                                                          # needs sage.rings.complex_double sage.symbolic
        (2.0, 3.0)

    A generator, or other iterable, may also be supplied as input.  Anything
    that can be converted to a :class:`~sage.structure.sequence.Sequence` is
    a possible input.  ::

        sage: type(i^2 for i in range(3))
        <... 'generator'>
        sage: v = vector(i^2 for i in range(3)); v
        (0, 1, 4)

    An empty list, without a ring given, will default to the integers. ::

        sage: x = vector([]); x
        ()
        sage: x.parent()
        Ambient free module of rank 0 over the principal ideal domain Integer Ring

    The ``immutable`` switch allows to create an immutable vector. ::

        sage: v = vector(QQ, {0:1/2, 4:-6, 7:0}, immutable=True); v
        (1/2, 0, 0, 0, -6, 0, 0, 0)
        sage: v.is_immutable()
        True

    The ``immutable`` switch works regardless of the type of valid input to the
    constructor. ::

        sage: v = vector(ZZ, 4, immutable=True)
        sage: v.is_immutable()
        True
        sage: w = vector(ZZ, [1,2,3])
        sage: v = vector(w, ZZ, immutable=True)
        sage: v.is_immutable()
        True
        sage: v = vector(QQ, w, immutable=True)
        sage: v.is_immutable()
        True

        sage: # needs numpy sage.symbolic
        sage: import numpy as np
        sage: w = np.array([1, 2, pi], float)
        sage: v = vector(w, immutable=True)
        sage: v.is_immutable()
        True
        sage: w = np.array([i, 2, 3], complex)
        sage: v = vector(w, immutable=True)
        sage: v.is_immutable()
        True

    TESTS:

    We check that :issue:`31470` is fixed::

        sage: k.<a> = GF(5^3)                                                           # needs sage.rings.finite_rings
        sage: S.<x> = k['x', k.frobenius_endomorphism()]                                # needs sage.rings.finite_rings
        sage: vector(S, 3)                                                              # needs sage.rings.finite_rings
        ...
        (0, 0, 0)
    """
    from sage.modules.free_module import FreeModule
    # We first efficiently handle the important special case of the zero vector
    # over a ring. See trac 11657.
    # !! PLEASE DO NOT MOVE THIS CODE LOWER IN THIS FUNCTION !!
    arg1_integer = isinstance(arg1, (int, Integer))
    if arg2 is None and arg1_integer and arg0 in Rings():
        v = FreeModule(arg0, arg1, bool(sparse)).zero_vector()
        if immutable:
            v.set_immutable()
        return v

    # The try...except is slightly faster than testing with hasattr first
    # quite a significant amount of time.
    try:
        arg0_vector_ = arg0._vector_
    except AttributeError:
        pass
    else:
        v = arg0_vector_(arg1)
        if immutable:
            v.set_immutable()
        return v

    try:
        arg1_vector_ = arg1._vector_
    except AttributeError:
        pass
    else:
        v = arg1_vector_(arg0)
        if immutable:
            v.set_immutable()
        return v

    # consider a possible degree specified in second argument
    degree = None
    maxindex = None
    if arg1_integer:
        if arg1 < 0:
            raise ValueError("cannot specify the degree of a vector as a negative integer (%s)" % arg1)
        if isinstance(arg2, dict):
            maxindex = max([-1]+[index for index in arg2])
            if not maxindex < arg1:
                raise ValueError("dictionary of entries has a key (index) exceeding the requested degree")
        # arg1 is now a legitimate degree
        # With no arg2, we can try to return a zero vector
        #   else we size-check arg2 and slide it into arg1
        degree = arg1
        if arg2 is None:
            if arg0 not in Rings():
                msg = "first argument must be base ring of zero vector, not {0}"
                raise TypeError(msg.format(arg0))
        else:
            if not isinstance(arg2, dict) and len(arg2) != degree:
                raise ValueError("incompatible degrees in vector constructor")
            arg1 = arg2

    # Analyze arg0 and arg1 to create a ring (R) and entries (v)
    if arg0 in Rings():
        R = arg0
        v = arg1
    elif arg1 in Rings():
        R = arg1
        v = arg0
    else:
        v = arg0
        R = None

    try:
        import numpy
        from numpy import ndarray
        if int(numpy.version.short_version[0]) > 1:
            numpy.set_printoptions(legacy="1.25")

    except ImportError:
        pass
    else:
        if isinstance(v, ndarray):
            if len(v.shape) != 1:
                raise TypeError("cannot convert %r-dimensional array to a vector" % len(v.shape))
            from sage.modules.free_module import VectorSpace
            if (R is None or isinstance(R, RealDoubleField)) and v.dtype.kind == 'f':
                from sage.rings.real_double import RDF
                V = VectorSpace(RDF, v.shape[0])
                from sage.modules.vector_real_double_dense import Vector_real_double_dense
                v = Vector_real_double_dense(V, v)
                if immutable:
                    v.set_immutable()
                return v
            if (R is None or isinstance(R, ComplexDoubleField)) and v.dtype.kind == 'c':
                from sage.rings.complex_double import CDF
                V = VectorSpace(CDF, v.shape[0])
                from sage.modules.vector_complex_double_dense import Vector_complex_double_dense
                v = Vector_complex_double_dense(V, v)
                if immutable:
                    v.set_immutable()
                return v
            # Use slower conversion via list
            v = list(v)

    if isinstance(v, dict):
        if degree is None:
            degree = max([-1]+[index for index in v])+1
        if sparse is None:
            sparse = True
    else:
        degree = None
        if sparse is None:
            sparse = False

    v, R = prepare(v, R, degree)

    M = FreeModule(R, len(v), bool(sparse))
    w = M(v)
    if immutable:
        w.set_immutable()
    return w


free_module_element = vector


def prepare(v, R, degree=None):
    r"""
    Convert an object describing elements of a vector into a list of entries in a common ring.

    INPUT:

    - ``v`` -- dictionary with nonnegative integers as keys,
      or a list or other object that can be converted by the ``Sequence``
      constructor
    - ``R`` -- a ring containing all the entries, possibly given as ``None``
    - ``degree`` -- a requested size for the list when the input is a dictionary,
      otherwise ignored

    OUTPUT: a pair

    The first item is a list of the values specified in the object ``v``.
    If the object is a dictionary , entries are placed in the list
    according to the indices that were their keys in the dictionary,
    and the remainder of the entries are zero.  The value of
    ``degree`` is assumed to be larger than any index provided
    in the dictionary and will be used as the number of entries
    in the returned list.

    The second item returned is a ring that contains all of
    the entries in the list. If ``R`` is given, the entries
    are coerced in.  Otherwise a common ring is found. For
    more details, see the
    :class:`~sage.structure.sequence.Sequence` object.  When ``v``
    has no elements and ``R`` is ``None``, the ring returned is
    the integers.

    EXAMPLES::

        sage: from sage.modules.free_module_element import prepare
        sage: prepare([1, 2/3, 5], None)
        ([1, 2/3, 5], Rational Field)

        sage: prepare([1, 2/3, 5], RR)
        ([1.00000000000000, 0.666666666666667, 5.00000000000000],
         Real Field with 53 bits of precision)

        sage: prepare({1: 4, 3: -2}, ZZ, 6)
        ([0, 4, 0, -2, 0, 0], Integer Ring)

        sage: prepare({3: 1, 5: 3}, QQ, 6)
        ([0, 0, 0, 1, 0, 3], Rational Field)

        sage: prepare([1, 2/3, '10', 5], RR)
        ([1.00000000000000, 0.666666666666667, 10.0000000000000, 5.00000000000000],
         Real Field with 53 bits of precision)

        sage: prepare({}, QQ, 0)
        ([], Rational Field)

        sage: prepare([1, 2/3, '10', 5], None)
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements

    Some objects can be converted to sequences even if they are not always
    thought of as sequences.  ::

        sage: c = CDF(2 + 3*I)                                                          # needs sage.symbolic
        sage: prepare(c, None)                                                          # needs sage.symbolic
        ([2.0, 3.0], Real Double Field)

    This checks a bug listed at :issue:`10595`. Without good evidence
    for a ring, the default is the integers. ::

        sage: prepare([], None)
        ([], Integer Ring)
    """
    if isinstance(v, dict):
        # convert to a list
        X = [0]*degree
        for key, value in v.iteritems():
            X[key] = value
        v = X
    # convert to a Sequence over common ring
    # default to ZZ on an empty list
    if R is None:
        try:
            if len(v) == 0:
                R = ZZ
        except TypeError:
            pass
    v = Sequence(v, universe=R, use_sage_types=True)
    ring = v.universe()
    if ring not in Rings():
        raise TypeError("unable to find a common ring for all elements")
    return v, ring


def zero_vector(arg0, arg1=None):
    r"""
    Return a vector or free module element with a specified number of zeros.

    CALL FORMATS:

        1. zero_vector(degree)

        2. zero_vector(ring, degree)

    INPUT:

    - ``degree`` -- the number of zero entries in the vector or
      free module element

    - ``ring`` -- (default: ``ZZ``) the base ring of the vector
      space or module containing the constructed zero vector

    OUTPUT:

    A vector or free module element with ``degree`` entries,
    all equal to zero and belonging to the ring if specified.
    If no ring is given, a free module element over ``ZZ``
    is returned.

    EXAMPLES:

    A zero vector over the field of rationals. ::

        sage: v = zero_vector(QQ, 5); v
        (0, 0, 0, 0, 0)
        sage: v.parent()
        Vector space of dimension 5 over Rational Field

    A free module zero element. ::

        sage: w = zero_vector(Integers(6), 3); w
        (0, 0, 0)
        sage: w.parent()
        Ambient free module of rank 3 over Ring of integers modulo 6

    If no ring is given, the integers are used. ::

        sage: u = zero_vector(9); u
        (0, 0, 0, 0, 0, 0, 0, 0, 0)
        sage: u.parent()
        Ambient free module of rank 9 over the principal ideal domain Integer Ring

    Non-integer degrees produce an error. ::

        sage: zero_vector(5.6)
        Traceback (most recent call last):
        ...
        TypeError: Attempt to coerce non-integral RealNumber to Integer

    Negative degrees also give an error. ::

        sage: zero_vector(-3)
        Traceback (most recent call last):
        ...
        ValueError: rank (=-3) must be nonnegative

    Garbage instead of a ring will be recognized as such. ::

        sage: zero_vector(x^2, 5)                                                       # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: first argument must be a ring
    """
    if arg1 is None:
        arg0 = ZZ(arg0)
        # default to a zero vector over the integers (ZZ) if no ring given
        return (ZZ**arg0).zero_vector()
    if arg0 in Rings():
        return (arg0**arg1).zero_vector()
    raise TypeError("first argument must be a ring")


def random_vector(ring, degree=None, *args, **kwds):
    r"""
    Return a vector (or module element) with random entries.

    INPUT:

    - ``ring`` -- (default: ``ZZ``) the base ring for the entries
    - ``degree`` -- nonnegative integer for the number of entries in the vector
    - ``sparse`` -- (default: ``False``) whether to use a sparse implementation
    - ``*args``, ``**kwds`` -- additional arguments and keywords are passed
      to the ``random_element()`` method of the ring

    OUTPUT:

    A vector, or free module element, with ``degree`` elements
    from ``ring``, chosen randomly from the ring according to
    the ring's ``random_element()`` method.

    .. NOTE::
        See below for examples of how random elements are
        generated by some common base rings.

    EXAMPLES:

    First, module elements over the integers.
    The default distribution is tightly clustered around -1, 0, 1.
    Uniform distributions can be specified by giving bounds, though
    the upper bound is never met.  See
    :meth:`sage.rings.integer_ring.IntegerRing_class.random_element`
    for several other variants. ::

        sage: random_vector(10).parent()
        Ambient free module of rank 10 over the principal ideal domain Integer Ring
        sage: random_vector(20).parent()
        Ambient free module of rank 20 over the principal ideal domain Integer Ring

        sage: v = random_vector(ZZ, 20, x=4)
        sage: all(i in range(4) for i in v)
        True

        sage: v = random_vector(ZZ, 20, x=-20, y=100)
        sage: all(i in range(-20, 100) for i in v)
        True

    If the ring is not specified, the default is the integers, and
    parameters for the random distribution may be passed without using
    keywords.  This is a random vector with 20 entries uniformly distributed
    between -20 and 100.  ::

        sage: random_vector(20, -20, 100).parent()
        Ambient free module of rank 20 over the principal ideal domain Integer Ring

    Now over the rationals.  Note that bounds on the numerator and
    denominator may be specified.  See
    :meth:`sage.rings.rational_field.RationalField.random_element`
    for documentation. ::

        sage: random_vector(QQ, 10).parent()
        Vector space of dimension 10 over Rational Field

        sage: v = random_vector(QQ, 10, num_bound=15, den_bound=5)
        sage: v.parent()
        Vector space of dimension 10 over Rational Field
        sage: all(q.numerator() <= 15 and q.denominator() <= 5 for q in v)
        True

    Inexact rings may be used as well.  The reals have
    uniform distributions, with the range `(-1,1)` as
    the default.  More at:
    :meth:`sage.rings.real_mpfr.RealField_class.random_element` ::

        sage: v = random_vector(RR, 5)
        sage: v.parent()
        Vector space of dimension 5 over Real Field with 53 bits of precision
        sage: all(-1 <= r <= 1 for r in v)
        True

        sage: v = random_vector(RR, 5, min=8, max=14)
        sage: v.parent()
        Vector space of dimension 5 over Real Field with 53 bits of precision
        sage: all(8 <= r <= 14 for r in v)
        True

    Any ring with a ``random_element()`` method may be used. ::

        sage: F = FiniteField(23)
        sage: hasattr(F, 'random_element')
        True
        sage: v = random_vector(F, 10)
        sage: v.parent()
        Vector space of dimension 10 over Finite Field of size 23

    The default implementation is a dense representation, equivalent to
    setting ``sparse=False``. ::

        sage: v = random_vector(10)
        sage: v.is_sparse()
        False

        sage: w = random_vector(ZZ, 20, sparse=True)
        sage: w.is_sparse()
        True

    The elements are chosen using the ring's ``random_element`` method::

        sage: from sage.misc.randstate import current_randstate
        sage: seed = current_randstate().seed()
        sage: set_random_seed(seed)
        sage: v1 = random_vector(ZZ, 20, distribution='1/n')
        sage: v2 = random_vector(ZZ, 15, x=-1000, y=1000)
        sage: v3 = random_vector(QQ, 10)
        sage: v4 = random_vector(FiniteField(17), 10)
        sage: v5 = random_vector(RR, 10)
        sage: set_random_seed(seed)
        sage: w1 = vector(ZZ.random_element(distribution='1/n') for _ in range(20))
        sage: w2 = vector(ZZ.random_element(x=-1000, y=1000) for _ in range(15))
        sage: w3 = vector(QQ.random_element() for _ in range(10))
        sage: [v1, v2, v3] == [w1, w2, w3]
        True
        sage: w4 = vector(FiniteField(17).random_element() for _ in range(10))
        sage: v4 == w4
        True
        sage: w5 = vector(RR.random_element() for _ in range(10))
        sage: v5 == w5
        True

    Inputs get checked before constructing the vector. ::

        sage: random_vector('junk')
        Traceback (most recent call last):
        ...
        TypeError: degree of a random vector must be an integer, not None

        sage: random_vector('stuff', 5)
        Traceback (most recent call last):
        ...
        TypeError: elements of a vector, or module element, must come from a ring, not stuff

        sage: random_vector(ZZ, -9)
        Traceback (most recent call last):
        ...
        ValueError: degree of a random vector must be nonnegative, not -9
    """
    if isinstance(ring, (Integer, int)):
        if degree is not None:
            arglist = list(args)
            arglist.insert(0, degree)
            args = tuple(arglist)
        degree = ring
        ring = ZZ
    if not isinstance(degree, (Integer, int)):
        raise TypeError("degree of a random vector must be an integer, not %s" % degree)
    if degree < 0:
        raise ValueError("degree of a random vector must be nonnegative, not %s" % degree)
    if ring not in Rings():
        raise TypeError("elements of a vector, or module element, must come from a ring, not %s" % ring)
    if not hasattr(ring, "random_element"):
        raise AttributeError("cannot create a random vector since there is no random_element() method for %s" % ring )
    sparse = kwds.pop('sparse', False)
    entries = [ring.random_element(*args, **kwds) for _ in range(degree)]
    return vector(ring, degree, entries, sparse)


cdef class FreeModuleElement(Vector):   # abstract base class
    """
    An element of a generic free module.
    """
    def __init__(self, parent):
        """
        EXAMPLES::

            sage: v = sage.modules.free_module_element.FreeModuleElement(QQ^3)
            sage: type(v)
            <class 'sage.modules.free_module_element.FreeModuleElement'>
        """
        self._parent = parent
        self._degree = parent.degree()
        self._is_immutable = 0

    # specified in ModulesWithBasis.ElementMethods.monomial_coefficients
    def monomial_coefficients(self, copy=True):
        r"""
        Return a dictionary whose keys are indices of basis elements
        in the support of ``self`` and whose values are the
        corresponding coefficients.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``self`` is internally
          represented by a dictionary ``d``, then make a copy of ``d``.
          If ``False``, then this can cause undesired behavior by
          mutating ``d``.

        EXAMPLES::

            sage: V = ZZ^3
            sage: v = V([1, 0, 5])
            sage: v.monomial_coefficients()
            {0: 1, 2: 5}

        Check that it works for submodules (:issue:`34455`)::

            sage: V = ZZ^3
            sage: U = V.submodule([[1, 2, 3], [1, 1, 1]])
            sage: U
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: u = U([2, 3, 4])
            sage: u.monomial_coefficients()
            {0: 2, 1: 3}
        """
        base_ring = self.parent().base_ring()
        if self.parent().is_ambient() and base_ring == self.parent().coordinate_ring():
            return self.dict(copy=copy)
        coordinates = self.parent().coordinate_vector(self)
        # coordinate_vector returns coefficients in the fraction field.
        # convert back to the base ring.
        return {index: base_ring(value)
                for index, value in enumerate(coordinates)
                if value}

    def _giac_init_(self):
        """
        EXAMPLES::

            sage: v = vector(ZZ, 4, range(4))
            sage: giac(v) + v                                                           # needs sage.libs.giac
            [0,2,4,6]

        ::

            sage: v = vector(QQ, 3, [2/3, 0, 5/4])
            sage: giac(v)                                                               # needs sage.libs.giac
            [2/3,0,5/4]

        ::

            sage: P.<x> = ZZ[]
            sage: v = vector(P, 3, [x^2 + 2, 2*x + 1, -2*x^2 + 4*x])
            sage: giac(v)                                                               # needs sage.libs.giac
            [sageVARx^2+2,2*sageVARx+1,-2*sageVARx^2+4*sageVARx]
        """
        return self.list()

    def __pari__(self):
        """
        Convert ``self`` to a PARI vector.

        OUTPUT: a PARI ``gen`` of type ``t_VEC``

        EXAMPLES::

            sage: v = vector(range(4))
            sage: v.__pari__()                                                          # needs sage.libs.pari
            [0, 1, 2, 3]
            sage: v.__pari__().type()                                                   # needs sage.libs.pari
            't_VEC'

        A list of vectors::

            sage: L = [vector(i^n for i in range(4)) for n in [1,3,5]]
            sage: pari(L)                                                               # needs sage.libs.pari
            [[0, 1, 2, 3], [0, 1, 8, 27], [0, 1, 32, 243]]
        """
        from sage.libs.pari.all import pari
        return pari(self.list())

    def _pari_init_(self):
        """
        Give a string which, when evaluated in GP, gives a PARI
        representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: v = vector(range(4))
            sage: v._pari_init_()                                                       # needs sage.libs.pari
            '[0,1,2,3]'

        Create the multiplication table of `GF(4)` using GP::

            sage: # needs sage.libs.pari
            sage: k.<a> = GF(4, impl='pari_ffelt')
            sage: v = gp(vector(list(k)))
            sage: v
            [0, 1, a, a + 1]
            sage: v.mattranspose() * v
            [0, 0, 0, 0; 0, 1, a, a + 1; 0, a, a + 1, 1; 0, a + 1, 1, a]
        """
        # Elements in vectors are always Sage Elements, so they should
        # have a _pari_init_() method.
        L = [x._pari_init_() for x in self.list()]
        return "[" + ",".join(L) + "]"

    def _magma_init_(self, magma):
        r"""
        Convert ``self`` to Magma.

        EXAMPLES::

            sage: F = FreeModule(ZZ, 2, inner_product_matrix=matrix(ZZ, 2, 2, [1, 0, 0, -1]))
            sage: v = F([1, 2])
            sage: M = magma(v); M                       # optional - magma
            (1 2)
            sage: M.Type()                              # optional - magma
            ModTupRngElt
            sage: M.Parent()                            # optional - magma
            Full RSpace of degree 2 over Integer Ring
            Inner Product Matrix:
            [ 1  0]
            [ 0 -1]
            sage: M.sage()                              # optional - magma
            (1, 2)
            sage: M.sage() == v                         # optional - magma
            True
            sage: M.sage().parent() is v.parent()       # optional - magma
            True

        ::

            sage: v = vector(QQ, [1, 2, 5/6])
            sage: M = magma(v); M                       # optional - magma
            (  1   2 5/6)
            sage: M.Type()                              # optional - magma
            ModTupFldElt
            sage: M.Parent()                            # optional - magma
            Full Vector space of degree 3 over Rational Field
            sage: M.sage()                              # optional - magma
            (1, 2, 5/6)
            sage: M.sage() == v                         # optional - magma
            True
            sage: M.sage().parent() is v.parent()       # optional - magma
            True
        """
        # Get a reference to Magma version of parent.
        R = magma(self.parent())
        # Get list of coefficients.
        v = ','.join(a._magma_init_(magma) for a in self.list())
        return '%s![%s]' % (R.name(), v)

    def numpy(self, dtype=object):
        """
        Convert ``self`` to a numpy array.

        INPUT:

        - ``dtype`` -- the `numpy dtype <http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_
          of the returned array

        EXAMPLES::

            sage: # needs numpy
            sage: v = vector([1,2,3])
            sage: v.numpy()
            array([1, 2, 3], dtype=object)
            sage: v.numpy() * v.numpy()
            array([1, 4, 9], dtype=object)

            sage: vector(QQ, [1, 2, 5/6]).numpy()                                       # needs numpy
            array([1, 2, 5/6], dtype=object)

        By default, the ``object`` `dtype <http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_ is used.
        Alternatively, the desired dtype can be passed in as a parameter::

            sage: # needs numpy
            sage: v = vector(QQ, [1, 2, 5/6])
            sage: v.numpy()
            array([1, 2, 5/6], dtype=object)
            sage: v.numpy(dtype=float)
            array([1.        , 2.        , 0.83333333])
            sage: v.numpy(dtype=int)
            array([1, 2, 0])
            sage: import numpy
            sage: v.numpy(dtype=numpy.uint8)
            array([1, 2, 0], dtype=uint8)

        Passing a dtype of None will let numpy choose a native type, which can
        be more efficient but may have unintended consequences::

            sage: # needs numpy
            sage: v.numpy(dtype=None)
            array([1.        , 2.        , 0.83333333])

            sage: w = vector(ZZ, [0, 1, 2^63 -1]); w
            (0, 1, 9223372036854775807)
            sage: wn = w.numpy(dtype=None); wn                                          # needs numpy
            array([                  0,                   1, 9223372036854775807]...)
            sage: wn.dtype                                                              # needs numpy
            dtype('int64')
            sage: w.dot_product(w)
            85070591730234615847396907784232501250
            sage: wn.dot(wn)        # overflow                                          # needs numpy
            2

        Numpy can give rather obscure errors; we wrap these to give a bit of context::

            sage: vector([1, 1/2, QQ['x'].0]).numpy(dtype=float)                        # needs numpy
            Traceback (most recent call last):
            ...
            ValueError: Could not convert vector over Univariate Polynomial Ring in x
            over Rational Field to numpy array of type <... 'float'>:
            setting an array element with a sequence.
        """
        import numpy
        from numpy import array
        if int(numpy.version.short_version[0]) > 1:
            numpy.set_printoptions(legacy="1.25")

        try:
            return array(self, dtype=dtype)
        except ValueError as e:
            raise ValueError(
                "Could not convert vector over %s to numpy array of type %s: %s" % (self.coordinate_ring(), dtype, e))

    def __hash__(self):
        """
        Return hash of this vector.  Only immutable vectors are hashable.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: v = vector([1,2/3,pi])
            sage: v.__hash__()
            Traceback (most recent call last):
            ...
            TypeError: mutable vectors are unhashable
            sage: v.set_immutable()
            sage: v.__hash__()   # random output
        """
        if not self._is_immutable:
            raise TypeError("mutable vectors are unhashable")
        return hash(tuple(self))

    def _vector_(self, R=None):
        r"""
        Return ``self`` as a vector.

        EXAMPLES::

            sage: v = vector(ZZ, [2, 12, 22])
            sage: vector(v)
            (2, 12, 22)
            sage: vector(GF(7), v)
            (2, 5, 1)
            sage: vector(v, ZZ['x', 'y'])
            (2, 12, 22)

            sage: vector(vector((1, 6.8)))
            (1.00000000000000, 6.80000000000000)
            sage: vector(vector(SR, (1, sqrt(2)) ) )                                    # needs sage.symbolic
            (1, sqrt(2))
        """
        if R is None:
            return self
        return self.change_ring(R)

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(vector(RR, [pi, e, 0.5]), verify=True)                     # needs sage.symbolic
            # Verified
            vector(RR, [3.1415926535897931, 2.7182818284590451, 0.5])
            sage: sage_input(vector(GF(5), [1, 2, 3, 4, 5]), verify=True)
            # Verified
            vector(GF(5), [1, 2, 3, 4, 0])
            sage: sage_input(vector([0, 0, 0, 1, 0, 0, 0], sparse=True), verify=True)
            # Verified
            vector(ZZ, {3:1, 6:0})
            sage: sage_input(vector(ZZ, []), verify=True)
            # Verified
            vector(ZZ, [])
            sage: sage_input(vector(RealField(27), [], sparse=True), verify=True)
            # Verified
            vector(RealField(27), {})
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: vector(ZZ, [42, 389])._sage_input_(SageInputBuilder(), False)
            {call: {atomic:vector}({atomic:ZZ}, {list: ({atomic:42}, {atomic:389})})}
            sage: vector(RDF, {1:pi, 1000:e})._sage_input_(SageInputBuilder(), False)   # needs sage.symbolic
            {call: {atomic:vector}({atomic:RDF}, {dict: {{atomic:1}:{atomic:3.1415926535897931}, {atomic:1000}:{atomic:2.718281828459045...}}})}
        """
        # Not a lot of room for prettiness here.
        # We always specify the ring, because that lets us use coerced=2
        # on the elements, which is sometimes much prettier than
        # the coerced=False we would get otherwise.
        if self.is_sparse_c():
            items = [(n, sib(e, 2))
                     for n,e in self.dict().items()]
            items.sort()
            if len(self):
                # we may need to add an extra element on the end to
                # set the size.  (There doesn't seem to be a better way
                # to do it.)
                if len(items) == 0 or len(self)-1 > items[-1][0]:
                    items.append((len(self)-1, sib.int(0)))
            items_dict = sib.dict([(sib.int(n), e) for n,e in items])

            return sib.name('vector')(self.base_ring(), items_dict)
        else:
            return sib.name('vector')(self.base_ring(),
                                      [sib(e, 2) for e in self])

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        r"""
        Return a numerical approximation of ``self`` with ``prec`` bits
        (or decimal ``digits``) of precision, by approximating all
        entries.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- which algorithm to use to compute the
          approximation of the entries (the accepted algorithms depend
          on the object)

        If neither ``prec`` nor ``digits`` is given, the default
        precision is 53 bits (roughly 16 digits).

        EXAMPLES::

            sage: v = vector(RealField(212), [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: numerical_approx(v)
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision
            sage: numerical_approx(v, digits=3)
            (1.00, 2.00, 3.00)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 14 bits of precision

        Both functional and object-oriented usage is possible.  ::

            sage: u = vector(QQ, [1/2, 1/3, 1/4])
            sage: u.n()
            (0.500000000000000, 0.333333333333333, 0.250000000000000)
            sage: u.numerical_approx()
            (0.500000000000000, 0.333333333333333, 0.250000000000000)
            sage: n(u)
            (0.500000000000000, 0.333333333333333, 0.250000000000000)
            sage: N(u)
            (0.500000000000000, 0.333333333333333, 0.250000000000000)
            sage: numerical_approx(u)
            (0.500000000000000, 0.333333333333333, 0.250000000000000)

        Precision (bits) and digits (decimal) may be specified.
        When both are given, ``prec`` wins.  ::

            sage: u = vector(QQ, [1/2, 1/3, 1/4])
            sage: n(u, prec=15)
            (0.5000, 0.3333, 0.2500)
            sage: n(u, digits=5)
            (0.50000, 0.33333, 0.25000)
            sage: n(u, prec=30, digits=100)
            (0.50000000, 0.33333333, 0.25000000)

        These are some legacy doctests that were part of various specialized
        versions of the numerical approximation routine that were removed as
        part of :issue:`12195`.  ::

            sage: v = vector(ZZ, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision

            sage: v = vector(RDF, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v = vector(CDF, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Complex Field with 53 bits of precision

            sage: v = vector(Integers(8), [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision

            sage: v = vector(QQ, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision

        ::

            sage: v = vector(GF(2), [1,2,3])
            sage: v.n()
            (1.00000000000000, 0.000000000000000, 1.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 0.0000000000000000000000, 1.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision

        TESTS:

        Sparse vectors have a similar method that works efficiently for
        the sparse case.  We test that it is working as it should.  ::

            sage: v = vector(QQ, [1/2, 0, 0, 1/3, 0, 0, 0, 1/4], sparse=True)
            sage: u = v.numerical_approx(digits=4)
            sage: u.is_sparse()
            True
            sage: u
            (0.5000, 0.0000, 0.0000, 0.3333, 0.0000, 0.0000, 0.0000, 0.2500)
        """
        if prec is None:
            prec = digits_to_bits(digits)
        return vector([e.numerical_approx(prec, algorithm=algorithm) for e in self])

    def row(self):
        r"""
        Return a matrix with a single row and the same entries as the vector ``self``.

        OUTPUT:

        A matrix over the same ring as the vector (or free module element), with
        a single row.  The entries of the row are identical to those of the vector,
        and in the same order.

        EXAMPLES::

            sage: v = vector(ZZ, [1,2,3])
            sage: w = v.row(); w
            [1 2 3]
            sage: w.parent()
            Full MatrixSpace of 1 by 3 dense matrices over Integer Ring

            sage: x = vector(FiniteField(13), [2,4,8,16])
            sage: x.row()
            [2 4 8 3]

        There is more than one way to get one-row matrix from a vector,
        but the ``row`` method is more efficient than making a column and
        then taking a transpose.  Notice that supplying a vector to the
        matrix constructor demonstrates Sage's preference for rows. ::

            sage: x = vector(RDF, [sin(i*pi/20) for i in range(10)])                    # needs sage.symbolic
            sage: x.row() == matrix(x)
            True
            sage: x.row() == x.column().transpose()
            True

        Sparse or dense implementations are preserved. ::

            sage: d = vector(RR, [1.0, 2.0, 3.0])
            sage: s = vector(CDF, {2: 5.0+6.0*I})                                       # needs sage.symbolic
            sage: dm = d.row()
            sage: sm = s.row()                                                          # needs sage.symbolic
            sage: all([d.is_dense(), dm.is_dense(), s.is_sparse(), sm.is_sparse()])     # needs sage.symbolic
            True

        TESTS:

        The :meth:`~sage.matrix.matrix1.Matrix.row` method will return
        a specified row of a matrix as a vector.  So here are a couple
        of round-trips. ::

            sage: A = matrix(ZZ, [[1,2,3]])
            sage: A == A.row(0).row()
            True
            sage: v = vector(ZZ, [4,5,6])
            sage: v == v.row().row(0)
            True

        And a very small corner case. ::

            sage: v = vector(ZZ, [])
            sage: w = v.row()
            sage: w.parent()
            Full MatrixSpace of 1 by 0 dense matrices over Integer Ring
        """
        from sage.matrix.args import MatrixArgs
        ma = MatrixArgs(self._parent._base, 1, self.degree(),
                list(self), sparse=self.is_sparse())
        return ma.matrix()

    def column(self):
        r"""
        Return a matrix with a single column and the same entries as the vector ``self``.

        OUTPUT:

        A matrix over the same ring as the vector (or free module element), with
        a single column.  The entries of the column are identical to those of the
        vector, and in the same order.

        EXAMPLES::

            sage: v = vector(ZZ, [1,2,3])
            sage: w = v.column(); w
            [1]
            [2]
            [3]
            sage: w.parent()
            Full MatrixSpace of 3 by 1 dense matrices over Integer Ring

            sage: x = vector(FiniteField(13), [2,4,8,16])
            sage: x.column()
            [2]
            [4]
            [8]
            [3]

        There is more than one way to get one-column matrix from a vector.
        The ``column`` method is about equally efficient to making a row and
        then taking a transpose.  Notice that supplying a vector to the
        matrix constructor demonstrates Sage's preference for rows. ::

            sage: x = vector(RDF, [sin(i*pi/20) for i in range(10)])                    # needs sage.libs.pari sage.symbolic
            sage: x.column() == matrix(x).transpose()
            True
            sage: x.column() == x.row().transpose()
            True

        Sparse or dense implementations are preserved. ::

            sage: d = vector(RR, [1.0, 2.0, 3.0])
            sage: s = vector(CDF, {2: 5.0+6.0*I})                                       # needs sage.symbolic
            sage: dm = d.column()
            sage: sm = s.column()                                                       # needs sage.symbolic
            sage: all([d.is_dense(), dm.is_dense(), s.is_sparse(), sm.is_sparse()])     # needs sage.symbolic
            True

        TESTS:

        The :meth:`~sage.matrix.matrix1.Matrix.column` method will return
        a specified column of a matrix as a vector.  So here are a couple
        of round-trips. ::

            sage: A = matrix(ZZ, [[1],[2],[3]])
            sage: A == A.column(0).column()
            True
            sage: v = vector(ZZ, [4,5,6])
            sage: v == v.column().column(0)
            True

        And a very small corner case. ::

            sage: v = vector(ZZ, [])
            sage: w = v.column()
            sage: w.parent()
            Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        """
        from sage.matrix.args import MatrixArgs
        ma = MatrixArgs(self._parent._base, self.degree(), 1,
                [(x,) for x in self], sparse=self.is_sparse())
        return ma.matrix()

    def __copy__(self):
        """
        Make a copy of this vector.

        EXAMPLES::

            sage: v = vector([1..5]); v
            (1, 2, 3, 4, 5)
            sage: w = copy(v)
            sage: v == w
            True
            sage: v is w
            False

        ::

            sage: v = vector([1..5], sparse=True); v
            (1, 2, 3, 4, 5)
            sage: copy(v)
            (1, 2, 3, 4, 5)
        """
        if self.is_sparse():
            return self.parent()(self.dict())
        else:
            return self.parent()(self.list())

    def subs(self, in_dict=None, **kwds):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('a,b,d,e')
            (a, b, d, e)
            sage: v = vector([a, b, d, e])
            sage: v.substitute(a=1)
            (1, b, d, e)
            sage: v.subs(a=b, b=d)
            (b, d, d, e)
        """
        return self.parent()([ a.subs(in_dict, **kwds) for a in self.list() ])

    def change_ring(self, R):
        """
        Change the base ring of this vector.

        EXAMPLES::

            sage: v = vector(QQ['x,y'], [1..5]); v.change_ring(GF(3))
            (1, 2, 0, 1, 2)

        TESTS:

        Check for :issue:`29630`::

            sage: v = vector(QQ, 4, {0:1}, sparse=True)
            sage: v.change_ring(AA).is_sparse()                                         # needs sage.rings.number_field
            True
        """
        if self.base_ring() is R:
            return self
        M = self._parent.change_ring(R)
        if M.is_sparse():
            return M(self.dict(), coerce=True)
        return M(self.list(), coerce=True)

    def coordinate_ring(self):
        """
        Return the ring from which the coefficients of this vector come.

        This is different from :meth:`base_ring`, which returns the ring
        of scalars.

        EXAMPLES::

            sage: M = (ZZ^2) * (1/2)
            sage: v = M([0,1/2])
            sage: v.base_ring()
            Integer Ring
            sage: v.coordinate_ring()
            Rational Field
        """
        return self._parent.coordinate_ring()

    def additive_order(self):
        """
        Return the additive order of ``self``.

        EXAMPLES::

            sage: v = vector(Integers(4), [1,2])
            sage: v.additive_order()
            4

        ::

            sage: v = vector([1,2,3])
            sage: v.additive_order()
            +Infinity

        ::

            sage: v = vector(Integers(30), [6, 15]); v
            (6, 15)
            sage: v.additive_order()
            10
            sage: 10*v
            (0, 0)
        """
        cdef list v = []
        cdef Py_ssize_t i
        for i in range(self._degree):
            ord = self[i].additive_order()
            if isinstance(ord, AnInfinity):
                return ord
            v.append(ord)
        from sage.arith.functions import lcm
        return lcm(v)

    def items(self):
        """
        Return an iterator over ``self``.

        EXAMPLES::

            sage: v = vector([1,2/3,pi])                                                # needs sage.symbolic
            sage: v.items()                                                             # needs sage.symbolic
            <...generator object at ...>
            sage: list(v.items())                                                       # needs sage.symbolic
            [(0, 1), (1, 2/3), (2, pi)]

        TESTS:

        Using iteritems as an alias::

            sage: list(v.iteritems())                                                   # needs sage.symbolic
            [(0, 1), (1, 2/3), (2, pi)]
        """
        cdef dict d = self.dict(copy=False)
        yield from d.iteritems()

    iteritems = items

    def __abs__(self):
        """
        Return the square root of the sum of the squares of the entries of
        this vector.

        EXAMPLES::

            sage: v = vector([1..5]); abs(v)                                            # needs sage.symbolic
            sqrt(55)
            sage: v = vector(RDF, [1..5]); abs(v)
            7.416198487095663
        """
        return sum([x**2 for x in self.list()]).sqrt()

    def norm(self, p=__two__):
        r"""
        Return the `p`-norm of ``self``.

        INPUT:

        - ``p`` -- (default: 2) ``p`` can be a real number greater than 1,
          infinity (``oo`` or ``Infinity``), or a symbolic expression:

          - `p=1`: the taxicab (Manhattan) norm
          - `p=2`: the usual Euclidean norm (the default)
          - `p=\infty`: the maximum entry (in absolute value)

        .. NOTE::

            See also :func:`sage.misc.functional.norm`

        EXAMPLES::

            sage: v = vector([1,2,-3])
            sage: v.norm(5)                                                             # needs sage.symbolic
            276^(1/5)

        The default is the usual Euclidean norm.  ::

            sage: v.norm()                                                              # needs sage.symbolic
            sqrt(14)
            sage: v.norm(2)                                                             # needs sage.symbolic
            sqrt(14)

        The infinity norm is the maximum size (in absolute value)
        of the entries.  ::

            sage: v.norm(Infinity)
            3
            sage: v.norm(oo)
            3

        Real or symbolic values may be used for ``p``.  ::

            sage: v=vector(RDF,[1,2,3])
            sage: v.norm(5)
            3.077384885394063

            sage: # needs sage.symbolic
            sage: v.norm(pi/2)    # abs tol 1e-15
            4.216595864704748
            sage: _ = var('a b c d p'); v = vector([a, b, c, d])
            sage: v.norm(p)
            (abs(a)^p + abs(b)^p + abs(c)^p + abs(d)^p)^(1/p)

        Notice that the result may be a symbolic expression, owing to
        the necessity of taking a square root (in the default case).
        These results can be converted to numerical values if needed. ::

            sage: v = vector(ZZ, [3,4])
            sage: nrm = v.norm(); nrm
            5
            sage: nrm.parent()
            Rational Field

            sage: # needs sage.symbolic
            sage: v = vector(QQ, [3, 5])
            sage: nrm = v.norm(); nrm
            sqrt(34)
            sage: nrm.parent()
            Symbolic Ring
            sage: numeric = N(nrm); numeric
            5.83095189484...
            sage: numeric.parent()
            Real Field with 53 bits of precision

        TESTS:

        The value of ``p`` must be greater than, or
        equal to, one. ::

            sage: v = vector(QQ, [1,2])
            sage: v.norm(0.99)
            Traceback (most recent call last):
            ...
            ValueError: 0.990000000000000 is not greater than or equal to 1

        Norm works with Python integers (see :issue:`13502`). ::

            sage: v = vector(QQ, [1,2])
            sage: v.norm(int(2))                                                        # needs sage.symbolic
            sqrt(5)
        """
        abs_self = [abs(x) for x in self]
        if p == Infinity:
            return max(abs_self)
        if p < 1:
            raise ValueError("%s is not greater than or equal to 1" % p)

        s = sum(a ** p for a in abs_self)
        return s**(__one__/p)

    cpdef _richcmp_(left, right, int op):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: v = vector(SR, [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(SR, [-1,x,pi,0])
            sage: bool(w < v)
            True
            sage: bool(w > v)
            False

        TESTS::

            sage: F.<y> = PolynomialRing(QQ, 'y')
            sage: type(vector(F, [0]*4, sparse=True))
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
            sage: vector(F, [0,0,0,y]) == vector(F, [0,0,0,y])
            True
            sage: vector(F, [0,0,0,0]) == vector(F, [0,2,0,y])
            False

        Verify that :issue:`33697` is fixed::

            sage: # needs sage.symbolic
            sage: v = vector(SR, [x])
            sage: w = vector(SR, [1])
            sage: v == w
            False
            sage: assume(x > 0)
            sage: v == w
            False
            sage: forget()
        """
        cdef Py_ssize_t i
        for i in range(left._degree):
            lx = left[i]
            rx = right[i]
            if not(lx == rx):
                return richcmp_not_equal(lx, rx, op)
        return rich_to_bool(op, 0)

    def __getitem__(self, i):
        """
        Return `i`-th entry or slice of ``self``.

        EXAMPLES::

            sage: v = sage.modules.free_module_element.FreeModuleElement(QQ^3)
            sage: v.__getitem__(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        cdef Py_ssize_t d = self._degree
        cdef Py_ssize_t start, stop, step, slicelength
        cdef Py_ssize_t n
        cdef list values
        if isinstance(i, slice):
            PySlice_GetIndicesEx(i, d, &start, &stop, &step, &slicelength)
            values = []
            for n in range(slicelength):
                values.append(self.get_unsafe(start + n*step))
            from sage.modules.free_module import FreeModule
            M = FreeModule(self.coordinate_ring(), slicelength, sparse=self.is_sparse())
            return M(values, coerce=False, copy=False)
        else:
            n = i
            if n < 0:
                n += d
            if n < 0 or n >= d:
                raise IndexError("vector index out of range")
            return self.get_unsafe(n)

    cdef get_unsafe(self, Py_ssize_t i):
        """
        Cython function to get the `i`-th entry of this vector.

        Used as building block for a generic ``__getitem__``.
        """
        raise NotImplementedError

    def get(self, i):
        """
        Like ``__getitem__`` but without bounds checking:
        `i` must satisfy ``0 <= i < self.degree``.

        EXAMPLES::

            sage: vector(SR, [1/2,2/5,0]).get(0)                                        # needs sage.symbolic
            1/2
        """
        return self.get_unsafe(i)

    def __setitem__(self, i, value):
        """
        Set the `i`-th entry or slice of ``self`` to ``value``.

        EXAMPLES::

            sage: v = sage.modules.free_module_element.FreeModuleElement(QQ^3)
            sage: v[0] = 5
            Traceback (most recent call last):
            ...
            NotImplementedError

        For derived classes, this works::

            sage: v = vector([1,2/3,8])
            sage: v[0] = 5
            sage: v
            (5, 2/3, 8)
        """
        if self._is_immutable:
            raise ValueError("vector is immutable; please change a copy instead (use copy())")
        cdef Py_ssize_t d = self._degree
        cdef Py_ssize_t start, stop, step, slicelength
        cdef Py_ssize_t n
        cdef list values
        R = self.coordinate_ring()
        if isinstance(i, slice):
            PySlice_GetIndicesEx(i, d, &start, &stop, &step, &slicelength)
            values = [R(x) for x in value]
            if len(values) != slicelength:
                raise IndexError("slice assignment would change dimension")
            for n in range(slicelength):
                self.set_unsafe(start + n*step, values[n])
        else:
            n = i
            if n < 0:
                n += d
            if n < 0 or n >= d:
                raise IndexError("vector index out of range")
            self.set_unsafe(n, R(value))

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        Cython function to set the `i`-th entry of this vector to
        ``value``.

        Used as building block for a generic ``__setitem__``.
        """
        raise NotImplementedError

    def set(self, i, value):
        """
        Like ``__setitem__`` but without type or bounds checking:
        `i` must satisfy ``0 <= i < self.degree`` and ``value`` must be
        an element of the coordinate ring.

        EXAMPLES::

            sage: v = vector(SR, [1/2,2/5,0]); v                                        # needs sage.symbolic
            (1/2, 2/5, 0)
            sage: v.set(2, pi); v                                                       # needs sage.symbolic
            (1/2, 2/5, pi)
        """
        assert value.parent() is self.coordinate_ring()
        self.set_unsafe(i, value)

    def __invert__(self):
        """
        Invert v, which makes no sense, and is hence is not implemented.

        EXAMPLES::

            sage: vector([1,2/3,pi]).__invert__()                                       # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __len__(self):
        """
        EXAMPLES::

            sage: len(sage.modules.free_module_element.FreeModuleElement(QQ^2010))
            2010
        """
        return self._degree

    def __mod__(self, p):
        """
        EXAMPLES::

            sage: V = vector(ZZ, [5, 9, 13, 15])
            sage: V % 7
            (5, 2, 6, 1)
            sage: parent(V % 7)
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        return self.parent()([x % p for x in self.list()], copy=False, coerce=False, check=False)

    def Mod(self, p):
        """
        EXAMPLES::

            sage: V = vector(ZZ, [5, 9, 13, 15])
            sage: V.Mod(7)
            (5, 2, 6, 1)
            sage: parent(V.Mod(7))
            Vector space of dimension 4 over Ring of integers modulo 7
        """
        return self.change_ring(self.base_ring().quotient_ring(p))

    def list(self, copy=True):
        """
        Return list of elements of ``self``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); whether returned list is a
          copy that is safe to change (ignored)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: v = vector([x,y,z], sparse=True)
            sage: type(v)
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
            sage: a = v.list(); a
            [x, y, z]
            sage: a[0] = x*y; v
            (x, y, z)

        The optional argument ``copy`` is ignored::

            sage: a = v.list(copy=False); a
            [x, y, z]
            sage: a[0] = x*y; v
            (x, y, z)
        """
        cdef Py_ssize_t i
        return [self[i] for i in range(self._degree)]

    def list_from_positions(self, positions):
        """
        Return list of elements chosen from this vector using the
        given positions of this vector.

        INPUT:

        - ``positions`` -- iterable of integers

        EXAMPLES::

            sage: v = vector([1, 2/3, pi])                                              # needs sage.symbolic
            sage: v.list_from_positions([0,0,0,2,1])                                    # needs sage.symbolic
            [1, 1, 1, pi, 2/3]
        """
        cdef Py_ssize_t i
        return [self[i] for i in positions]

    def lift(self):
        """
        Lift ``self`` to the cover ring.

        OUTPUT:

        Return a lift of ``self`` to the covering ring of the base ring `R`,
        which is by definition the ring returned by calling
        :meth:`~sage.rings.quotient_ring.QuotientRing_nc.cover_ring`
        on `R`, or just `R` itself if the
        :meth:`~sage.rings.quotient_ring.QuotientRing_nc.cover_ring`
        method is not defined.

        EXAMPLES::

            sage: V = vector(Integers(7), [5, 9, 13, 15]) ; V
            (5, 2, 6, 1)
            sage: V.lift()
            (5, 2, 6, 1)
            sage: parent(V.lift())
            Ambient free module of rank 4 over the principal ideal domain Integer Ring

        If the base ring does not have a cover method, return a copy of the vector::

            sage: W = vector(QQ, [1, 2, 3])
            sage: W1 = W.lift()
            sage: W is W1
            False
            sage: parent(W1)
            Vector space of dimension 3 over Rational Field
        """
        try:
            return self.change_ring(self.base_ring().cover_ring())
        except AttributeError:
            from copy import copy
            return copy(self)

    def lift_centered(self):
        r"""
        Lift to a congruent, centered vector.

        INPUT:

        - ``self`` A vector with coefficients in `Integers(n)`.

        OUTPUT:

        - The unique integer vector `v` such that foreach `i`,
          `Mod(v[i],n) = Mod(self[i],n)` and `-n/2 < v[i] \leq n/2`.

        EXAMPLES::

            sage: V = vector(Integers(7), [5, 9, 13, 15]) ; V
            (5, 2, 6, 1)
            sage: V.lift_centered()
            (-2, 2, -1, 1)
            sage: parent(V.lift_centered())
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        R = self.base_ring().cover_ring()
        l = [foo.lift_centered() for foo in self]
        P = self.parent().change_ring(R)
        return P(l)

    def __pos__(self):
        """
        Always returns ``self``, since ``+self == self``.

        EXAMPLES::

            sage: v = vector([1,2/3,8])
            sage: v.__pos__()
            (1, 2/3, 8)
            sage: +v
            (1, 2/3, 8)
        """
        return self

    def __pow__(self, n, dummy):
        """
        Raises a :exc:`NotImplementedError`, since powering doesn't make
        sense for vectors.

        EXAMPLES::

            sage: v = vector([1,2/3,8])
            sage: v^2
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_(self):
        """
        String representation of a vector.

        EXAMPLES::

            sage: vector(QQ, [])._repr_()
            '()'
            sage: vector(QQ, range(5))._repr_()
            '(0, 1, 2, 3, 4)'

        Symbolic are not displayed using ASCII art.

        ::

            sage: x = var('x')                                                          # needs sage.symbolic
            sage: v = vector([x/(2*x)+sqrt(2)+var('theta')^3,x/(2*x)]); v               # needs sage.symbolic
            (theta^3 + sqrt(2) + 1/2, 1/2)
            sage: v._repr_()                                                            # needs sage.symbolic
            '(theta^3 + sqrt(2) + 1/2, 1/2)'
        """
        cdef Py_ssize_t d = self._degree
        if d == 0:
            return "()"
        # compute column widths
        S = [repr(x) for x in self.list(copy=False)]
        # width = max([len(x) for x in S])
        s = "("
        for i in range(d):
            if i == d-1:
                sep = ""
            else:
                sep=", "
            entry = S[i]
            #if i > 0:
            #    entry = " "*(width-len(entry)) + entry
            s = s + entry + sep
        s = s + ")"
        return s

    def _maple_init_(self):
        """
        EXAMPLES::

            sage: v = vector(ZZ, 4, range(4))
            sage: maple(v)  # optional - maple
            Vector[row](4, [0,1,2,3])

        ::

            sage: v = vector(QQ, 3, [2/3, 0, 5/4])
            sage: maple(v)  # optional - maple
            Vector[row](3, [2/3,0,5/4])

        ::

            sage: P.<x> = ZZ[]
            sage: v = vector(P, 3, [x^2 + 2, 2*x + 1, -2*x^2 + 4*x])
            sage: maple(v)  # optional - maple
            Vector[row](3, [x^2+2,2*x+1,-2*x^2+4*x])
        """
        return "Vector[row](%s)"%(str(self.list()))

    def degree(self):
        """
        Return the degree of this vector, which is simply the number
        of entries.

        EXAMPLES::

            sage: sage.modules.free_module_element.FreeModuleElement(QQ^389).degree()
            389
            sage: vector([1,2/3,8]).degree()
            3
        """
        return self._degree

    def denominator(self):
        """
        Return the least common multiple of the denominators of the
        entries of ``self``.

        EXAMPLES::

            sage: v = vector([1/2,2/5,3/14])
            sage: v.denominator()
            70
            sage: 2*5*7
            70

        ::

            sage: M = (ZZ^2)*(1/2)
            sage: M.basis()[0].denominator()
            2

        TESTS:

        The following was fixed in :issue:`8800`::

            sage: M = GF(5)^3
            sage: v = M((4,0,2))
            sage: v.denominator()
            1
        """
        # It may be that the coordinates do not have a denominator
        # (but if one coordinate has it, they all should have it)
        d = self.coordinate_ring().one()
        try:
            d = d.denominator()
        except AttributeError:
            return d
        for y in self.list():
            d = d.lcm(y.denominator())
        return d

    def dict(self, copy=True):
        """
        Return dictionary of nonzero entries of ``self``.

        More precisely, this returns a dictionary whose keys are indices
        of basis elements in the support of ``self`` and whose values are
        the corresponding coefficients.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``self`` is internally
          represented by a dictionary ``d``, then make a copy of ``d``.
          If ``False``, then this can cause undesired behavior by
          mutating ``d``.

        OUTPUT: Python dictionary

        EXAMPLES::

            sage: v = vector([0,0,0,0,1/2,0,3/14])
            sage: v.dict()
            {4: 1/2, 6: 3/14}
            sage: sorted(v.support())
            [4, 6]

        In some cases, when ``copy=False``, we get back a dangerous
        reference::

            sage: v = vector({0:5, 2:3/7}, sparse=True)
            sage: v.dict(copy=False)
            {0: 5, 2: 3/7}
            sage: v.dict(copy=False)[0] = 18
            sage: v
            (18, 0, 3/7)
        """
        cdef dict e = {}
        cdef Py_ssize_t i
        for i in range(self._degree):
            c = self[i]
            if c:
                e[i] = c
        return e

    #############################
    # Plotting
    #############################
    def plot(self, plot_type=None, start=None, **kwds):
        """
        INPUT:

        - ``plot_type`` -- (default: 'arrow' if v has 3 or fewer components,
            otherwise 'step') type of plot. Options are:

            - 'arrow' to draw an arrow

            - 'point' to draw a point at the coordinates specified by the
              vector

            - 'step' to draw a step function representing the coordinates
              of the vector.

          Both 'arrow' and 'point' raise exceptions if the vector has
          more than 3 dimensions.

        - ``start`` -- (default: origin in correct dimension) may be a tuple,
          list, or vector

        EXAMPLES:

        The following both plot the given vector::

            sage: v = vector(RDF, (1,2))
            sage: A = plot(v)                                                           # needs sage.plot
            sage: B = v.plot()                                                          # needs sage.plot
            sage: A + B  # should just show one vector                                  # needs sage.plot
            Graphics object consisting of 2 graphics primitives

        Examples of the plot types::

            sage: # needs sage.plot
            sage: A = plot(v, plot_type='arrow')
            sage: B = plot(v, plot_type='point', color='green', size=20)
            sage: C = plot(v, plot_type='step') # calls v.plot_step()
            sage: A+B+C
            Graphics object consisting of 3 graphics primitives

        You can use the optional arguments for :meth:`plot_step`::

            sage: eps = 0.1
            sage: plot(v, plot_type='step', eps=eps, xmax=5, hue=0)                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Three-dimensional examples::

            sage: v = vector(RDF, (1,2,1))
            sage: plot(v)  # defaults to an arrow plot                                  # needs sage.plot
            Graphics3d Object

        ::

            sage: plot(v, plot_type='arrow')                                            # needs sage.plot
            Graphics3d Object

        ::

            sage: from sage.plot.plot3d.shapes2 import frame3d                          # needs sage.plot
            sage: plot(v, plot_type='point')+frame3d((0,0,0), v.list())                 # needs sage.plot
            Graphics3d Object

        ::

            sage: plot(v, plot_type='step')  # calls v.plot_step()                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: plot(v, plot_type='step', eps=eps, xmax=5, hue=0)                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        With greater than three coordinates, it defaults to a step plot::

            sage: v = vector(RDF, (1,2,3,4))
            sage: plot(v)                                                               # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        One dimensional vectors are plotted along the horizontal axis of
        the coordinate plane::

            sage: plot(vector([1]))                                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        An optional start argument may also be specified by a tuple, list, or vector::

            sage: u = vector([1,2]); v = vector([2,5])
            sage: plot(u, start=v)                                                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        TESTS::

            sage: u = vector([1,1]); v = vector([2,2,2]); z=(3,3,3)
            sage: plot(u)  #test when start=None                                        # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: # needs sage.plot
            sage: plot(u, start=v) #test when coordinate dimension mismatch exists
            Traceback (most recent call last):
            ...
            ValueError: vector coordinates are not of the same dimension
            sage: P = plot(v, start=z)        # test when start coordinates are passed as a tuple
            sage: P = plot(v, start=list(z))  # test when start coordinates are passed as a list
        """
        # Give sensible defaults based on the vector length
        if plot_type is None:
            if len(self)<=3:
                plot_type='arrow'
            else:
                plot_type='step'

        coords = self.list()

        if start is None:
            start = [0]*len(coords)
        elif len(start)!=len(coords):
            raise ValueError("vector coordinates are not of the same dimension")
        else:
            start = list(start)

        if plot_type == 'arrow' or plot_type == 'point':
            dimension = len(coords)
            if dimension == 3:
                from sage.plot.plot3d.shapes2 import line3d, point3d

                if plot_type == 'arrow':
                    return line3d([start, [(u+v) for u,v in zip(coords, start)]], arrow_head=True, **kwds)
                else:
                    return point3d(coords, **kwds)
            elif dimension < 3:
                if dimension < 2:
                    # pad to make 2-dimensional
                    coords.extend([0]*(2-dimension))
                    start.extend([0]*(2-dimension))

                from sage.plot.all import arrow, point
                if plot_type == 'arrow':
                    return arrow(start, [(u+v) for u,v in zip(coords, start)], **kwds)
                else:
                    return point(coords, **kwds)
            else:
                raise ValueError("arrow and point plots require vectors with 3 or fewer components")

        elif plot_type == 'step':
            return self.plot_step(**kwds)
        else:
            raise NotImplementedError("plot_type was unrecognized")

    def plot_step(self, xmin=0, xmax=1, eps=None, res=None,
             connect=True, **kwds):
        r"""
        INPUT:

        - ``xmin`` -- (default: 0) start x position to start
          plotting

        - ``xmax`` -- (default: 1) stop x position to stop
          plotting

        - ``eps`` -- (default: determined by xmax) we view this
          vector as defining a function at the points xmin, xmin + eps, xmin
          + 2\*eps, ...,

        - ``res`` -- (default: all points) total number of
          points to include in the graph

        - ``connect`` -- boolean (default: ``True``); if ``True`` draws a line,
          otherwise draw a list of points

        EXAMPLES::

            sage: eps = 0.1
            sage: v = vector(RDF, [sin(n*eps) for n in range(100)])
            sage: v.plot_step(eps=eps, xmax=5, hue=0)                                   # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        import math
        if res is None:
            res = self.degree()
        if eps is None:
            eps = float(xmax - xmin)/res
        v = []
        x = xmin
        for i in range(0, self.degree(), int(math.ceil(self.degree()/res))):
            y = float(self[i])
            if x > xmax:
                break
            v.append((x,y))
            x += eps
            v.append((x,y))
        from sage.plot.all import line, points
        if connect:
            return line(v, **kwds)
        else:
            return points(v, **kwds)

    cpdef _dot_product_coerce_(left, Vector right):
        """
        Return the dot product of left and right.

        This function works even if the parents are different, the
        degrees have to match however.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: v = vector(RDF, [0,1,2])
            sage: w = vector(R, [x,0,0])
            sage: p = v._dot_product_coerce_(w)
            sage: p
            0
            sage: parent(p)
            Univariate Polynomial Ring in x over Real Double Field

        Zero-dimensional vectors also work correctly::

            sage: v = vector(RDF, [])
            sage: w = vector(R, [])
            sage: parent(v._dot_product_coerce_(w))
            Univariate Polynomial Ring in x over Real Double Field
        """
        if left._degree == 0:
            return (left.coordinate_ring().zero()
                    * right.coordinate_ring().zero())
        cdef list a = left.list(copy=False)
        cdef list b = right.list(copy=False)
        cdef Py_ssize_t i
        z = a[0] * b[0]
        for i in range(1, left._degree):
            z += a[i] * b[i]
        return z

    def dot_product(self, right):
        r"""
        Return the dot product of ``self`` and ``right``, which is the
        sum of the product of the corresponding entries.

        INPUT:

        - ``right`` -- a vector of the same degree as ``self``.
          It does not need to belong to the same parent as ``self``,
          so long as the necessary products and sums are defined.

        OUTPUT:

        If ``self`` and ``right`` are the vectors `\vec{x}` and `\vec{y}`,
        of degree `n`, then this method returns

        .. MATH::

            \sum_{i=1}^{n}x_iy_i

        .. NOTE::

            The :meth:`inner_product` is a more general version of
            this method, and the :meth:`hermitian_inner_product`
            method may be more appropriate if your vectors
            have complex entries.

        EXAMPLES::

            sage: V = FreeModule(ZZ, 3)
            sage: v = V([1,2,3])
            sage: w = V([4,5,6])
            sage: v.dot_product(w)
            32

        ::

            sage: R.<x> = QQ[]
            sage: v = vector([x,x^2,3*x]); w = vector([2*x,x,3+x])
            sage: v*w
            x^3 + 5*x^2 + 9*x
            sage: (x*2*x) + (x^2*x) + (3*x*(3+x))
            x^3 + 5*x^2 + 9*x
            sage: w*v
            x^3 + 5*x^2 + 9*x

        The vectors may be from different vector spaces,
        provided the necessary operations make sense.
        Notice that coercion will generate a result of
        the same type, even if the order of the
        arguments is reversed.::

            sage: v = vector(ZZ, [1,2,3])
            sage: w = vector(FiniteField(3), [0,1,2])
            sage: ip = w.dot_product(v); ip
            2
            sage: ip.parent()
            Finite Field of size 3

            sage: ip = v.dot_product(w); ip
            2
            sage: ip.parent()
            Finite Field of size 3

        The dot product of a vector with itself is the 2-norm, squared. ::

            sage: v = vector(QQ, [3, 4, 7])
            sage: v.dot_product(v) - v.norm()^2                                         # needs sage.symbolic
            0

        TESTS:

        The second argument must be a free module element. ::

            sage: v = vector(QQ, [1,2])
            sage: v.dot_product('junk')
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert str to sage.modules.free_module_element.FreeModuleElement

        The degrees of the arguments must match. ::

            sage: v = vector(QQ, [1,2])
            sage: w = vector(QQ, [1,2,3])
            sage: v.dot_product(w)
            Traceback (most recent call last):
            ...
            ArithmeticError: degrees (2 and 3) must be the same

        Check that vectors with different base rings play out nicely (:issue:`3103`)::

            sage: vector(CDF, [2, 2]) * vector(ZZ, [1, 3])
            8.0

        Zero-dimensional vectors work::

            sage: v = vector(ZZ, [])
            sage: v.dot_product(v)
            0

        TESTS:

        Check for :issue:`33814`::

            sage: rings = [ZZ, QQ, RDF, ZZ['x']]
            sage: rings += [RR]                                                         # needs sage.rings.real_mpfr
            sage: rings += [GF(2), GF(3)]
            sage: rings += [GF(4)]                                                      # needs sage.rings.finite_rings
            sage: for R in rings:
            ....:     _ = (R**0)().dot_product((R**0)())
        """
        cdef FreeModuleElement r = <FreeModuleElement?>right
        if self._parent is r._parent:
            # If the parents are equal, the degree is also equal
            if self._degree == 0:
                return self._parent.coordinate_ring().zero()
            return self._dot_product_(r)
        if self._degree != r._degree:
            raise ArithmeticError("degrees (%s and %s) must be the same"%(self.degree(), right.degree()))
        # Base rings are not equal => use dot product with coercion
        return self._dot_product_coerce_(r)

    def cross_product(self, right):
        """
        Return the cross product of ``self`` and ``right``, which is only defined
        for vectors of length 3 or 7.

        INPUT:

        - ``right`` -- a vector of the same size as ``self``, either
          degree three or degree seven

        OUTPUT:

        The cross product (vector product) of ``self`` and ``right``,
        a vector of the same size of ``self`` and ``right``.

        This product is performed under the assumption that the basis
        vectors are orthonormal. See the method
        :meth:`~sage.manifolds.differentiable.vectorfield.VectorField.cross_product`
        of vector fields for more general cases.

        EXAMPLES::

            sage: v = vector([1,2,3]); w = vector([0,5,-9])
            sage: v.cross_product(v)
            (0, 0, 0)
            sage: u = v.cross_product(w); u
            (-33, 9, 5)
            sage: u.dot_product(v)
            0
            sage: u.dot_product(w)
            0

        The cross product is defined for degree seven vectors as well:
        see :wikipedia:`Cross_product`.
        The 3-D cross product is achieved using the quaternions,
        whereas the 7-D cross product is achieved using the octonions. ::

            sage: u = vector(QQ, [1, -1/3, 57, -9, 56/4, -4,1])
            sage: v = vector(QQ, [37, 55, -99/57, 9, -12, 11/3, 4/98])
            sage: u.cross_product(v)
            (1394815/2793, -2808401/2793, 39492/49, -48737/399, -9151880/2793, 62513/2793, -326603/171)

        The degree seven cross product is anticommutative. ::

            sage: u.cross_product(v) + v.cross_product(u)
            (0, 0, 0, 0, 0, 0, 0)

        The degree seven cross product is distributive across addition. ::

            sage: v = vector([-12, -8/9, 42, 89, -37, 60/99, 73])
            sage: u = vector([31, -42/7, 97, 80, 30/55, -32, 64])
            sage: w = vector([-25/4, 40, -89, -91, -72/7, 79, 58])
            sage: v.cross_product(u + w) - (v.cross_product(u) + v.cross_product(w))
            (0, 0, 0, 0, 0, 0, 0)

        The degree seven cross product respects scalar multiplication. ::

            sage: v = vector([2, 17, -11/5, 21, -6, 2/17, 16])
            sage: u = vector([-8, 9, -21, -6, -5/3, 12, 99])
            sage: (5*v).cross_product(u) - 5*(v.cross_product(u))
            (0, 0, 0, 0, 0, 0, 0)
            sage: v.cross_product(5*u) - 5*(v.cross_product(u))
            (0, 0, 0, 0, 0, 0, 0)
            sage: (5*v).cross_product(u) - (v.cross_product(5*u))
            (0, 0, 0, 0, 0, 0, 0)

        The degree seven cross product respects the scalar triple product. ::

            sage: v = vector([2,6,-7/4,-9/12,-7,12,9])
            sage: u = vector([22,-7,-9/11,12,15,15/7,11])
            sage: w = vector([-11,17,19,-12/5,44,21/56,-8])
            sage: v.dot_product(u.cross_product(w)) - w.dot_product(v.cross_product(u))
            0

        TESTS:

        Both vectors need to be of length three or both vectors need to be of length seven. ::

            sage: u = vector(range(7))
            sage: v = vector(range(3))
            sage: u.cross_product(v)
            Traceback (most recent call last):
            ...
            TypeError: Cross product only defined for vectors of length three or seven, not (7 and 3)

        AUTHOR:

        Billy Wonderly (2010-05-11), Added 7-D Cross Product
        """
        if not isinstance(right, FreeModuleElement):
            raise TypeError("right must be a free module element")
        r = right.list(copy=False)
        l = self.list(copy=False)
        if len(r) == 3 and len(l) == 3:
            return vector([l[1]*r[2] - l[2]*r[1],
                           l[2]*r[0] - l[0]*r[2],
                           l[0]*r[1] - l[1]*r[0]])

        elif len(r) == 7 and len(l) == 7:
            return vector([l[1]*r[3] - l[3]*r[1] + l[2]*r[6] - l[6]*r[2] + l[4]*r[5] - l[5]*r[4],
                           l[2]*r[4] - l[4]*r[2] + l[3]*r[0] - l[0]*r[3] + l[5]*r[6] - l[6]*r[5],
                           l[3]*r[5] - l[5]*r[3] + l[4]*r[1] - l[1]*r[4] + l[6]*r[0] - l[0]*r[6],
                           l[4]*r[6] - l[6]*r[4] + l[5]*r[2] - l[2]*r[5] + l[0]*r[1] - l[1]*r[0],
                           l[5]*r[0] - l[0]*r[5] + l[6]*r[3] - l[3]*r[6] + l[1]*r[2] - l[2]*r[1],
                           l[6]*r[1] - l[1]*r[6] + l[0]*r[4] - l[4]*r[0] + l[2]*r[3] - l[3]*r[2],
                           l[0]*r[2] - l[2]*r[0] + l[1]*r[5] - l[5]*r[1] + l[3]*r[4] - l[4]*r[3]])

        else:
            raise TypeError("Cross product only defined for vectors of length three or seven, not (%s and %s)" % (len(l), len(r)))

    def cross_product_matrix(self):
        r"""
        Return the matrix which describes a cross product
        between ``self`` and some other vector.

        This operation is sometimes written using the hat operator:
        see :wikipedia:`Hat_operator#Cross_product`.
        It is only defined for vectors of length 3 or 7.
        For a vector `v` the cross product matrix `\hat{v}`
        is a matrix which satisfies `\hat{v} \cdot w = v \times w`
        and also `w \cdot \hat{v} = w \times v` for all vectors `w`.
        The basis vectors are assumed to be orthonormal.

        OUTPUT: the cross product matrix of this vector

        EXAMPLES::

            sage: v = vector([1, 2, 3])
            sage: vh = v.cross_product_matrix()
            sage: vh
            [ 0 -3  2]
            [ 3  0 -1]
            [-2  1  0]
            sage: w = random_vector(3, x=1, y=100)
            sage: vh*w == v.cross_product(w)
            True
            sage: w*vh == w.cross_product(v)
            True
            sage: vh.is_alternating()
            True

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(previous_prime(2^32))
            sage: v = random_vector(F, 3)
            sage: w = random_vector(F, 3)
            sage: vh = v.cross_product_matrix()
            sage: vh*w == v.cross_product(w)
            True
            sage: w*vh == w.cross_product(v)
            True
            sage: vh.is_alternating()
            True
            sage: v = random_vector(F, 7)
            sage: w = random_vector(F, 7)
            sage: vh = v.cross_product_matrix()
            sage: vh*w == v.cross_product(w)
            True
            sage: w*vh == w.cross_product(v)
            True
            sage: vh.is_alternating()
            True
            sage: random_vector(F, 5).cross_product_matrix()
            Traceback (most recent call last):
            ...
            TypeError: Cross product only defined for vectors of length three or seven, not 5
        """
        from sage.matrix.matrix_space import MatrixSpace
        rank = self.parent().rank()
        R = self.base_ring()
        zero = R.zero()
        if rank == 3:
            MS = MatrixSpace(R, rank, rank, sparse=self.is_sparse())
            s = self.list(copy=False)
            return MS([
                [ zero, -s[2],  s[1]],
                [ s[2],  zero, -s[0]],
                [-s[1],  s[0],  zero]])
        elif rank == 7:
            MS = MatrixSpace(R, rank, rank, sparse=self.is_sparse())
            s = self.list(copy=False)
            return MS([
                [ zero, -s[3], -s[6],  s[1], -s[5],  s[4],  s[2]],
                [ s[3],  zero, -s[4], -s[0],  s[2], -s[6],  s[5]],
                [ s[6],  s[4],  zero, -s[5], -s[1],  s[3], -s[0]],
                [-s[1],  s[0],  s[5],  zero, -s[6], -s[2],  s[4]],
                [ s[5], -s[2],  s[1],  s[6],  zero, -s[0], -s[3]],
                [-s[4],  s[6], -s[3],  s[2],  s[0],  zero, -s[1]],
                [-s[2], -s[5],  s[0], -s[4],  s[3],  s[1],  zero]])
        else:
            raise TypeError("Cross product only defined for vectors of length three or seven, not {}".format(rank))

    def pairwise_product(self, right):
        """
        Return the pairwise product of ``self`` and ``right``, which is a vector of
        the products of the corresponding entries.

        INPUT:

        - ``right`` -- vector of the same degree as ``self``. It
          need not be in the same vector space as ``self``, as long as the
          coefficients can be multiplied.

        EXAMPLES::

            sage: V = FreeModule(ZZ, 3)
            sage: v = V([1,2,3])
            sage: w = V([4,5,6])
            sage: v.pairwise_product(w)
            (4, 10, 18)
            sage: sum(v.pairwise_product(w)) == v.dot_product(w)
            True

        ::

            sage: W = VectorSpace(GF(3), 3)
            sage: w = W([0,1,2])
            sage: w.pairwise_product(v)
            (0, 2, 0)
            sage: w.pairwise_product(v).parent()
            Vector space of dimension 3 over Finite Field of size 3

        Implicit coercion is well defined (regardless of order), so we
        get 2 even if we do the dot product in the other order.

        ::

            sage: v.pairwise_product(w).parent()
            Vector space of dimension 3 over Finite Field of size 3

        TESTS::

        ::

            sage: parent(vector(ZZ,[1,2]).pairwise_product(vector(ZZ,[1,2])))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(vector(ZZ,[1,2]).pairwise_product(vector(QQ,[1,2])))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2]).pairwise_product(vector(ZZ,[1,2])))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2]).pairwise_product(vector(QQ,[1,2])))
            Vector space of dimension 2 over Rational Field

        ::

            sage: parent(vector(QQ,[1,2,3,4]).pairwise_product(vector(ZZ['x'],[1,2,3,4])))
            Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x'],[1,2,3,4]).pairwise_product(vector(QQ,[1,2,3,4])))
            Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

        ::

            sage: parent(vector(QQ,[1,2,3,4]).pairwise_product(vector(ZZ['x']['y'],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain
             Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2,3,4]).pairwise_product(vector(QQ,[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain
             Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: parent(vector(QQ['x'],[1,2,3,4]).pairwise_product(vector(ZZ['x']['y'],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain
             Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2,3,4]).pairwise_product(vector(QQ['x'],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain
             Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: parent(vector(QQ['y'],[1,2,3,4]).pairwise_product(vector(ZZ['x']['y'],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain
             Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2,3,4]).pairwise_product(vector(QQ['y'],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain
             Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: parent(vector(ZZ['x'],[1,2,3,4]).pairwise_product(vector(ZZ['y'],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and
            'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ['x'],[1,2,3,4]).pairwise_product(vector(QQ['y'],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and
            'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ['x'],[1,2,3,4]).pairwise_product(vector(ZZ['y'],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and
            'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ['x'],[1,2,3,4]).pairwise_product(vector(QQ['y'],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and
            'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: v = vector({1: 1, 3: 2})  # test sparse vectors
            sage: w = vector({0: 6, 3: -4})
            sage: v.pairwise_product(w)
            (0, 0, 0, -8)
            sage: w.pairwise_product(v) == v.pairwise_product(w)
            True
        """
        if not isinstance(right, FreeModuleElement):
            raise TypeError("right must be a free module element")
        if self._parent is not (<FreeModuleElement>right)._parent:
            self, right = canonical_coercion(self, right)
        return self._pairwise_product_(right)

    def _variables(self):
        """
        Return the ordered variable of self, as defined by the basering.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: vector([x, y, 3])._variables()
            [x, y, z]
            sage: vector(SR, [x, y, 3])._variables()                                    # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: Unable to determine ordered variable names for Symbolic Ring
            sage: v(x, y, z) = (-y, x, 0)                                               # needs sage.symbolic
            sage: v._variables()                                                        # needs sage.symbolic
            [(x, y, z) |--> x, (x, y, z) |--> y, (x, y, z) |--> z]
        """
        R = self._parent.base_ring()
        try:
            var_names = R.variable_names()
        except ValueError:
            if hasattr(R, 'arguments'):
                var_names = R.arguments()
            else:
                raise ValueError("Unable to determine ordered variable names for %s" % R)
        return [R(x) for x in var_names]

    def div(self, variables=None):
        """
        Return the divergence of this vector function.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: vector([x, y, z]).div()
            3
            sage: vector([x*y, y*z, z*x]).div()
            x + y + z

            sage: R.<x,y,z,w> = QQ[]
            sage: vector([x*y, y*z, z*x]).div([x, y, z])
            x + y + z
            sage: vector([x*y, y*z, z*x]).div([z, x, y])
            0
            sage: vector([x*y, y*z, z*x]).div([x, y, w])
            y + z

            sage: vector(SR, [x*y, y*z, z*x]).div()                                     # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: Unable to determine ordered variable names for Symbolic Ring
            sage: vector(SR, [x*y, y*z, z*x]).div([x, y, z])                            # needs sage.symbolic
            x + y + z

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.divergence`
            of vector fields on Euclidean spaces (and more generally
            pseudo-Riemannian manifolds), in particular for computing the
            divergence in curvilinear coordinates.
        """
        if variables is None:
            variables = self._variables()
        if len(variables) != len(self):
            raise ValueError("number of variables must equal dimension of self")
        return sum(c.derivative(x) for (c, x) in zip(self, variables))

    def curl(self, variables=None):
        """
        Return the curl of this two-dimensional or three-dimensional
        vector function.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: vector([-y, x, 0]).curl()
            (0, 0, 2)
            sage: vector([y, -x, x*y*z]).curl()
            (x*z, -y*z, -2)
            sage: vector([y^2, 0, 0]).curl()
            (0, 0, -2*y)
            sage: (R^3).random_element().curl().div()
            0

        For rings where the variable order is not well defined, it must be
        defined explicitly::

            sage: v = vector(SR, [-y, x, 0])                                            # needs sage.symbolic
            sage: v.curl()                                                              # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: Unable to determine ordered variable names for Symbolic Ring
            sage: v.curl([x, y, z])                                                     # needs sage.symbolic
            (0, 0, 2)

        Note that callable vectors have well defined variable orderings::

            sage: v(x, y, z) = (-y, x, 0)                                               # needs sage.symbolic
            sage: v.curl()                                                              # needs sage.symbolic
            (x, y, z) |--> (0, 0, 2)

        In two dimensions, this returns a scalar value::

            sage: R.<x,y> = QQ[]
            sage: vector([-y, x]).curl()
            2

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.vectorfield.VectorField.curl`
            of vector fields on Euclidean spaces (and more generally
            pseudo-Riemannian manifolds), in particular for computing the curl
            in curvilinear coordinates.
        """
        if len(self) == 3:
            if variables is None:
                variables = self._variables()
            if len(variables) != 3:
                raise ValueError("exactly 3 variables must be provided")
            x, y, z = variables
            Fx, Fy, Fz = self
            return self.parent([Fz.derivative(y) - Fy.derivative(z),
                                Fx.derivative(z) - Fz.derivative(x),
                                Fy.derivative(x) - Fx.derivative(y)])

        if len(self) == 2:
            if variables is None:
                variables = self._variables()
            if len(variables) != 2:
                raise ValueError("exactly 2 variables must be provided")
            x, y = variables
            Fx, Fy = self
            return Fy.derivative(x) - Fx.derivative(y)

        raise TypeError("curl only defined for 2 or 3 dimensions")

    def element(self):
        """
        Simply return ``self``.  This is useful, since for many objects,
        ``self.element()`` returns a vector corresponding to ``self``.

        EXAMPLES::

            sage: v = vector([1/2,2/5,0]); v
            (1/2, 2/5, 0)
            sage: v.element()
            (1/2, 2/5, 0)
        """
        return self

    def monic(self):
        """
        Return this vector divided through by the first nonzero entry of
        this vector.

        EXAMPLES::

            sage: v = vector(QQ, [0, 4/3, 5, 1, 2])
            sage: v.monic()
            (0, 1, 15/4, 3/4, 3/2)
            sage: v = vector(QQ, [])
            sage: v.monic()
            ()
        """
        cdef Py_ssize_t i
        for i in range(self._degree):
            if self[i]:
                return (~self[i]) * self
        return self

    def normalized(self, p=__two__):
        """
        Return the input vector divided by the p-norm.

        INPUT:

        - ``p`` -- (default: 2) p value for the norm

        EXAMPLES::

            sage: v = vector(QQ, [4, 1, 3, 2])
            sage: v.normalized()                                                        # needs sage.symbolic
            (2/15*sqrt(30), 1/30*sqrt(30), 1/10*sqrt(30), 1/15*sqrt(30))
            sage: sum(v.normalized(1))
            1

        Note that normalizing the vector may change the base ring::

            sage: v.base_ring() == v.normalized().base_ring()                           # needs sage.symbolic
            False
            sage: u = vector(RDF, [-3, 4, 6, 9])
            sage: u.base_ring() == u.normalized().base_ring()
            True
        """
        return self / self.norm(p)

    def conjugate(self):
        r"""
        Return a vector where every entry has been replaced by its complex conjugate.

        OUTPUT:

        A vector of the same length, over the same ring,
        but with each entry replaced by the complex conjugate, as
        implemented by the ``conjugate()`` method for elements of
        the base ring, which is presently always complex conjugation.

        EXAMPLES::

            sage: v = vector(CDF, [2.3 - 5.4*I, -1.7 + 3.6*I])                          # needs sage.symbolic
            sage: w = v.conjugate(); w                                                  # needs sage.symbolic
            (2.3 + 5.4*I, -1.7 - 3.6*I)
            sage: w.parent()                                                            # needs sage.symbolic
            Vector space of dimension 2 over Complex Double Field

        Even if conjugation seems nonsensical over a certain ring, this
        method for vectors cooperates silently. ::

            sage: u = vector(ZZ, range(6))
            sage: u.conjugate()
            (0, 1, 2, 3, 4, 5)

        Sage implements a few specialized subfields of the complex numbers,
        such as the cyclotomic fields.  This example uses such a field
        containing a primitive 7-th root of unity named ``a``. ::

            sage: # needs sage.rings.number_field
            sage: F.<a> = CyclotomicField(7)
            sage: v = vector(F, [a^i for i in range(7)])
            sage: v
            (1, a, a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1)
            sage: v.conjugate()
            (1, -a^5 - a^4 - a^3 - a^2 - a - 1, a^5, a^4, a^3, a^2, a)

        Sparse vectors are returned as such. ::

            sage: # needs sage.symbolic
            sage: v = vector(CC, {1: 5 - 6*I, 3: -7*I}); v
            (0.000000000000000, 5.00000000000000 - 6.00000000000000*I, 0.000000000000000, -7.00000000000000*I)
            sage: v.is_sparse()
            True
            sage: vc = v.conjugate(); vc
            (0.000000000000000, 5.00000000000000 + 6.00000000000000*I, 0.000000000000000, 7.00000000000000*I)
            sage: vc.conjugate()
            (0.000000000000000, 5.00000000000000 - 6.00000000000000*I, 0.000000000000000, -7.00000000000000*I)

        TESTS::

            sage: n = 15
            sage: x = vector(CDF, [sin(i*pi/n)+cos(i*pi/n)*I for i in range(n)])        # needs sage.symbolic
            sage: x + x.conjugate() in RDF^n                                            # needs sage.symbolic
            True
            sage: I*(x - x.conjugate()) in RDF^n                                        # needs sage.symbolic
            True

        The parent of the conjugate is the same as that of the original vector.
        We test this by building a specialized vector space with a non-standard
        inner product, and constructing a test vector in this space. ::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: V = VectorSpace(CDF, 2, inner_product_matrix=[[2,1],[1,5]])
            sage: v = vector(CDF, [2-3*I, 4+5*I])
            sage: w = V(v)
            sage: w.parent()
            Ambient quadratic space of dimension 2 over Complex Double Field
            Inner product matrix:
            [2.0 1.0]
            [1.0 5.0]
            sage: w.conjugate().parent()
            Ambient quadratic space of dimension 2 over Complex Double Field
            Inner product matrix:
            [2.0 1.0]
            [1.0 5.0]
        """
        V = self.parent()
        R = self.base_ring()
        if self.is_sparse():
            # this could be a dictionary comprehension in Python 3
            entries = {}
            for index, entry in self.iteritems():
                entries[index] = entry.conjugate()
        else:
            entries = [entry.conjugate() for entry in self]
        return V(vector(R, self._degree, entries))

    def inner_product(self, right):
        r"""
        Return the inner product of ``self`` and ``right``,
        possibly using an inner product matrix from the parent of ``self``.

        INPUT:

        - ``right`` -- a vector of the same degree as ``self``

        OUTPUT:

        If the parent vector space does not have an inner product
        matrix defined, then this is the usual dot product
        (:meth:`dot_product`).  If ``self`` and ``right`` are
        considered as single column matrices, `\vec{x}` and `\vec{y}`,
        and `A` is the inner product matrix, then this method computes

        .. MATH::

            \left(\vec{x}\right)^tA\vec{y}

        where `t` indicates the transpose.

        .. NOTE::

            If your vectors have complex entries, the
            :meth:`hermitian_inner_product` may be more
            appropriate for your purposes.

        EXAMPLES::

            sage: v = vector(QQ, [1,2,3])
            sage: w = vector(QQ, [-1,2,-3])
            sage: v.inner_product(w)
            -6
            sage: v.inner_product(w) == v.dot_product(w)
            True

        The vector space or free module that is the parent to
        ``self`` can have an inner product matrix defined, which
        will be used by this method.  This matrix will be passed
        through to subspaces. ::

            sage: ipm = matrix(ZZ,[[2,0,-1], [0,2,0], [-1,0,6]])
            sage: M = FreeModule(ZZ, 3, inner_product_matrix=ipm)
            sage: v = M([1,0,0])
            sage: v.inner_product(v)
            2
            sage: K = M.span_of_basis([[0/2,-1/2,-1/2], [0,1/2,-1/2], [2,0,0]])
            sage: (K.0).inner_product(K.0)
            2
            sage: w = M([1,3,-1])
            sage: v = M([2,-4,5])
            sage: w.row()*ipm*v.column() == w.inner_product(v)
            True

        Note that the inner product matrix comes from the parent of ``self``.
        So if a vector is not an element of the correct parent, the result
        could be a source of confusion.  ::

            sage: V = VectorSpace(QQ, 2, inner_product_matrix=[[1,2],[2,1]])
            sage: v = V([12, -10])
            sage: w = vector(QQ, [10,12])
            sage: v.inner_product(w)
            88
            sage: w.inner_product(v)
            0
            sage: w = V(w)
            sage: w.inner_product(v)
            88

        .. NOTE::

            The use of an inner product matrix makes no restrictions on
            the nature of the matrix.  In particular, in this context it
            need not be Hermitian and positive-definite (as it is in the
            example above).

        TESTS:

        Most error handling occurs in the :meth:`dot_product` method.
        But with an inner product defined, this method will check
        that the input is a vector or free module element. ::

            sage: W = VectorSpace(RDF, 2, inner_product_matrix=matrix(RDF, 2, [1.0,2.0,3.0,4.0]))
            sage: v = W([2.0, 4.0])
            sage: v.inner_product(5)
            Traceback (most recent call last):
            ...
            TypeError: right must be a free module element
        """
        if self.parent().is_ambient() and self.parent()._inner_product_is_dot_product():
            return self.dot_product(right)
        if not isinstance(right, FreeModuleElement):
            raise TypeError("right must be a free module element")
        M = self.parent()
        if M.is_ambient() or M.uses_ambient_inner_product():
            A = M.ambient_module().inner_product_matrix()
            return A.linear_combination_of_rows(self).dot_product(right)
        else:
            A = M.inner_product_matrix()
            v = M.coordinate_vector(self)
            w = M.coordinate_vector(right)
            return A.linear_combination_of_rows(v).dot_product(w)

    def outer_product(self, right):
        r"""
        Return a matrix, the outer product of two vectors ``self`` and ``right``.

        INPUT:

        - ``right`` -- a vector (or free module element) of any size, whose
          elements are compatible (with regard to multiplication) with the
          elements of ``self``.

        OUTPUT:

        The outer product of two vectors `x` and `y` (respectively
        ``self`` and ``right``) can be described several ways.  If we
        interpret `x` as a `m\times 1` matrix and interpret `y` as a
        `1\times n` matrix, then the outer product is the `m\times n`
        matrix from the usual matrix product `xy`.  Notice how this
        is the "opposite" in some ways from an inner product (which
        would require `m=n`).

        If we just consider vectors, use each entry of `x` to create
        a scalar multiples of the vector `y` and use these vectors as
        the rows of a matrix.  Or use each entry of `y` to create a
        scalar multiples of `x` and use these vectors as the columns
        of a matrix.

        EXAMPLES::

            sage: u = vector(QQ, [1/2, 1/3, 1/4, 1/5])
            sage: v = vector(ZZ, [60, 180, 600])
            sage: u.outer_product(v)
            [ 30  90 300]
            [ 20  60 200]
            [ 15  45 150]
            [ 12  36 120]
            sage: M = v.outer_product(u); M
            [ 30  20  15  12]
            [ 90  60  45  36]
            [300 200 150 120]
            sage: M.parent()
            Full MatrixSpace of 3 by 4 dense matrices over Rational Field

        The more general :meth:`sage.matrix.matrix2.tensor_product` is an
        operation on a pair of matrices.  If we construct a pair of vectors
        as a column vector and a row vector, then an outer product and a
        tensor product are identical.  Thus `tensor_product` is a synonym
        for this method.  ::

            sage: u = vector(QQ, [1/2, 1/3, 1/4, 1/5])
            sage: v = vector(ZZ, [60, 180, 600])
            sage: u.tensor_product(v) == (u.column()).tensor_product(v.row())
            True

        The result is always a dense matrix, no matter if the two
        vectors are, or are not, dense.  ::

            sage: d = vector(ZZ,[4,5], sparse=False)
            sage: s = vector(ZZ, [1,2,3], sparse=True)
            sage: dd = d.outer_product(d)
            sage: ds = d.outer_product(s)
            sage: sd = s.outer_product(d)
            sage: ss = s.outer_product(s)
            sage: all([dd.is_dense(), ds.is_dense(), sd.is_dense(), dd.is_dense()])
            True

        Vectors with no entries do the right thing.  ::

            sage: v = vector(ZZ, [])
            sage: z = v.outer_product(v)
            sage: z.parent()
            Full MatrixSpace of 0 by 0 dense matrices over Integer Ring

        There is a fair amount of latitude in the value of the ``right``
        vector, and the matrix that results can have entries from a new
        ring large enough to contain the result. If you know better,
        you can sometimes bring the result down to a less general ring.  ::

            sage: R.<t> = ZZ[]
            sage: v = vector(R, [12, 24*t])
            sage: w = vector(QQ, [1/2, 1/3, 1/4])
            sage: op = v.outer_product(w); op
            [   6    4    3]
            [12*t  8*t  6*t]
            sage: op.base_ring()
            Univariate Polynomial Ring in t over Rational Field
            sage: m = op.change_ring(R); m
            [   6    4    3]
            [12*t  8*t  6*t]
            sage: m.base_ring()
            Univariate Polynomial Ring in t over Integer Ring

        But some inputs are not compatible, even if vectors. ::

            sage: w = vector(GF(5), [1,2])
            sage: v = vector(GF(7), [1,2,3,4])
            sage: z = w.outer_product(v)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
            'Full MatrixSpace of 2 by 1 dense matrices over Finite Field of size 5' and
            'Full MatrixSpace of 1 by 4 dense matrices over Finite Field of size 7'

        And some inputs don't make any sense at all. ::

            sage: w = vector(QQ, [5,10])
            sage: z = w.outer_product(6)
            Traceback (most recent call last):
            ...
            TypeError: right operand in an outer product must be a vector,
            not an element of Integer Ring
        """
        if not isinstance(right, FreeModuleElement):
            raise TypeError('right operand in an outer product must be a vector, not an element of %s' % right.parent())
        return self.column()*right.row()

    # tensor product is an alias in the special case of two vectors
    tensor_product = outer_product

    def hermitian_inner_product(self, right):
        r"""
        Return the dot product, but with the entries of the first vector
        conjugated beforehand.

        INPUT:

        - ``right`` -- a vector of the same degree as ``self``

        OUTPUT:

        If ``self`` and ``right`` are the vectors `\vec{x}` and
        `\vec{y}` of degree `n` then this routine computes

        .. MATH::

            \sum_{i=1}^{n}\overline{x}_i{y}_i

        where the bar indicates complex conjugation.

        .. NOTE::

            If your vectors do not contain complex entries, then
            :meth:`dot_product` will return the same result without
            the overhead of conjugating elements of ``self``.

            If you are not computing a weighted inner product, and
            your vectors do not have complex entries, then the
            :meth:`dot_product` will return the same result.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: v = vector(CDF, [2+3*I, 5-4*I])
            sage: w = vector(CDF, [6-4*I, 2+3*I])
            sage: v.hermitian_inner_product(w)
            -2.0 - 3.0*I

        Sage implements a few specialized fields over the complex numbers,
        such as cyclotomic fields and quadratic number fields.  So long as
        the base rings have a conjugate method, then the Hermitian inner
        product will be available. ::

            sage: # needs sage.rings.number_field
            sage: Q.<a> = QuadraticField(-7)
            sage: a^2
            -7
            sage: v = vector(Q, [3+a, 5-2*a])
            sage: w = vector(Q, [6, 4+3*a])
            sage: v.hermitian_inner_product(w)
            17*a - 4

        The Hermitian inner product should be additive in
        each argument (we only need to test one), linear
        in each argument (with conjugation on the first scalar),
        and anti-commutative. ::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: alpha = CDF(5.0 + 3.0*I)
            sage: u = vector(CDF, [2+4*I, -3+5*I, 2-7*I])
            sage: v = vector(CDF, [-1+3*I, 5+4*I, 9-2*I])
            sage: w = vector(CDF, [8+3*I, -4+7*I, 3-6*I])
            sage: (u+v).hermitian_inner_product(w) == u.hermitian_inner_product(w) + v.hermitian_inner_product(w)
            True
            sage: (alpha*u).hermitian_inner_product(w) == alpha.conjugate()*u.hermitian_inner_product(w)
            True
            sage: u.hermitian_inner_product(alpha*w) == alpha*u.hermitian_inner_product(w)
            True
            sage: u.hermitian_inner_product(v) == v.hermitian_inner_product(u).conjugate()
            True

        For vectors with complex entries, the Hermitian inner product
        has a more natural relationship with the 2-norm (which is the
        default for the :meth:`norm` method). The norm squared equals
        the Hermitian inner product of the vector with itself.  ::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: v = vector(CDF, [-0.66+0.47*I, -0.60+0.91*I, -0.62-0.87*I, 0.53+0.32*I])
            sage: abs(v.norm()^2 - v.hermitian_inner_product(v)) < 1.0e-10
            True

        TESTS:

        This method is built on the :meth:`dot_product` method,
        which allows for a wide variety of inputs.  Any error
        handling happens there. ::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: v = vector(CDF, [2+3*I])
            sage: w = vector(CDF, [5+2*I, 3+9*I])
            sage: v.hermitian_inner_product(w)
            Traceback (most recent call last):
            ...
            ArithmeticError: degrees (1 and 2) must be the same
        """
        return (self.conjugate()).dot_product(right)

    def is_dense(self):
        """
        Return ``True`` if this is a dense vector, which is just a
        statement about the data structure, not the number of nonzero
        entries.

        EXAMPLES::

            sage: vector([1/2, 2/5, 0]).is_dense()
            True
            sage: vector([1/2, 2/5, 0], sparse=True).is_dense()
            False
        """
        return self.is_dense_c()

    cdef bint is_dense_c(self) noexcept:
        return self.parent().is_dense()

    def is_sparse(self):
        """
        Return ``True`` if this is a sparse vector, which is just a
        statement about the data structure, not the number of nonzero
        entries.

        EXAMPLES::

            sage: vector([1/2, 2/5, 0]).is_sparse()
            False
            sage: vector([1/2, 2/5, 0], sparse=True).is_sparse()
            True
        """
        return self.is_sparse_c()

    cdef bint is_sparse_c(self) noexcept:
        return self.parent().is_sparse()

    def is_vector(self):
        """
        Return ``True``, since this is a vector.

        EXAMPLES::

            sage: vector([1/2, 2/5, 0]).is_vector()
            True
        """
        return True

    def _macaulay2_(self, macaulay2=None):
        r"""
        Convert this vector to a Macaulay2 vector.

        EXAMPLES::

            sage: vector(QQ, [1, 2, 3])._macaulay2_()   # optional - macaulay2
            | 1 |
            | 2 |
            | 3 |
            sage: _.ring()                              # optional - macaulay2
            QQ

        ::

            sage: R.<x,y> = QQ[]
            sage: macaulay2(vector(R, [1, x+y]))        # optional - macaulay2
            |  1  |
            | x+y |

        TESTS:

        Entries of the vector get promoted to the base ring::

            sage: R.<x,y> = QQ[]
            sage: v = macaulay2(vector(R, [1, 2]))      # optional - macaulay2
            sage: v.ring()._operator('===', R).sage()   # optional - macaulay2
            True
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default
        return (macaulay2(self.base_ring()).matrix([self.list()]).transpose()
                .vector())

    def _mathematica_init_(self):
        """
        Return string representation of this vector as a Mathematica list.

        EXAMPLES::

            sage: # optional - mathematica, needs sage.symbolic
            sage: vector((1,2,3), QQ)._mathematica_init_()
            '{1/1, 2/1, 3/1}'
            sage: mathematica(vector((1,2,3), QQ))
            {1, 2, 3}
            sage: a = vector(SR, 5, [1, x, x^2, sin(x), pi]); a
            (1, x, x^2, sin(x), pi)
            sage: a._mathematica_init_()
            '{1, x, (x)^(2), Sin[x], Pi}'
        """
        return '{' + ', '.join(x._mathematica_init_() for x in self.list()) + '}'

    def _sympy_(self):
        """
        Return a SymPy column vector (matrix) corresponding to ``self``.

        OUTPUT:

        - An instance of either an ``ImmutableMatrix`` or ``ImmutableSparseMatrix``,
          regardless of whether ``self`` is mutable or not.

        EXAMPLES::

            sage: v = vector([1, 2, 3]); v
            (1, 2, 3)
            sage: sv = v._sympy_(); sv                                                  # needs sympy
            Matrix([
            [1],
            [2],
            [3]])
            sage: type(sv)                                                              # needs sympy
            <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>

            sage: w = vector({1: 1, 5: -1}, sparse=True)
            sage: sw = w._sympy_(); sw                                                  # needs sympy
            Matrix([
            [ 0],
            [ 1],
            [ 0],
            [ 0],
            [ 0],
            [-1]])
            sage: type(sw)                                                              # needs sympy
            <class 'sympy.matrices.immutable.ImmutableSparseMatrix'>

        If ``self`` was immutable, then converting the result to Sage gives
        back ``self``::

            sage: immv = vector([1, 2, 3], immutable=True)
            sage: immv._sympy_()._sage_() is immv                                       # needs sympy
            True

        If ``self`` was mutable, then converting back to Sage creates a new
        matrix (column vector)::

            sage: sv._sage_()                                                           # needs sympy
            [1]
            [2]
            [3]
            sage: sv._sage_() is v                                                      # needs sympy
            False
            sage: sv._sage_() == v                                                      # needs sympy
            False
        """
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        from sympy.matrices import ImmutableMatrix, ImmutableSparseMatrix
        if self.is_sparse():
            matrix = ImmutableSparseMatrix(self._degree, 1,
                                           {(i, 0): v
                                            for i, v in self.dict(copy=False).items()})
        else:
            matrix = ImmutableMatrix(self._degree, 1,
                                     self.list(copy=False))
        if self.is_immutable():
            matrix._sage_object = self
        return matrix

    def nonzero_positions(self):
        """
        Return the sorted list of integers ``i`` such that ``self[i] != 0``.

        EXAMPLES::

            sage: vector([-1,0,3,0,0,0,0.01]).nonzero_positions()
            [0, 2, 6]
        """
        v = self.list()
        cdef Py_ssize_t i
        return [i for i in range(self._degree) if v[i]]

    def support(self):   # do not override.
        """
        Return the integers ``i`` such that ``self[i] != 0``.
        This is the same as the ``nonzero_positions`` function.

        EXAMPLES::

            sage: vector([-1,0,3,0,0,0,0.01]).support()
            [0, 2, 6]
        """
        return self.nonzero_positions()

    cpdef int hamming_weight(self) noexcept:
        """
        Return the number of positions ``i`` such that ``self[i] != 0``.

        EXAMPLES::

            sage: vector([-1,0,3,0,0,0,0.01]).hamming_weight()
            3
        """
        cdef Py_ssize_t res = 0
        for x in iter(self.list()):
            if not x.is_zero():
                res += 1
        return res

    def _latex_(self):
        r"""
        Return a latex representation of the vector ``self``.

        OUTPUT:

        If ``self`` is the free module element (1,2,3,4),
        then a string with the following latex is returned:
        "\left(1,\,2,\,3,\,4\right)" (without the quotes).
        The vector is enclosed in parentheses by default,
        but the delimiters can be changed using the command
        ``latex.vector_delimiters(...)`` as in the example below.

        EXAMPLES::

            sage: v = vector(QQ, [1,2,3])
            sage: latex(v)
            \left(1,\,2,\,3\right)

        This is an example of how to change the delimiters.
        You have the power to mix and match, though it is
        probably not advisable.  For more detail see
        :meth:`~sage.misc.latex.Latex.vector_delimiters`.

            sage: latex.vector_delimiters('[', '\\rangle')
            sage: w = vector(CDF, [1,2,3])
            sage: latex(w)
            \left[1.0,\,2.0,\,3.0\right\rangle
        """
        from sage.misc.latex import latex
        vector_delimiters = latex.vector_delimiters()
        s = '\\left' + vector_delimiters[0]
        s += r',\,'.join(latex(a) for a in self.list())
        return s + '\\right' + vector_delimiters[1]

    def dense_vector(self):
        """
        Return dense version of ``self``.  If ``self`` is dense, just return
        ``self``; otherwise, create and return correspond dense vector.

        EXAMPLES::

            sage: vector([-1,0,3,0,0,0]).dense_vector().is_dense()
            True
            sage: vector([-1,0,3,0,0,0],sparse=True).dense_vector().is_dense()
            True
            sage: vector([-1,0,3,0,0,0],sparse=True).dense_vector()
            (-1, 0, 3, 0, 0, 0)
        """
        if self.is_dense():
            return self
        else:
            return self.parent().ambient_module().dense_module()(self.list())

    def sparse_vector(self):
        """
        Return sparse version of ``self``.  If ``self`` is sparse, just return
        self; otherwise, create and return correspond sparse vector.

        EXAMPLES::

            sage: vector([-1,0,3,0,0,0]).sparse_vector().is_sparse()
            True
            sage: vector([-1,0,3,0,0,0]).sparse_vector().is_sparse()
            True
            sage: vector([-1,0,3,0,0,0]).sparse_vector()
            (-1, 0, 3, 0, 0, 0)
        """
        if self.is_sparse():
            return self
        else:
            return self.parent().ambient_module().sparse_module()(self.list())

    def apply_map(self, phi, R=None, sparse=None):
        """
        Apply the given map phi (an arbitrary Python function or callable
        object) to this free module element. If R is not given,
        automatically determine the base ring of the resulting element.

        INPUT:

        - ``sparse`` -- boolean; will control whether the result is sparse.
          By default, the result is sparse iff self is sparse.

        - ``phi`` -- arbitrary Python function or callable object

        - ``R`` -- (optional) ring

        OUTPUT: a free module element over R

        EXAMPLES::

            sage: m = vector([1,x,sin(x+1)])                                            # needs sage.symbolic
            sage: m.apply_map(lambda x: x^2)                                            # needs sage.symbolic
            (1, x^2, sin(x + 1)^2)
            sage: m.apply_map(sin)                                                      # needs sage.symbolic
            (sin(1), sin(x), sin(sin(x + 1)))

        ::

            sage: m = vector(ZZ, 9, range(9))
            sage: k.<a> = GF(9)                                                         # needs sage.rings.finite_rings
            sage: m.apply_map(k)                                                        # needs sage.rings.finite_rings
            (0, 1, 2, 0, 1, 2, 0, 1, 2)

        In this example, we explicitly specify the codomain.

        ::

            sage: s = GF(3)
            sage: f = lambda x: s(x)
            sage: n = m.apply_map(f, k); n                                              # needs sage.rings.finite_rings
            (0, 1, 2, 0, 1, 2, 0, 1, 2)
            sage: n.parent()                                                            # needs sage.rings.finite_rings
            Vector space of dimension 9 over Finite Field in a of size 3^2

        If your map sends 0 to a nonzero value, then your resulting
        vector is not mathematically sparse::

            sage: v = vector([0] * 6 + [1], sparse=True); v
            (0, 0, 0, 0, 0, 0, 1)
            sage: v2 = v.apply_map(lambda x: x+1); v2
            (1, 1, 1, 1, 1, 1, 2)

        but it's still represented with a sparse data type::

            sage: parent(v2)
            Ambient sparse free module of rank 7 over the principal ideal domain Integer Ring

        This data type is inefficient for dense vectors, so you may
        want to specify sparse=False::

            sage: v2 = v.apply_map(lambda x: x+1, sparse=False); v2
            (1, 1, 1, 1, 1, 1, 2)
            sage: parent(v2)
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        Or if you have a map that will result in mostly zeroes, you may
        want to specify sparse=True::

            sage: v = vector(srange(10))
            sage: v2 = v.apply_map(lambda x: 0 if x else 1, sparse=True); v2
            (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: parent(v2)
            Ambient sparse free module of rank 10 over the principal ideal domain Integer Ring

        TESTS::

            sage: m = vector(SR,[])                                                     # needs sage.symbolic
            sage: m.apply_map(lambda x: x*x) == m                                       # needs sage.symbolic
            True

        Check that we don't unnecessarily apply phi to 0 in the sparse case::

            sage: m = vector(ZZ, range(1, 4), sparse=True)
            sage: m.apply_map(lambda x: 1/x)
            (1, 1/2, 1/3)

            sage: parent(vector(RDF, (), sparse=True).apply_map(lambda x: x, sparse=True))
            Sparse vector space of dimension 0 over Real Double Field
            sage: parent(vector(RDF, (), sparse=True).apply_map(lambda x: x, sparse=False))
            Vector space of dimension 0 over Real Double Field
            sage: parent(vector(RDF, (), sparse=False).apply_map(lambda x: x, sparse=True))
            Sparse vector space of dimension 0 over Real Double Field
            sage: parent(vector(RDF, (), sparse=False).apply_map(lambda x: x, sparse=False))
            Vector space of dimension 0 over Real Double Field

        Check that the bug in :issue:`14558` has been fixed::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(9)
            sage: v = vector([a, 0, 0, 0], sparse=True)
            sage: f = F.hom([a**3])
            sage: v.apply_map(f)
            (2*a + 1, 0, 0, 0)
        """
        if sparse is None:
            sparse = self.is_sparse()

        if self._degree == 0:
            if sparse == self.is_sparse():
                from copy import copy
                return copy(self)
            elif sparse:
                return self.sparse_vector()
            else:
                return self.dense_vector()

        v = None

        if self.is_sparse():
            zero_res = 0
            if len(self.dict(copy=False)) < self._degree:
                # OK, we have some zero entries.
                zero_res = phi(self.base_ring()(0))
                if not zero_res.is_zero():
                    # And phi maps 0 to a nonzero value.
                    v = [zero_res] * self._degree
                    for i,z in self.dict(copy=False).items():
                        v[i] = phi(z)

            if v is None:
                # phi maps 0 to 0 (or else we don't have any zeroes at all)
                v = dict([(i,phi(z)) for i,z in self.dict(copy=False).items()])
                # add a zero at the last position, if it is not already set.
                # This will help the constructor to determine the right degree.
                v.setdefault(self._degree-1, zero_res)
        else:
            v = [phi(z) for z in self.list()]

        if R is None:
            return vector(v, sparse=sparse)
        else:
            return vector(R, v, sparse=sparse)

    def _derivative(self, var=None):
        """
        Differentiate with respect to var by differentiating each element
        with respect to var.

        .. seealso:

           :meth:`derivative`

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: v = vector([1,x,x^2])
            sage: v._derivative(x)
            (0, 1, 2*x)
            sage: type(v._derivative(x)) == type(v)
            True
            sage: v = vector([1,x,x^2], sparse=True)
            sage: v._derivative(x)
            (0, 1, 2*x)
            sage: type(v._derivative(x)) == type(v)
            True

        If no variables are specified and the vector contains callable
        symbolic expressions, then calculate the matrix derivative
        (i.e., the Jacobian matrix)::

            sage: # needs sage.symbolic
            sage: T(r,theta) = [r*cos(theta), r*sin(theta)]
            sage: T
            (r, theta) |--> (r*cos(theta), r*sin(theta))
            sage: T.diff() # matrix derivative
            [   (r, theta) |--> cos(theta) (r, theta) |--> -r*sin(theta)]
            [   (r, theta) |--> sin(theta)  (r, theta) |--> r*cos(theta)]
            sage: diff(T) # matrix derivative again
            [   (r, theta) |--> cos(theta) (r, theta) |--> -r*sin(theta)]
            [   (r, theta) |--> sin(theta)  (r, theta) |--> r*cos(theta)]
            sage: T.diff().det() # Jacobian
            (r, theta) |--> r*cos(theta)^2 + r*sin(theta)^2
        """
        if var is None:
            if isinstance(self.coordinate_ring(), sage.rings.abc.CallableSymbolicExpressionRing):
                from sage.calculus.all import jacobian
                return jacobian(self, self.coordinate_ring().arguments())
            else:
                raise ValueError("No differentiation variable specified.")

        return self.apply_map(lambda x: x.derivative(var))

    def derivative(self, *args):
        """
        Derivative with respect to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        :meth:`diff` is an alias of this function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: v = vector([1,x,x^2])
            sage: v.derivative(x)
            (0, 1, 2*x)
            sage: type(v.derivative(x)) == type(v)
            True
            sage: v = vector([1,x,x^2], sparse=True)
            sage: v.derivative(x)
            (0, 1, 2*x)
            sage: type(v.derivative(x)) == type(v)
            True
            sage: v.derivative(x,x)
            (0, 0, 2)
        """
        from sage.misc.derivative import multi_derivative
        return multi_derivative(self, args)

    diff = derivative

    def integral(self, *args, **kwds):
        """
        Return a symbolic integral of the vector, component-wise.

        :meth:`integrate` is an alias of the function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: t = var('t')
            sage: r = vector([t,t^2,sin(t)])
            sage: r.integral(t)
            (1/2*t^2, 1/3*t^3, -cos(t))
            sage: integrate(r, t)
            (1/2*t^2, 1/3*t^3, -cos(t))
            sage: r.integrate(t, 0, 1)
            (1/2, 1/3, -cos(1) + 1)
        """
        from sage.misc.functional import integral
        return self.apply_map(lambda x: integral(x,*args, **kwds))

    integrate = integral

    def nintegral(self, *args, **kwds):
        """
        Return a numeric integral of the vector, component-wise, and
        the result of the nintegral command on each component of the
        input.

        :meth:`nintegrate` is an alias of the function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: t = var('t')
            sage: r = vector([t,t^2,sin(t)])
            sage: vec, answers = r.nintegral(t,0,1)
            sage: vec  # abs tol 1e-15
            (0.5, 0.3333333333333334, 0.4596976941318602)
            sage: type(vec)
            <class 'sage.modules.vector_real_double_dense.Vector_real_double_dense'>
            sage: answers
            [(0.5, 5.55111512312578...e-15, 21, 0),
             (0.3333333333333..., 3.70074341541719...e-15, 21, 0),
             (0.45969769413186..., 5.10366964392284...e-15, 21, 0)]

            sage: # needs sage.symbolic
            sage: r = vector([t,0,1], sparse=True)
            sage: r.nintegral(t, 0, 1)
            ((0.5, 0.0, 1.0),
             {0: (0.5, 5.55111512312578...e-15, 21, 0),
              2: (1.0, 1.11022302462515...e-14, 21, 0)})
        """
        # If Cython supported lambda functions, we would just do
        # return self.apply_map(lambda x: x.nintegral(*args, **kwds) for x in self)

        if self.is_sparse():
            v = [(i,z.nintegral(*args,**kwds)) for i,z in self.dict(copy=False).items()]
            answers = dict([(i,a[0]) for i,a in v])
            v=dict(v)
        else:
            v = [z.nintegral(*args,**kwds) for z in self.list()]
            answers = [a[0] for a in v]

        return (vector(answers,sparse=self.is_sparse()), v)

    nintegrate = nintegral

    def concatenate(self, other, *, ring=None):
        r"""
        Return the result of concatenating this vector with a sequence
        of elements given by another iterable.

        If the optional keyword argument ``ring`` is passed, this method
        will return a vector over the specified ring (or fail). If no
        base ring is given, the base ring is determined automatically by
        the :func:`vector` constructor.

        EXAMPLES::

            sage: v = vector([1, 2, 3])
            sage: w = vector([4, 5])
            sage: v.concatenate(w)
            (1, 2, 3, 4, 5)
            sage: v.parent()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: w.parent()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: v.concatenate(w).parent()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring

        Forcing a base ring is possible using the ``ring`` argument::

            sage: v.concatenate(w, ring=QQ)
            (1, 2, 3, 4, 5)
            sage: v.concatenate(w, ring=QQ).parent()
            Vector space of dimension 5 over Rational Field

        ::

            sage: v.concatenate(w, ring=Zmod(3))
            (1, 2, 0, 1, 2)

        The method accepts arbitrary iterables of elements which can
        be coerced to a common base ring::

            sage: v.concatenate(range(4,8))
            (1, 2, 3, 4, 5, 6, 7)
            sage: v.concatenate(range(4,8)).parent()
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        ::

            sage: w2 = [4, QQbar(-5).sqrt()]
            sage: v.concatenate(w2)
            (1, 2, 3, 4, 2.236...*I)
            sage: v.concatenate(w2).parent()
            Vector space of dimension 5 over Algebraic Field
            sage: w2 = vector(w2)
            sage: v.concatenate(w2)
            (1, 2, 3, 4, 2.236...*I)
            sage: v.concatenate(w2).parent()
            Vector space of dimension 5 over Algebraic Field

        ::

            sage: w2 = polygen(QQ)^4 + 5
            sage: v.concatenate(w2)
            (1, 2, 3, 5, 0, 0, 0, 1)
            sage: v.concatenate(w2).parent()
            Vector space of dimension 8 over Rational Field
            sage: v.concatenate(w2, ring=ZZ)
            (1, 2, 3, 5, 0, 0, 0, 1)
            sage: v.concatenate(w2, ring=ZZ).parent()
            Ambient free module of rank 8 over the principal ideal domain Integer Ring

        ::

            sage: v.concatenate(GF(9).gens())
            (1, 2, 0, z2)
            sage: v.concatenate(GF(9).gens()).parent()
            Vector space of dimension 4 over Finite Field in z2 of size 3^2
        """
        from itertools import chain
        coeffs = chain(self, other)
        if ring is not None:
            return vector(ring, coeffs)
        return vector(coeffs)


#############################################
# Generic dense element
#############################################

@cython.binding(True)
def make_FreeModuleElement_generic_dense(parent, entries, degree):
    """
    EXAMPLES::

        sage: sage.modules.free_module_element.make_FreeModuleElement_generic_dense(QQ^3, [1,2,-3/7], 3)
        (1, 2, -3/7)
    """
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef FreeModuleElement_generic_dense v
    v = FreeModuleElement_generic_dense.__new__(FreeModuleElement_generic_dense)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    return v


@cython.binding(True)
def make_FreeModuleElement_generic_dense_v1(parent, entries, degree, is_mutable):
    """
    EXAMPLES::

        sage: v = sage.modules.free_module_element.make_FreeModuleElement_generic_dense_v1(QQ^3, [1,2,-3/7], 3, True); v
        (1, 2, -3/7)
        sage: v[0] = 10; v
        (10, 2, -3/7)
        sage: v = sage.modules.free_module_element.make_FreeModuleElement_generic_dense_v1(QQ^3, [1,2,-3/7], 3, False); v
        (1, 2, -3/7)
        sage: v[0] = 10
        Traceback (most recent call last):
        ...
        ValueError: vector is immutable; please change a copy instead (use copy())
    """
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v2
    # and changed the reduce method below.
    cdef FreeModuleElement_generic_dense v
    v = FreeModuleElement_generic_dense.__new__(FreeModuleElement_generic_dense)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    v._is_immutable = not is_mutable
    return v


cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    """
    A generic dense element of a free module.

    TESTS::

        sage: V = ZZ^3
        sage: loads(dumps(V)) == V
        True
        sage: v = V.0
        sage: loads(dumps(v)) == v
        True
        sage: v = (QQ['x']^3).0
        sage: loads(dumps(v)) == v
        True

    ::

        sage: v = vector([1,2/3,pi])                                                    # needs sage.symbolic
        sage: v == v
        True

    ::

        sage: v = vector(RR, [1,2/3,pi])                                                # needs sage.symbolic
        sage: v.set_immutable()
        sage: isinstance(hash(v), int)
        True
    """
    cdef _new_c(self, object v):
        """
        Create a new dense free module element with minimal overhead and
        no type checking.

        INPUT:

        - ``v`` -- list which is used as the new entries (without copying)
        """
        cdef type t = type(self)
        cdef FreeModuleElement_generic_dense x = t.__new__(t)
        x._is_immutable = 0
        x._parent = self._parent
        x._entries = v
        x._degree = self._degree
        return x

    cdef bint is_dense_c(self) noexcept:
        return 1

    cdef bint is_sparse_c(self) noexcept:
        return 0

    def __copy__(self):
        """
        Return a copy of this generic dense vector.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: v = vector([-1,0,3,pi])
            sage: type(v)
            <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
            sage: v.__copy__()
            (-1, 0, 3, pi)
            sage: v.__copy__() is v
            False

            sage: copy(v)                                                               # needs sage.symbolic
            (-1, 0, 3, pi)
            sage: copy(v) == v                                                          # needs sage.symbolic
            True
            sage: copy(v) is v                                                          # needs sage.symbolic
            False
        """
        return self._new_c(list(self._entries))

    def __init__(self, parent, entries, coerce=True, copy=True):
        """
        EXAMPLES::

            sage: type(vector(RR, [-1,0,2/3,pi,oo]))                                    # needs sage.symbolic
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_dense'>

        We can initialize with lists, tuples and derived types::

            sage: from sage.modules.free_module_element import FreeModuleElement_generic_dense
            sage: FreeModuleElement_generic_dense(RR^5, [-1,0,2/3,pi,oo])               # needs sage.symbolic
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: FreeModuleElement_generic_dense(RR^5, (-1,0,2/3,pi,oo))               # needs sage.symbolic
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: FreeModuleElement_generic_dense(RR^5, Sequence([-1,0,2/3,pi,oo]))     # needs sage.symbolic
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: FreeModuleElement_generic_dense(RR^0, 0)
            ()

        TESTS:

        Disabling coercion can lead to illegal objects::

            sage: FreeModuleElement_generic_dense(RR^5, [-1,0,2/3,pi,oo], coerce=False)             # needs sage.symbolic
            (-1, 0, 2/3, pi, +Infinity)

        We test the ``copy`` flag::

            sage: # needs sage.symbolic
            sage: from sage.modules.free_module_element import FreeModuleElement_generic_dense
            sage: L = [RR(x) for x in (-1,0,2/3,pi,oo)]
            sage: FreeModuleElement_generic_dense(RR^5, tuple(L), coerce=False, copy=False)
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: v = FreeModuleElement_generic_dense(RR^5, L, coerce=False, copy=False)
            sage: L[4] = 42.0
            sage: v  # last entry changed since we didn't copy
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, 42.0000000000000)

        ::

            sage: # needs sage.symbolic
            sage: L = [RR(x) for x in (-1,0,2/3,pi,oo)]
            sage: v = FreeModuleElement_generic_dense(RR^5, L, coerce=False, copy=True)
            sage: L[4] = 42.0
            sage: v  # last entry did not change
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)

        Check that :issue:`11751` is fixed::

            sage: K.<x> = QQ[]
            sage: M = K^1
            sage: N = M.span([[1/x]]); N
            Free module of degree 1 and rank 1 over Univariate Polynomial Ring in x over Rational Field
            Echelon basis matrix:
            [1/x]
            sage: N([1/x]) # this used to fail prior to #11751
            (1/x)
            sage: N([1/x^2])
            Traceback (most recent call last):
            ...
            TypeError: element [1/x^2] is not in free module

        ::

            sage: L=K^2
            sage: R=L.span([[x,0],[0,1/x]], check=False, already_echelonized=True)
            sage: R.basis()[0][0].parent()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
            sage: R=L.span([[x,x^2]])
            sage: R.basis()[0][0].parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        FreeModuleElement.__init__(self, parent)
        R = self.base_ring()
        if not entries:
            entries = [R.zero()]*self._degree
        else:
            if type(entries) is not list:
                if not isinstance(entries, (list, tuple)):
                    raise TypeError("entries must be a list or tuple, not %s" % type(entries))
                copy = True  # ensure we have a true Python list

            if len(entries) != self._degree:
                raise TypeError("entries must be a list of length %s" % self.degree())
            if coerce:
                coefficient_ring = parent.coordinate_ring()
                try:
                    entries = [coefficient_ring(x) for x in entries]
                except TypeError:
                    raise TypeError("Unable to coerce entries (=%s) to coefficients in %s" % (entries, coefficient_ring))
            elif copy:
                entries = list(entries)  # make a copy/convert to list
        self._entries = entries

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _add_(left, right):
        """
        Add left and right.

        EXAMPLES::

            sage: v = vector([1,2/3,pi]); w = vector([-2/3,pi^2,1])                     # needs sage.symbolic
            sage: v._add_(w)                                                            # needs sage.symbolic
            (1/3, pi^2 + 2/3, pi + 1)
        """
        cdef list a = left._entries
        cdef list b = (<FreeModuleElement_generic_dense>right)._entries
        v = [(<RingElement> a[i])._add_(<RingElement> b[i]) for i in range(left._degree)]
        return left._new_c(v)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _sub_(left, right):
        """
        Subtract right from left.

        EXAMPLES::

            sage: V = QQ^5
            sage: W = V.span([V.1, V.2])
            sage: W.0 - V.0
            (-1, 1, 0, 0, 0)
            sage: V.0 - W.0
            (1, -1, 0, 0, 0)
        """
        cdef list a = left._entries
        cdef list b = (<FreeModuleElement_generic_dense>right)._entries
        v = [(<RingElement> a[i])._sub_(<RingElement> b[i]) for i in range(left._degree)]
        return left._new_c(v)

    cpdef _rmul_(self, Element left):
        """
        EXAMPLES::

            sage: V = ZZ['x']^5
            sage: 5 * V.0
            (5, 0, 0, 0, 0)
        """
        if left._parent is self._parent._base:
            v = [left._mul_(<RingElement>x) for x in self._entries]
        else:
            v = [left * x for x in self._entries]
        return self._new_c(v)

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES::

            sage: v = vector([-1,0,3,pi])                                               # needs sage.symbolic
            sage: v._lmul_(2/3)                                                         # needs sage.symbolic
            (-2/3, 0, 2, 2/3*pi)
            sage: v * (2/3)                                                             # needs sage.symbolic
            (-2/3, 0, 2, 2/3*pi)
        """
        if right._parent is self._parent._base:
            v = [(<RingElement>x)._mul_(right) for x in self._entries]
        else:
            v = [x * right for x in self._entries]
        return self._new_c(v)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _pairwise_product_(left, Vector right):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = vector([x,x^2,3*x]); w = vector([2*x,x,3+x])
            sage: v.pairwise_product(w)
            (2*x^2, x^3, 3*x^2 + 9*x)
            sage: w.pairwise_product(v)
            (2*x^2, x^3, 3*x^2 + 9*x)
        """
        if right._parent is not left._parent:
            right = left.parent().ambient_module()(right)
        cdef list a = left._entries
        cdef list b = (<FreeModuleElement_generic_dense>right)._entries
        v = [(<RingElement> a[i])._mul_(<RingElement> b[i]) for i in range(left._degree)]
        return left._new_c(v)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: v = vector([-1,0,3,pi])                                               # needs sage.symbolic
            sage: v.__reduce__()                                                        # needs sage.symbolic
            (<cyfunction make_FreeModuleElement_generic_dense_v1 at ...>,
             (Vector space of dimension 4 over Symbolic Ring, [-1, 0, 3, pi], 4, True))
        """
        return (make_FreeModuleElement_generic_dense_v1, (self._parent, self._entries,
                                                          self._degree, not self._is_immutable))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: v = vector(RR, [-1,0,2/3,pi])                                         # needs sage.symbolic
            sage: v.get(3)                                                              # needs sage.symbolic
            3.14159265358979

        ::

            sage: v = vector([RR(1), RR(2)]); v
            (1.00000000000000, 2.00000000000000)
            sage: v[0]
            1.00000000000000
            sage: v[-1]
            2.00000000000000
            sage: v[4]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
            sage: v[-4]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range

        ::

            sage: v = vector(QQ['x,y'], [1,2, 'x*y'])
            sage: v
            (1, 2, x*y)
            sage: v[1:]
            (2, x*y)
        """
        return self._entries[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: v = vector(RR, [-1, 0, 2/3, pi])                                      # needs sage.symbolic
            sage: v.set(3, RR(1))                                                       # needs sage.symbolic
            sage: v                                                                     # needs sage.symbolic
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 1.00000000000000)
        """
        self._entries[i] = value

    def list(self, copy=True):
        """
        Return list of elements of ``self``.

        INPUT:

        - ``copy`` -- boolean; return list of underlying entries

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: v = vector([x,y,z])
            sage: type(v)
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_dense'>
            sage: a = v.list(); a
            [x, y, z]
            sage: a[0] = x*y; v
            (x, y, z)
            sage: a = v.list(copy=False); a
            [x, y, z]
            sage: a[0] = x*y; v
            (x*y, y, z)
        """
        if copy:
            return list(self._entries)
        else:
            return self._entries

    def __call__(self, *args, **kwargs):
        """
        Calling a free module element returns the result of calling each
        component.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x, y = var('x,y')
            sage: f = x^2 + y^2
            sage: g = f.gradient()
            sage: g
            (2*x, 2*y)
            sage: type(g)
            <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
            sage: g(y=2, x=3)
            (6, 4)
            sage: f(x,y) = x^2 + y^2
            sage: g = f.gradient()
            sage: g(3,2)
            (6, 4)
            sage: g(x=3, y=2)
            (6, 4)
        """
        return vector([e(*args, **kwargs) for e in self])

    def function(self, *args):
        """
        Return a vector over a callable symbolic expression ring.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x, y = var('x,y')
            sage: v = vector([x, y, x*sin(y)])
            sage: w = v.function([x,y]); w
            (x, y) |--> (x, y, x*sin(y))
            sage: w.coordinate_ring()
            Callable function ring with arguments (x, y)
            sage: w(1,2)
            (1, 2, sin(2))
            sage: w(2,1)
            (2, 1, 2*sin(1))
            sage: w(y=1,x=2)
            (2, 1, 2*sin(1))

        ::

            sage: # needs sage.symbolic
            sage: x,y = var('x,y')
            sage: v = vector([x, y, x*sin(y)])
            sage: w = v.function([x]); w
            x |--> (x, y, x*sin(y))
            sage: w.coordinate_ring()
            Callable function ring with argument x
            sage: w(4)
            (4, y, 4*sin(y))
        """
        from sage.symbolic.callable import CallableSymbolicExpressionRing
        return vector(CallableSymbolicExpressionRing(args), self.list())


#############################################
# Generic sparse element
#############################################

@cython.binding(True)
def make_FreeModuleElement_generic_sparse(parent, entries, degree):
    """
    EXAMPLES::

        sage: v = sage.modules.free_module_element.make_FreeModuleElement_generic_sparse(QQ^3, {2:5/2}, 3); v
        (0, 0, 5/2)
    """
    cdef FreeModuleElement_generic_sparse v
    v = FreeModuleElement_generic_sparse.__new__(FreeModuleElement_generic_sparse)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    return v


@cython.binding(True)
def make_FreeModuleElement_generic_sparse_v1(parent, entries, degree, is_mutable):
    """
    EXAMPLES::

        sage: v = sage.modules.free_module_element.make_FreeModuleElement_generic_sparse_v1(QQ^3, {2:5/2}, 3, False); v
        (0, 0, 5/2)
        sage: v.is_mutable()
        False
    """
    cdef FreeModuleElement_generic_sparse v
    v = FreeModuleElement_generic_sparse.__new__(FreeModuleElement_generic_sparse)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    v._is_immutable = not is_mutable
    return v


cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    """
    A generic sparse free module element is a dictionary with keys ints
    i and entries in the base ring.

    TESTS::

        sage: v = vector([1,2/3,pi], sparse=True)                                       # needs sage.symbolic
        sage: v.set_immutable()                                                         # needs sage.symbolic
        sage: isinstance(hash(v), int)                                                  # needs sage.symbolic
        True

    Pickling works::

        sage: v = FreeModule(ZZ, 3, sparse=True).0
        sage: loads(dumps(v)) == v
        True
        sage: v = FreeModule(Integers(8)['x,y'], 5, sparse=True).1
        sage: loads(dumps(v)) - v
        (0, 0, 0, 0, 0)

    ::

        sage: a = vector([-1,0,1/1],sparse=True); b = vector([-1/1,0,0],sparse=True)
        sage: a.parent()
        Sparse vector space of dimension 3 over Rational Field
        sage: b - a
        (0, 0, -1)
        sage: (b-a).dict()
        {2: -1}
    """
    cdef _new_c(self, object v):
        """
        Create a new sparse free module element with minimal overhead and
        no type checking.

        INPUT:

        - ``v`` -- dictionary which is used as the new entries (without
          copying)
        """
        cdef type t = type(self)
        cdef FreeModuleElement_generic_sparse x = t.__new__(t)
        x._is_immutable = 0
        x._parent = self._parent
        x._entries = v
        x._degree = self._degree
        return x

    cdef bint is_dense_c(self) noexcept:
        return 0

    cdef bint is_sparse_c(self) noexcept:
        return 1

    def __copy__(self):
        """
        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: v.__copy__()                                                          # needs sage.symbolic
            (1, 2/3, pi)
        """
        return self._new_c(dict(self._entries))

    def __init__(self, parent, entries=0, coerce=True, copy=True):
        """
        EXAMPLES::

            sage: v = sage.modules.free_module_element.FreeModuleElement_generic_sparse(VectorSpace(QQ,3,sparse=True), {1:5/4}); v
            (0, 5/4, 0)
            sage: v.is_sparse()
            True

        We can initialize with dicts, lists, tuples and derived types::

            sage: from sage.modules.free_module_element import FreeModuleElement_generic_sparse
            sage: def S(R,n):
            ....:     return FreeModule(R, n, sparse=True)

            sage: # needs sage.symbolic
            sage: FreeModuleElement_generic_sparse(S(RR,5), {0:-1, 2:2/3, 3:pi, 4:oo})
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: FreeModuleElement_generic_sparse(S(RR,5), [-1,0,2/3,pi,oo])
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: FreeModuleElement_generic_sparse(S(RR,5), (-1,0,2/3,pi,oo))
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)
            sage: FreeModuleElement_generic_sparse(S(RR,5), Sequence([-1,0,2/3,pi,oo]))
            (-1.00000000000000, 0.000000000000000, 0.666666666666667, 3.14159265358979, +infinity)

            sage: FreeModuleElement_generic_sparse(S(RR,0), 0)
            ()

            sage: from collections import defaultdict
            sage: D = defaultdict(RR)
            sage: D[0] = -1
            sage: FreeModuleElement_generic_sparse(S(RR,5), D)
            (-1.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000)

        TESTS:

        Test that :issue:`11751` is fixed::

            sage: K.<x> = QQ[]
            sage: M = FreeModule(K, 1, sparse=True)
            sage: N = M.span([{0:1/x}]); N
            Sparse free module of degree 1 and rank 1 over Univariate Polynomial Ring in x over Rational Field
            Echelon basis matrix:
            [1/x]
            sage: N({0:1/x}) # this used to fail prior to #11751
            (1/x)
            sage: N({0:1/x^2})
            Traceback (most recent call last):
            ...
            TypeError: element {0: 1/x^2} is not in free module

        ::

            sage: L = FreeModule(K, 2, sparse=True)
            sage: R = L.span([{0:x, 1:0}, {0:0, 1:1/x}], check=False, already_echelonized=True)
            sage: R.basis()[0][0].parent()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
            sage: R = L.span([{0:x, 1:x^2}])
            sage: R.basis()[0][0].parent()
            Univariate Polynomial Ring in x over Rational Field

        Test that :issue:`17101` is fixed::

            sage: # needs sage.rings.real_interval_field
            sage: v = vector([RIF(-1, 1)], sparse=True)
            sage: v.is_zero()
            False

        We correctly initialize values which become 0 only after coercion::

            sage: v = FreeModuleElement_generic_sparse(S(GF(3), 6), [1,2,3,4,5,6])
            sage: v.nonzero_positions()
            [0, 1, 3, 4]
        """
        # WARNING: In creation, we do not check that the indices i satisfy
        #     0 <= i < degree
        # or even that the indices are integers.
        FreeModuleElement.__init__(self, parent)
        cdef Py_ssize_t i
        cdef dict entries_dict, e
        if not entries:
            entries_dict = {}
        else:
            if type(entries) is not dict:
                if isinstance(entries, dict):
                    # Convert derived type to dict
                    copy = True
                    entries_dict = <dict> entries
                elif isinstance(entries, (list, tuple)):
                    if len(entries) != self._degree:
                        raise TypeError("entries has the wrong length")
                    entries_dict = {}
                    for i in range(self._degree):
                        x = entries[i]
                        if x:
                            entries_dict[i] = x
                    copy = False
                else:
                    raise TypeError("entries must be a dict, list or tuple, not %s" % type(entries))
            else:
                entries_dict = <dict> entries
            if coerce:
                coefficient_ring = parent.coordinate_ring()
                e = entries_dict
                entries_dict = {}
                try:
                    for k, x in (<dict> e).iteritems():
                        x = coefficient_ring(x)
                        if x:
                            entries_dict[k] = x
                except TypeError:
                    raise TypeError("unable to coerce value (=%s) of entries dict (=%s) to %s" % (x, entries, coefficient_ring))
            elif copy:
                entries_dict = dict(entries_dict)  # make a copy/convert to dict
        self._entries = entries_dict

    cpdef _add_(left, right):
        """
        Add left and right.

        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: v._add_(v)                                                            # needs sage.symbolic
            (2, 4/3, 2*pi)
        """
        cdef dict v = dict((<FreeModuleElement_generic_sparse>right)._entries)
        for i, a in left._entries.iteritems():
            if i in v:
                sum = (<RingElement>a)._add_(<RingElement> v[i])
                if sum:
                    v[i] = sum
                else:
                    del v[i]
            elif a:
                v[i] = a
        return left._new_c(v)

    cpdef _sub_(left, right):
        """
        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: v._sub_(v)                                                            # needs sage.symbolic
            (0, 0, 0)
        """
        cdef dict v = dict(left._entries)   # dict to make a copy
        for i, a in (<FreeModuleElement_generic_sparse>right)._entries.iteritems():
            if i in v:
                diff = (<RingElement> v[i])._sub_(<RingElement>a)
                if diff:
                    v[i] = diff
                else:
                    del v[i]
            elif a:
                v[i] = -a
        return left._new_c(v)

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: v._lmul_(SR(3))                                                       # needs sage.symbolic
            (3, 2, 3*pi)
        """
        cdef dict v = {}
        if right:
            for i, a in self._entries.iteritems():
                prod = (<RingElement>a)._mul_(right)
                if prod:
                    v[i] = prod
        return self._new_c(v)

    cpdef _rmul_(self, Element left):
        """
        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: v._rmul_(SR(3))                                                       # needs sage.symbolic
            (3, 2, 3*pi)
        """
        cdef dict v = {}
        if left:
            for i, a in self._entries.iteritems():
                prod = left._mul_(a)
                if prod:
                    v[i] = prod
        return self._new_c(v)

    cpdef _dot_product_coerce_(left, Vector right):
        """
        Return the dot product of left and right.

        EXAMPLES::

            sage: v = vector([1,2,0], sparse=True); w = vector([0,5,-9], sparse=True)
            sage: v * w
            10
            sage: w * v
            10

        Over different rings::

            sage: R.<x> = ZZ[]
            sage: v = vector(RDF, [0,1,2], sparse=True)
            sage: w = vector(R, [x,0,0], sparse=True)
            sage: p = v._dot_product_coerce_(w)
            sage: p
            0
            sage: parent(p)
            Univariate Polynomial Ring in x over Real Double Field

        Zero-dimensional vectors also work correctly::

            sage: v = vector(RDF, [], sparse=True)
            sage: w = vector(R, [], sparse=True)
            sage: parent(v._dot_product_coerce_(w))
            Univariate Polynomial Ring in x over Real Double Field

        TESTS:

        Check that :issue:`19377` is fixed::

            sage: w = vector(ZZ, (1,2,3), sparse=False)
            sage: v = vector(ZZ, (1,2,3), sparse=True)
            sage: v._dot_product_coerce_(w)
            14
        """
        cdef dict e
        try:
            e = (<FreeModuleElement_generic_sparse?>right)._entries
        except TypeError:
            e = right.dict()
        z = left.base_ring().zero()
        if left.base_ring() is not right.base_ring():
            z *= right.base_ring().zero()
        for i, a in left._entries.iteritems():
            if i in e:
                z += a * e[i]
        return z

    cpdef _pairwise_product_(left, Vector right):
        """
        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True); w = vector([-2/3,pi^2,1],sparse=True)        # needs sage.symbolic
            sage: v._pairwise_product_(w)                                               # needs sage.symbolic
            (-2/3, 2/3*pi^2, pi)
        """
        # Component wise vector * vector multiplication.
        cdef dict e = (<FreeModuleElement_generic_sparse>right)._entries
        cdef dict v = {}
        for i, a in left._entries.iteritems():
            if i in e:
                prod = (<RingElement>a)._mul_(<RingElement> e[i])
                if prod:
                    v[i] = prod
        return left._new_c(v)

    cpdef _richcmp_(left, right, int op):
        """
        Compare two sparse free module elements.

        Free module elements are compared in lexicographic order on
        the underlying list of coefficients. Two free module elements
        are equal if their coefficients are the same. (This is true
        even if one is sparse and one is dense.)

        TESTS::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: w = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: w == v                                                                # needs sage.symbolic
            True

        Check that the bug in :issue:`13929` has been fixed::

            sage: V = FreeModule(GF(3), 2, sparse=True)
            sage: a = V([0,1])
            sage: b = V([1,0])
            sage: a < b
            True
        """
        a = sorted((<FreeModuleElement_generic_sparse>left)._entries.iteritems())
        b = sorted((<FreeModuleElement_generic_sparse>right)._entries.iteritems())

        return richcmp([(-x, y) for x, y in a], [(-x, y) for x, y in b], op)

    def items(self):
        """
        Return an iterator over the entries of ``self``.

        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: next(v.items())                                                       # needs sage.symbolic
            (0, 1)
            sage: list(v.items())                                                       # needs sage.symbolic
            [(0, 1), (1, 2/3), (2, pi)]

        TESTS:

        Using iteritems as an alias::

            sage: list(v.iteritems())                                                   # needs sage.symbolic
            [(0, 1), (1, 2/3), (2, pi)]
        """
        return iter(self._entries.iteritems())

    iteritems = items

    def __reduce__(self):
        """
        EXAMPLES::

            sage: v = vector([1,2/3,pi], sparse=True)                                   # needs sage.symbolic
            sage: v.__reduce__()                                                        # needs sage.symbolic
            (<cyfunction make_FreeModuleElement_generic_sparse_v1 at ...>,
             (Sparse vector space of dimension 3 over Symbolic Ring, {0: 1, 1: 2/3, 2: pi}, 3, True))
        """
        return (make_FreeModuleElement_generic_sparse_v1, (self._parent, self._entries,
                                                           self._degree, not self._is_immutable))

    @cython.cdivision(True)
    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: v = vector(RR, range(6), sparse=True); v
            (0.000000000000000, 1.00000000000000, 2.00000000000000, 3.00000000000000, 4.00000000000000, 5.00000000000000)
            sage: v[1]
            1.00000000000000
            sage: v[-1]
            5.00000000000000
            sage: v[9]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
            sage: v[-7]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
            sage: v[::2]
            (0.000000000000000, 2.00000000000000, 4.00000000000000)
            sage: v[5:2:-1]
            (5.00000000000000, 4.00000000000000, 3.00000000000000)

        All these operations with zero vectors should be very fast::

            sage: v = vector(RR, 10^9, sparse=True)
            sage: v[123456789]
            0.000000000000000
            sage: w = v[::-1]
            sage: v[::-250000000]
            (0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000)
            sage: v[123456789:123456798:3]
            (0.000000000000000, 0.000000000000000, 0.000000000000000)
        """
        cdef Py_ssize_t d = self._degree
        cdef Py_ssize_t start, stop, step, slicelength
        cdef Py_ssize_t min, max, mod
        cdef Py_ssize_t k, n
        cdef dict newentries
        if isinstance(i, slice):
            PySlice_GetIndicesEx(i, d, &start, &stop, &step, &slicelength)
            if step > 0:
                min = start
                max = stop-1
            else:
                min = stop+1
                max = start
            mod = start % step
            # Loop over the old dict and convert old index n to new
            # index k in slice
            newentries = {}
            for n, x in self._entries.iteritems():
                if min <= n <= max and n % step == mod:
                    k = (n - start) // step
                    newentries[k] = x
            from sage.modules.free_module import FreeModule
            M = FreeModule(self.coordinate_ring(), slicelength, sparse=True)
            return M(newentries, coerce=False, copy=False)

        n = i
        if n < 0:
            n += d
        if n < 0 or n >= d:
            raise IndexError("vector index out of range")
        return self.get_unsafe(n)

    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: v = vector([-1,0,2/3,pi], sparse=True)                                # needs sage.symbolic
            sage: v.get(1)                                                              # needs sage.symbolic
            0
            sage: v.get(2)                                                              # needs sage.symbolic
            2/3

        For this class, 0 is returned if the access is out of bounds::

            sage: v.get(10)                                                             # needs sage.symbolic
            0
        """
        try:
            return self._entries[i]
        except KeyError:
            return self.coordinate_ring().zero()

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: V = VectorSpace(GF(17), 10000000, sparse=True)
            sage: w = V(0)
            sage: w[39893] = 20
            sage: w[39893]
            3
            sage: w[39000:39003] = [4, 5, 6]; w[39000:39003]
            (4, 5, 6)
            sage: parent(w[39893])
            Finite Field of size 17
            sage: w[39893] = sqrt(2)                                                    # needs sage.rings.finite_rings sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: self must be a numeric expression

        ::

            sage: # needs sage.symbolic
            sage: v = vector([1,2/3,pi], sparse=True)
            sage: v.set(1, pi^3)
            sage: v
            (1, pi^3, pi)
            sage: v.set(2, SR(0))
            sage: v
            (1, pi^3, 0)

        This assignment is illegal::

            sage: v.set(10, pi)                                                         # needs sage.symbolic

        This lack of bounds checking causes trouble later::

            sage: v                                                                     # needs sage.symbolic
            <repr(<sage.modules.free_module.FreeModule_ambient_field_with_category.element_class at 0x...>) failed:
             IndexError: list assignment index out of range>
        """
        if value:
            self._entries[i] = value
        else:
            self._entries.pop(i, None)

    def denominator(self):
        """
        Return the least common multiple of the denominators of the
        entries of ``self``.

        EXAMPLES::

            sage: v = vector([1/2,2/5,3/14], sparse=True)
            sage: v.denominator()
            70
        """
        # It may be that the coordinates do not have a denominator
        # (but if one coordinate has it, they all should have it)
        d = self.coordinate_ring().one()
        try:
            d = d.denominator()
        except AttributeError:
            return d
        for y in self._entries.itervalues():
            d = d.lcm(y.denominator())
        return d

    def dict(self, copy=True):
        """
        Return dictionary of nonzero entries of ``self``.

        More precisely, this returns a dictionary whose keys are indices
        of basis elements in the support of ``self`` and whose values are
        the corresponding coefficients.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``self`` is internally
          represented by a dictionary ``d``, then make a copy of ``d``.
          If ``False``, then this can cause undesired behavior by
          mutating ``d``.

        OUTPUT: Python dictionary

        EXAMPLES::

            sage: v = vector([0,0,0,0,1/2,0,3/14], sparse=True)
            sage: v.dict()
            {4: 1/2, 6: 3/14}
            sage: sorted(v.support())
            [4, 6]
        """
        if copy:
            return dict(self._entries)
        else:
            return self._entries

    def list(self, copy=True):
        """
        Return list of elements of ``self``.

        INPUT:

        - ``copy`` -- ignored for sparse vectors

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: M = FreeModule(R, 3, sparse=True) * (1/x)
            sage: v = M([-x^2, 3/x, 0])
            sage: type(v)
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>
            sage: a = v.list()
            sage: a
            [-x^2, 3/x, 0]
            sage: [parent(c) for c in a]
            [Fraction Field of Univariate Polynomial Ring in x over Rational Field,
             Fraction Field of Univariate Polynomial Ring in x over Rational Field,
             Fraction Field of Univariate Polynomial Ring in x over Rational Field]
        """
        z = self._parent.coordinate_ring().zero()
        cdef list v = [z] * self._degree
        for i, a in self._entries.iteritems():
            v[i] = a
        return v

    def nonzero_positions(self):
        """
        Return the list of numbers ``i`` such that ``self[i] != 0``.

        EXAMPLES::

            sage: v = vector({1: 1, 3: -2})
            sage: w = vector({1: 4, 3: 2})
            sage: v+w
            (0, 5, 0, 0)
            sage: (v+w).nonzero_positions()
            [1]
        """
        return sorted(self._entries)

    cpdef int hamming_weight(self) noexcept:
        """
        Return the number of positions ``i`` such that ``self[i] != 0``.

        EXAMPLES::

            sage: v = vector({1: 1, 3: -2})
            sage: w = vector({1: 4, 3: 2})
            sage: v+w
            (0, 5, 0, 0)
            sage: (v+w).hamming_weight()
            1
        """
        return len(self._entries)

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        """
        Return a numerical approximation of ``self`` with ``prec`` bits
        (or decimal ``digits``) of precision, by approximating all
        entries.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- which algorithm to use to compute the
          approximation of the entries (the accepted algorithms depend
          on the object)

        If neither ``prec`` nor ``digits`` is given, the default
        precision is 53 bits (roughly 16 digits).

        EXAMPLES::

            sage: v = vector(RealField(200), [1,2,3], sparse=True)
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Sparse vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Sparse vector space of dimension 3 over Real Field with 75 bits of precision
        """
        if prec is None:
            prec = digits_to_bits(digits)
        return vector({k: v.numerical_approx(prec, algorithm=algorithm)
                for k, v in self._entries.iteritems()}, sparse=True)
