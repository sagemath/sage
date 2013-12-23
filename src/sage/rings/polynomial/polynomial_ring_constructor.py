r"""
Constructors for polynomial rings

This module provides the function :func:`PolynomialRing`, which constructs
rings of univariate and multivariate polynomials, and implements caching to
prevent the same ring being created in memory multiple times (which is
wasteful and breaks the general assumption in Sage that parents are unique).

There is also a function :func:`BooleanPolynomialRing_constructor`, used for
constructing Boolean polynomial rings, which are not technically polynomial
rings but rather quotients of them (see module
:mod:`sage.rings.polynomial.pbori` for more details); and a deprecated
constructor :func:`MPolynomialRing` (now subsumed by the generic
:meth:`PolynomialRing`.
"""

#################################################################
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
######################################################################


from sage.structure.parent_gens import normalize_names
from sage.structure.element import is_Element
import sage.rings.ring as ring
import sage.rings.padics.padic_base_leaves as padic_base_leaves

from sage.rings.integer import Integer
from sage.rings.finite_rings.constructor import is_FiniteField
from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing

from sage.misc.cachefunc import weak_cached_function

from sage.categories.fields import Fields
_Fields = Fields()
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
_UFD = UniqueFactorizationDomains()
from sage.categories.integral_domains import IntegralDomains
_ID = IntegralDomains()
from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()

import sage.misc.weak_dict
_cache = sage.misc.weak_dict.WeakValueDictionary()

def PolynomialRing(base_ring, arg1=None, arg2=None,
                   sparse=False, order='degrevlex',
                   names=None, name=None,
                   var_array=None,
                   implementation=None):
    r"""
    Return the globally unique univariate or multivariate polynomial
    ring with given properties and variable name or names.

    There are five ways to call the polynomial ring constructor:

    1. ``PolynomialRing(base_ring, name,    sparse=False)``
    2. ``PolynomialRing(base_ring, names,   order='degrevlex')``
    3. ``PolynomialRing(base_ring, name, n, order='degrevlex')``
    4. ``PolynomialRing(base_ring, n, name, order='degrevlex')``
    5. ``PolynomialRing(base_ring, n, var_array=var_array, order='degrevlex')``

    The optional arguments sparse and order *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:

    - ``base_ring`` -- a commutative ring
    - ``name`` -- a string
    - ``names`` -- a list or tuple of names, or a comma separated string
    - ``var_array`` -- a list or tuple of names, or a comma separated string
    - ``n`` -- an integer
    - ``sparse`` -- bool (default: False), whether or not elements are sparse
    - ``order`` -- string or
      :class:`~sage.rings.polynomial.term_order.TermOrder` object, e.g.,

      - ``'degrevlex'`` (default) -- degree reverse lexicographic
      - ``'lex'``  -- lexicographic
      - ``'deglex'`` -- degree lexicographic
      - ``TermOrder('deglex',3) + TermOrder('deglex',3)`` -- block ordering

    - ``implementation`` -- string or None; selects an implementation in cases
      where Sage includes multiple choices (currently `\ZZ[x]` can be
      implemented with 'NTL' or 'FLINT'; default is 'FLINT')

    .. NOTE::

        The following rules were introduced in :trac:`9944`, in order
        to preserve the "unique parent assumption" in Sage (i.e., if two
        parents evaluate equal then they should actually be identical).

        - In the multivariate case, a dense representation is not supported.
          Hence, the argument ``sparse=False`` is silently ignored in that case.
        - If the given implementation does not exist for rings with the given
          number of generators and the given sparsity, then an error results.

    OUTPUT:

    ``PolynomialRing(base_ring, name, sparse=False)`` returns a univariate
    polynomial ring; also, PolynomialRing(base_ring, names, sparse=False)
    yields a univariate polynomial ring, if names is a list or tuple
    providing exactly one name. All other input formats return a
    multivariate polynomial ring.

    UNIQUENESS and IMMUTABILITY: In Sage there is exactly one
    single-variate polynomial ring over each base ring in each choice
    of variable, sparseness, and implementation.  There is also exactly
    one multivariate polynomial ring over each base ring for each
    choice of names of variables and term order.  The names of the
    generators can only be temporarily changed after the ring has been
    created.  Do this using the localvars context:

        EXAMPLES of VARIABLE NAME CONTEXT::

            sage: R.<x,y> = PolynomialRing(QQ,2); R
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: f = x^2 - 2*y^2

        You can't just globally change the names of those variables.
        This is because objects all over Sage could have pointers to
        that polynomial ring.
        ::

            sage: R._assign_names(['z','w'])
            Traceback (most recent call last):
            ...
            ValueError: variable names cannot be changed after object creation.

        However, you can very easily change the names within a ``with`` block::

            sage: with localvars(R, ['z','w']):
            ...     print f
            ...
            z^2 - 2*w^2

        After the ``with`` block the names revert to what they were before.
        ::

            sage: print f
            x^2 - 2*y^2

    SQUARE BRACKETS NOTATION: You can alternatively create a single or
    multivariate polynomial ring over a ring `R` by writing ``R['varname']`` or
    ``R['var1,var2,var3,...']``.  This square brackets notation doesn't allow
    for setting any of the optional arguments.

    EXAMPLES:

    1. ``PolynomialRing(base_ring, name, sparse=False)``

       ::

        sage: PolynomialRing(QQ, 'w')
        Univariate Polynomial Ring in w over Rational Field

       Use the diamond brackets notation to make the variable
       ready for use after you define the ring::

        sage: R.<w> = PolynomialRing(QQ)
        sage: (1 + w)^3
        w^3 + 3*w^2 + 3*w + 1

       You must specify a name::

        sage: PolynomialRing(QQ)
        Traceback (most recent call last):
        ...
        TypeError: You must specify the names of the variables.

        sage: R.<abc> = PolynomialRing(QQ, sparse=True); R
        Sparse Univariate Polynomial Ring in abc over Rational Field

        sage: R.<w> = PolynomialRing(PolynomialRing(GF(7),'k')); R
        Univariate Polynomial Ring in w over Univariate Polynomial Ring in k over Finite Field of size 7

       The square bracket notation::

        sage: R.<y> = QQ['y']; R
        Univariate Polynomial Ring in y over Rational Field
        sage: y^2 + y
        y^2 + y

       In fact, since the diamond brackets on the left determine the
       variable name, you can omit the variable from the square brackets::

        sage: R.<zz> = QQ[]; R
        Univariate Polynomial Ring in zz over Rational Field
        sage: (zz + 1)^2
        zz^2 + 2*zz + 1

       This is exactly the same ring as what PolynomialRing returns::

        sage: R is PolynomialRing(QQ,'zz')
        True

       However, rings with different variables are different::

        sage: QQ['x'] == QQ['y']
        False

       Sage has two implementations of univariate polynomials over the
       integers, one based on NTL and one based on FLINT.  The default
       is FLINT. Note that FLINT uses a "more dense" representation for
       its polynomials than NTL, so in particular, creating a polynomial
       like 2^1000000 * x^1000000 in FLINT may be unwise.
       ::

        sage: ZxNTL = PolynomialRing(ZZ, 'x', implementation='NTL'); ZxNTL
        Univariate Polynomial Ring in x over Integer Ring (using NTL)
        sage: ZxFLINT = PolynomialRing(ZZ, 'x', implementation='FLINT'); ZxFLINT
        Univariate Polynomial Ring in x over Integer Ring
        sage: ZxFLINT is ZZ['x']
        True
        sage: ZxFLINT is PolynomialRing(ZZ, 'x')
        True
        sage: xNTL = ZxNTL.gen()
        sage: xFLINT = ZxFLINT.gen()
        sage: xNTL.parent()
        Univariate Polynomial Ring in x over Integer Ring (using NTL)
        sage: xFLINT.parent()
        Univariate Polynomial Ring in x over Integer Ring

       There is a coercion from the non-default to the default
       implementation, so the values can be mixed in a single
       expression::

        sage: (xNTL + xFLINT^2)
        x^2 + x

       The result of such an expression will use the default, i.e.,
       the FLINT implementation::

        sage: (xNTL + xFLINT^2).parent()
        Univariate Polynomial Ring in x over Integer Ring

    2. ``PolynomialRing(base_ring, names,   order='degrevlex')``

       ::

        sage: R = PolynomialRing(QQ, 'a,b,c'); R
        Multivariate Polynomial Ring in a, b, c over Rational Field

        sage: S = PolynomialRing(QQ, ['a','b','c']); S
        Multivariate Polynomial Ring in a, b, c over Rational Field

        sage: T = PolynomialRing(QQ, ('a','b','c')); T
        Multivariate Polynomial Ring in a, b, c over Rational Field

       All three rings are identical. ::

        sage: (R is S) and  (S is T)
        True

       There is a unique polynomial ring with each term order::

        sage: R = PolynomialRing(QQ, 'x,y,z', order='degrevlex'); R
        Multivariate Polynomial Ring in x, y, z over Rational Field
        sage: S = PolynomialRing(QQ, 'x,y,z', order='invlex'); S
        Multivariate Polynomial Ring in x, y, z over Rational Field
        sage: S is PolynomialRing(QQ, 'x,y,z', order='invlex')
        True
        sage: R == S
        False

       Note that a univariate polynomial ring is returned, if the list
       of names is of length one. If it is of length zero, a multivariate
       polynomial ring with no variables is returned.

       ::

        sage: PolynomialRing(QQ,["x"])
        Univariate Polynomial Ring in x over Rational Field
        sage: PolynomialRing(QQ,[])
        Multivariate Polynomial Ring in no variables over Rational Field

    3. ``PolynomialRing(base_ring, name, n, order='degrevlex')``

       If you specify a single name as a string and a number of
       variables, then variables labeled with numbers are created.

       ::

        sage: PolynomialRing(QQ, 'x', 10)
        Multivariate Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 over Rational Field

        sage: PolynomialRing(GF(7), 'y', 5)
        Multivariate Polynomial Ring in y0, y1, y2, y3, y4 over Finite Field of size 7

        sage: PolynomialRing(QQ, 'y', 3, sparse=True)
        Multivariate Polynomial Ring in y0, y1, y2 over Rational Field

       Note that a multivariate polynomial ring is returned when an
       explicit number is given.

       ::

        sage: PolynomialRing(QQ,"x",1)
        Multivariate Polynomial Ring in x over Rational Field
        sage: PolynomialRing(QQ,"x",0)
        Multivariate Polynomial Ring in no variables over Rational Field

       It is easy in Python to create fairly arbitrary variable names.  For
       example, here is a ring with generators labeled by the first 100
       primes::

        sage: R = PolynomialRing(ZZ, ['x%s'%p for p in primes(100)]); R
        Multivariate Polynomial Ring in x2, x3, x5, x7, x11, x13, x17, x19, x23, x29, x31, x37, x41, x43, x47, x53, x59, x61, x67, x71, x73, x79, x83, x89, x97 over Integer Ring

       By calling the
       :meth:`~sage.structure.category_object.CategoryObject.inject_variables`
       method, all those variable names are available for interactive use::

        sage: R.inject_variables()
        Defining x2, x3, x5, x7, x11, x13, x17, x19, x23, x29, x31, x37, x41, x43, x47, x53, x59, x61, x67, x71, x73, x79, x83, x89, x97
        sage: (x2 + x41 + x71)^2
        x2^2 + 2*x2*x41 + x41^2 + 2*x2*x71 + 2*x41*x71 + x71^2

    5. ``PolynomialRing(base_ring, n, m, var_array=var_array, order='degrevlex')``

       This creates an array of variables where each variables begins with an
       entry in ``var_array`` and is indexed from 0 to ``n-1``.

        sage: PolynomialRing(ZZ, 3, var_array=['x','y'])
        Multivariate Polynomial Ring in x0, y0, x1, y1, x2, y2 over Integer Ring
        sage: PolynomialRing(ZZ, 3, var_array='a,b')
        Multivariate Polynomial Ring in a0, b0, a1, b1, a2, b2 over Integer Ring

       If ``var_array`` is a single string, this creates an `m \times n`
       array of variables::

        sage: PolynomialRing(ZZ, 2, 3, var_array='m')
        Multivariate Polynomial Ring in m00, m01, m02, m10, m11, m12 over Integer Ring

       If ``var_array`` is a single string and `m` is not specified, this
       creates an `n \times n` array of variables::

        sage: PolynomialRing(ZZ, 2, var_array='m')
        Multivariate Polynomial Ring in m00, m01, m10, m11 over Integer Ring

    TESTS:

    We test here some changes introduced in :trac:`9944`.

    If there is no dense implementation for the given number of
    variables, then requesting a dense ring results yields the
    corresponding sparse ring::

        sage: R.<x,y> = QQ[]
        sage: S.<x,y> = PolynomialRing(QQ, sparse=False)
        sage: R is S
        True

    If the requested implementation is not known or not supported for
    the given number of variables and the given sparsity, then an
    error results::

        sage: R.<x> = PolynomialRing(ZZ, implementation='Foo')
        Traceback (most recent call last):
        ...
        ValueError: Unknown implementation Foo for ZZ[x]
        sage: R.<x,y> = PolynomialRing(ZZ, implementation='FLINT')
        Traceback (most recent call last):
        ...
        ValueError: The FLINT implementation is not known for multivariate polynomial rings

    The following corner case used to result in a warning message from
    ``libSingular``, and the generators of the resulting polynomial
    ring were not zero::

        sage: R = Integers(1)['x','y']
        sage: R.0 == 0
        True

    We verify that :trac:`13187` is fixed::

        sage: var('t')
        t
        sage: PolynomialRing(ZZ, name=t) == PolynomialRing(ZZ, name='t')
        True

    We verify that polynomials with interval coefficients from
    :trac:`7712` and :trac:`13760` are fixed::

        sage: P.<y,z> = PolynomialRing(RealIntervalField(2))
        sage: Q.<x> = PolynomialRing(P)
        sage: C = (y-x)^3
        sage: C(y/2)
        1.?*y^3
        sage: R.<x,y> = PolynomialRing(RIF,2)
        sage: RIF(-2,1)*x
        0.?e1*x

    """
    import sage.rings.polynomial.polynomial_ring as m

    if not var_array is None:
        # Make sure arg1 always corresponds to n and arg2 to m
        if arg2 is not None:
            arg1, arg2 = arg2, arg1
        if isinstance(var_array, str):
            if not ',' in var_array:
                if arg2 is None:
                    arg2 = arg1
                arg1 = ['%s%s%s'%(var_array, i, j) for i in range(arg2) for j in range(arg1)]
            else:
                var_array = var_array.split(',')
                if arg2 is not None and len(var_array) != arg2:
                    raise IndexError('the number of names must equal the number of base generators')
                arg1 = ['%s%s'%(x, i) for i in range(arg1) for x in var_array]
        elif isinstance(var_array, (list,tuple)):
            if arg2 is not None and len(var_array) != arg2:
                raise IndexError('the number of names must equal the number of base generators')
            arg1 = ['%s%s'%(x, i) for i in range(arg1) for x in var_array]
        else:
            raise TypeError("invalid input to PolynomialRing function; please see the docstring for that function")
        arg2 = len(arg1)
    else:
        if isinstance(arg1, (int, long, Integer)):
            arg1, arg2 = arg2, arg1
        if not names is None:
            arg1 = names
        elif not name is None:
            arg1 = name

    if is_Element(arg1) and not isinstance(arg1, (int, long, Integer)):
        arg1 = repr(arg1)
    if is_Element(arg2) and not isinstance(arg2, (int, long, Integer)):
        arg2 = repr(arg2)

    if not m.ring.is_Ring(base_ring):
        raise TypeError('base_ring must be a ring')

    if arg1 is None:
        raise TypeError("You must specify the names of the variables.")

    R = None
    if isinstance(arg1, (list, tuple)):
        arg1 = [str(x) for x in arg1]
    if isinstance(arg2, (list, tuple)):
        arg2 = [str(x) for x in arg2]
    if isinstance(arg2, (int, long, Integer)):
        # 3. PolynomialRing(base_ring, names, n, order='degrevlex'):
        if not isinstance(arg1, (list, tuple, str)):
            raise TypeError("You *must* specify the names of the variables.")
        n = int(arg2)
        names = arg1
        R = _multi_variate(base_ring, names, n, sparse, order, implementation)

    elif isinstance(arg1, str) or (isinstance(arg1, (list,tuple)) and len(arg1) == 1):
        if not ',' in arg1:
            # 1. PolynomialRing(base_ring, name, sparse=False):
            if not arg2 is None:
                raise TypeError("if second arguments is a string with no commas, then there must be no other non-optional arguments")
            name = arg1
            R = _single_variate(base_ring, name, sparse, implementation)
        else:
            # 2-4. PolynomialRing(base_ring, names, order='degrevlex'):
            if not arg2 is None:
                raise TypeError("invalid input to PolynomialRing function; please see the docstring for that function")
            names = arg1.split(',')
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order, implementation)
    elif isinstance(arg1, (list, tuple)):
            # PolynomialRing(base_ring, names (list or tuple), order='degrevlex'):
            names = arg1
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order, implementation)

    if arg1 is None and arg2 is None:
        raise TypeError("you *must* specify the indeterminates (as not None).")
    if R is None:
        raise TypeError("invalid input (%s, %s, %s) to PolynomialRing function; please see the docstring for that function"%(
            base_ring, arg1, arg2))

    return R

def MPolynomialRing(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future version of
    Sage. Please use PolynomialRing instead.

    If you have questions regarding this function and its replacement,
    please send your comments to sage-support@googlegroups.com.
    """
    from sage.misc.superseded import deprecation
    deprecation(6500, "MPolynomialRing is deprecated, use PolynomialRing instead!")
    return PolynomialRing(*args, **kwds)

def _get_from_cache(key):
    try:
        return _cache[key] #()
    except TypeError as msg:
        raise TypeError('key = %s\n%s'%(key,msg))
    except KeyError:
        return None

def _save_in_cache(key, R):
    try:
         _cache[key] = R
    except TypeError as msg:
        raise TypeError('key = %s\n%s'%(key,msg))


def _single_variate(base_ring, name, sparse, implementation):
    import sage.rings.polynomial.polynomial_ring as m
    name = normalize_names(1, name)
    key = (base_ring, name, sparse, implementation if not sparse else None)
    R = _get_from_cache(key)
    if not R is None:
        return R

    if isinstance(base_ring, ring.CommutativeRing):
        if is_IntegerModRing(base_ring) and not sparse:
            n = base_ring.order()
            if n.is_prime():
                R = m.PolynomialRing_dense_mod_p(base_ring, name, implementation=implementation)
            elif n > 1:
                R = m.PolynomialRing_dense_mod_n(base_ring, name, implementation=implementation)
            else:  # n == 1!
                R = m.PolynomialRing_integral_domain(base_ring, name)   # specialized code breaks in this case.

        elif is_FiniteField(base_ring) and not sparse:
            R = m.PolynomialRing_dense_finite_field(base_ring, name, implementation=implementation)

        elif isinstance(base_ring, padic_base_leaves.pAdicFieldCappedRelative):
            R = m.PolynomialRing_dense_padic_field_capped_relative(base_ring, name)

        elif isinstance(base_ring, padic_base_leaves.pAdicRingCappedRelative):
            R = m.PolynomialRing_dense_padic_ring_capped_relative(base_ring, name)

        elif isinstance(base_ring, padic_base_leaves.pAdicRingCappedAbsolute):
            R = m.PolynomialRing_dense_padic_ring_capped_absolute(base_ring, name)

        elif isinstance(base_ring, padic_base_leaves.pAdicRingFixedMod):
            R = m.PolynomialRing_dense_padic_ring_fixed_mod(base_ring, name)

        elif base_ring.is_field(proof = False):
            R = m.PolynomialRing_field(base_ring, name, sparse)

        elif base_ring.is_integral_domain(proof = False):
            R = m.PolynomialRing_integral_domain(base_ring, name, sparse, implementation)
        else:
            R = m.PolynomialRing_commutative(base_ring, name, sparse)
    else:
        R = m.PolynomialRing_general(base_ring, name, sparse)

    if hasattr(R, '_implementation_names'):
        for name in R._implementation_names:
            real_key = key[0:3] + (name,)
            _save_in_cache(real_key, R)
    else:
        _save_in_cache(key, R)
    return R

def _multi_variate(base_ring, names, n, sparse, order, implementation):
#    if not sparse:
#        raise ValueError, "A dense representation of multivariate polynomials is not supported"
    sparse = False
    # "True" would be correct, since there is no dense implementation of
    # multivariate polynomials. However, traditionally, "False" is used in the key,
    # even though it is meaningless.

    if implementation is not None:
        raise ValueError("The %s implementation is not known for multivariate polynomial rings"%implementation)

    names = normalize_names(n, names)

    import sage.rings.polynomial.multi_polynomial_ring as m
    from sage.rings.polynomial.term_order import TermOrder

    order = TermOrder(order, n)

    key = (base_ring, names, n, sparse, order)
    R = _get_from_cache(key)
    if not R is None:
        return R

    from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
    if m.integral_domain.is_IntegralDomain(base_ring):
        if n < 1:
            R = m.MPolynomialRing_polydict_domain(base_ring, n, names, order)
        else:
            try:
                R = MPolynomialRing_libsingular(base_ring, n, names, order)
            except ( TypeError, NotImplementedError ):
                R = m.MPolynomialRing_polydict_domain(base_ring, n, names, order)
    else:
        if not base_ring.is_zero():
            try:
                R = MPolynomialRing_libsingular(base_ring, n, names, order)
            except ( TypeError, NotImplementedError ):
                R = m.MPolynomialRing_polydict(base_ring, n, names, order)
        else:
            R = m.MPolynomialRing_polydict(base_ring, n, names, order)
    _save_in_cache(key, R)
    return R

#########################################################
# Choice of a category
from sage import categories
from sage.categories.algebras import Algebras
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.category import JoinCategory
# Some fixed categories, in order to avoid the function call overhead
_EuclideanDomains = categories.euclidean_domains.EuclideanDomains()
_UniqueFactorizationDomains = categories.unique_factorization_domains.UniqueFactorizationDomains()
_IntegralDomains = categories.integral_domains.IntegralDomains()
_Rings = category = categories.rings.Rings()

@weak_cached_function
def polynomial_default_category(base_ring,multivariate):
    """
    Choose an appropriate category for a polynomial ring.

    INPUT:

    - ``base_ring``: The ring over which the polynomial ring shall be defined.
    - ``multivariate``: Will the polynomial ring be multivariate?

    EXAMPLES::

        sage: QQ['t'].category() is Category.join([EuclideanDomains(), CommutativeAlgebras(QQ)])
        True
        sage: QQ['s','t'].category() is Category.join([UniqueFactorizationDomains(), CommutativeAlgebras(QQ)])
        True
        sage: QQ['s']['t'].category() is Category.join([UniqueFactorizationDomains(), CommutativeAlgebras(QQ['s'])])
        True

    """
    if base_ring in _Fields:
        if multivariate:
            return JoinCategory((_UniqueFactorizationDomains,CommutativeAlgebras(base_ring)))
        return JoinCategory((_EuclideanDomains,CommutativeAlgebras(base_ring)))
    if base_ring in _UFD: #base_ring.is_unique_factorization_domain():
        return JoinCategory((_UniqueFactorizationDomains,CommutativeAlgebras(base_ring)))
    if base_ring in _ID: #base_ring.is_integral_domain():
        return JoinCategory((_IntegralDomains,CommutativeAlgebras(base_ring)))
    if base_ring in _CommutativeRings: #base_ring.is_commutative():
        return CommutativeAlgebras(base_ring)
    return Algebras(base_ring)

def BooleanPolynomialRing_constructor(n=None, names=None, order="lex"):
    """
    Construct a boolean polynomial ring with the following
    parameters:

    INPUT:

    - ``n`` -- number of variables (an integer > 1)
    - ``names`` -- names of ring variables, may be a string or list/tuple of strings
    - ``order`` -- term order (default: lex)

    EXAMPLES::

        sage: R.<x, y, z> = BooleanPolynomialRing() # indirect doctest
        sage: R
        Boolean PolynomialRing in x, y, z

        sage: p = x*y + x*z + y*z
        sage: x*p
        x*y*z + x*y + x*z

        sage: R.term_order()
        Lexicographic term order

        sage: R = BooleanPolynomialRing(5,'x',order='deglex(3),deglex(2)')
        sage: R.term_order()
        Block term order with blocks:
        (Degree lexicographic term order of length 3,
         Degree lexicographic term order of length 2)

        sage: R = BooleanPolynomialRing(3,'x',order='degneglex')
        sage: R.term_order()
        Degree negative lexicographic term order

        sage: BooleanPolynomialRing(names=('x','y'))
        Boolean PolynomialRing in x, y

        sage: BooleanPolynomialRing(names='x,y')
        Boolean PolynomialRing in x, y

    TESTS::

        sage: P.<x,y> = BooleanPolynomialRing(2,order='deglex')
        sage: x > y
        True

        sage: P.<x0, x1, x2, x3> = BooleanPolynomialRing(4,order='deglex(2),deglex(2)')
        sage: x0 > x1
        True
        sage: x2 > x3
        True
    """

    if isinstance(n, str):
        names = n
        n = 0
    if n is None and names is not None:
        n = 0

    names = normalize_names(n, names)
    if n is 0:
        n = len(names)

    from sage.rings.polynomial.term_order import TermOrder

    order = TermOrder(order, n)

    key = ("pbori", names, n, order)
    R = _get_from_cache(key)
    if not R is None:
        return R

    from sage.rings.polynomial.pbori import BooleanPolynomialRing
    R = BooleanPolynomialRing(n, names, order)

    _save_in_cache(key, R)
    return R

#########################################################################################
# END (Factory function for making polynomial rings)
#########################################################################################
