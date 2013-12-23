# -*- coding: utf-8 -*-
r"""
Number Fields

AUTHORS:

- William Stein (2004, 2005): initial version

- Steven Sivek (2006-05-12): added support for relative extensions

- William Stein (2007-09-04): major rewrite and documentation

- Robert Bradshaw (2008-10): specified embeddings into ambient fields

- Simon King (2010-05): Improve coercion from GAP

- Jeroen Demeyer (2010-07, 2011-04): Upgrade PARI (#9343, #10430, #11130)

- Robert Harron (2012-08): added is_CM(), complex_conjugation(), and
  maximal_totally_real_subfield()

- Christian Stump (2012-11): added conversion to universal cyclotomic field

.. note::

   Unlike in PARI/GP, class group computations *in Sage* do *not* by default
   assume the Generalized Riemann Hypothesis. To do class groups computations
   not provably correctly you must often pass the flag ``proof=False`` to
   functions or call the function ``proof.number_field(False)``. It can easily
   take 1000's of times longer to do computations with ``proof=True`` (the
   default).

This example follows one in the Magma reference manual::

    sage: K.<y> = NumberField(x^4 - 420*x^2 + 40000)
    sage: z = y^5/11; z
    420/11*y^3 - 40000/11*y
    sage: R.<y> = PolynomialRing(K)
    sage: f = y^2 + y + 1
    sage: L.<a> = K.extension(f); L
    Number Field in a with defining polynomial y^2 + y + 1 over its base field
    sage: KL.<b> = NumberField([x^4 - 420*x^2 + 40000, x^2 + x + 1]); KL
    Number Field in b0 with defining polynomial x^4 - 420*x^2 + 40000 over its base field

We do some arithmetic in a tower of relative number fields::

    sage: K.<cuberoot2> = NumberField(x^3 - 2)
    sage: L.<cuberoot3> = K.extension(x^3 - 3)
    sage: S.<sqrt2> = L.extension(x^2 - 2)
    sage: S
    Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field
    sage: sqrt2 * cuberoot3
    cuberoot3*sqrt2
    sage: (sqrt2 + cuberoot3)^5
    (20*cuberoot3^2 + 15*cuberoot3 + 4)*sqrt2 + 3*cuberoot3^2 + 20*cuberoot3 + 60
    sage: cuberoot2 + cuberoot3
    cuberoot3 + cuberoot2
    sage: cuberoot2 + cuberoot3 + sqrt2
    sqrt2 + cuberoot3 + cuberoot2
    sage: (cuberoot2 + cuberoot3 + sqrt2)^2
    (2*cuberoot3 + 2*cuberoot2)*sqrt2 + cuberoot3^2 + 2*cuberoot2*cuberoot3 + cuberoot2^2 + 2
    sage: cuberoot2 + sqrt2
    sqrt2 + cuberoot2
    sage: a = S(cuberoot2); a
    cuberoot2
    sage: a.parent()
    Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field

.. warning::

   Doing arithmetic in towers of relative fields that depends on
   canonical coercions is currently VERY SLOW. It is much better to
   explicitly coerce all elements into a common field, then do
   arithmetic with them there (which is quite fast).
"""
#*****************************************************************************
#       Copyright (C) 2004, 2005, 2006, 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent_gens import localvars
from sage.misc.cachefunc import cached_method

import sage.libs.ntl.all as ntl
import sage.interfaces.gap
import sage.rings.arith

import sage.rings.complex_field
import sage.rings.real_mpfr
import sage.rings.real_mpfi
import sage.rings.complex_double
import sage.rings.real_double
import sage.rings.real_lazy

from sage.rings.finite_rings.integer_mod import mod
from sage.misc.functional import is_odd
from sage.misc.misc_c import prod
from sage.categories.homset import End

import sage.rings.ring
from sage.misc.latex import latex_variable_name
from sage.misc.misc import union

from unit_group import UnitGroup
from class_group import ClassGroup
from class_group import SClassGroup

from sage.structure.element import is_Element
from sage.structure.sequence import Sequence

import sage.structure.parent_gens

from sage.structure.proof.proof import get_flag
import maps
import number_field_morphisms
from itertools import count, izip


def is_NumberFieldHomsetCodomain(codomain):
    """
    Returns whether ``codomain`` is a valid codomain for a number
    field homset. This is used by NumberField._Hom_ to determine
    whether the created homsets should be a
    :class:`sage.rings.number_field.morphism.NumberFieldHomset`.

    EXAMPLES:

    This currently accepts any parent (CC, RR, ...) in :class:`Fields`::

        sage: from sage.rings.number_field.number_field import is_NumberFieldHomsetCodomain
        sage: is_NumberFieldHomsetCodomain(QQ)
        True
        sage: is_NumberFieldHomsetCodomain(NumberField(x^2 + 1, 'x'))
        True
        sage: is_NumberFieldHomsetCodomain(ZZ)
        False
        sage: is_NumberFieldHomsetCodomain(3)
        False
        sage: is_NumberFieldHomsetCodomain(MatrixSpace(QQ, 2))
        False
        sage: is_NumberFieldHomsetCodomain(InfinityRing)
        False

    Question: should, for example, QQ-algebras be accepted as well?

    Caveat: Gap objects are not (yet) in :class:`Fields`, and therefore
    not accepted as number field homset codomains::

        sage: is_NumberFieldHomsetCodomain(gap.Rationals)
        False
    """
    from sage.categories.fields import Fields
    return codomain in Fields()

from sage.rings.number_field.morphism import RelativeNumberFieldHomomorphism_from_abs

def proof_flag(t):
    """
    Used for easily determining the correct proof flag to use.

    Returns t if t is not None, otherwise returns the system-wide
    proof-flag for number fields (default: True).

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import proof_flag
        sage: proof_flag(True)
        True
        sage: proof_flag(False)
        False
        sage: proof_flag(None)
        True
        sage: proof_flag("banana")
        'banana'
    """
    return get_flag(t, "number_field")


import weakref

from sage.misc.latex import latex

import sage.rings.arith as arith
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.infinity as infinity
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.polynomial.polynomial_element as polynomial_element
import sage.rings.complex_field
import sage.groups.abelian_gps.abelian_group
import sage.rings.complex_interval_field

from sage.structure.parent_gens import ParentWithGens
import number_field_element
import number_field_element_quadratic
from number_field_ideal import is_NumberFieldIdeal, NumberFieldFractionalIdeal
from sage.libs.pari.all import pari, pari_gen

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()
RIF = sage.rings.real_mpfi.RealIntervalField()
CIF = sage.rings.complex_interval_field.ComplexIntervalField()
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.real_lazy import RLF, CLF


_nf_cache = {}
def NumberField(polynomial, name=None, check=True, names=None, cache=True,
                embedding=None, latex_name=None,
                assume_disc_small=False,
                maximize_at_primes=None):
    r"""
    Return *the* number field defined by the given irreducible
    polynomial and with variable with the given name. If check is True
    (the default), also verify that the defining polynomial is
    irreducible and over `\QQ`.

    INPUT:

        - ``polynomial`` - a polynomial over `\QQ` or a number field,
          or a list of polynomials.
        - ``name`` - a string (default: 'a'), the name of the generator
        - ``check`` - bool (default: True); do type checking and
          irreducibility checking.
        - ``embedding`` - image of the generator in an ambient field
          (default: None)
        - ``assume_disc_small`` -- (default: False); if True, assume
          that no square of a prime greater than PARI's primelimit
          (which should be 500000); only applies for absolute fields
          at present.
        - ``maximize_at_primes`` -- (default: None) None or a list of
          primes; if not None, then the maximal order is computed by
          maximizing only at the primes in this list, which completely
          avoids having to factor the discriminant, but of course can
          lead to wrong results; only applies for absolute fields at
          present.

    EXAMPLES::

        sage: z = QQ['z'].0
        sage: K = NumberField(z^2 - 2,'s'); K
        Number Field in s with defining polynomial z^2 - 2
        sage: s = K.0; s
        s
        sage: s*s
        2
        sage: s^2
        2

    Constructing a relative number field::

        sage: K.<a> = NumberField(x^2 - 2)
        sage: R.<t> = K[]
        sage: L.<b> = K.extension(t^3+t+a); L
        Number Field in b with defining polynomial t^3 + t + a over its base field
        sage: L.absolute_field('c')
        Number Field in c with defining polynomial x^6 + 2*x^4 + x^2 - 2
        sage: a*b
        a*b
        sage: L(a)
        a
        sage: L.lift_to_base(b^3 + b)
        -a

    Constructing another number field::

        sage: k.<i> = NumberField(x^2 + 1)
        sage: R.<z> = k[]
        sage: m.<j> = NumberField(z^3 + i*z + 3)
        sage: m
        Number Field in j with defining polynomial z^3 + i*z + 3 over its base field

    Number fields are globally unique::

        sage: K.<a>= NumberField(x^3-5)
        sage: a^3
        5
        sage: L.<a>= NumberField(x^3-5)
        sage: K is L
        True

    Having different defining polynomials makes the fields different::

        sage: x = polygen(QQ, 'x'); y = polygen(QQ, 'y')
        sage: k.<a> = NumberField(x^2 + 3)
        sage: m.<a> = NumberField(y^2 + 3)
        sage: k
        Number Field in a with defining polynomial x^2 + 3
        sage: m
        Number Field in a with defining polynomial y^2 + 3

    In case of conflict of the generator name with the name given by the preparser, the name given by the preparser takes precedence::

        sage: K.<b> = NumberField(x^2 + 5, 'a'); K
        Number Field in b with defining polynomial x^2 + 5

    One can also define number fields with specified embeddings, may be
    used for arithmetic and deduce relations with other number fields
    which would not be valid for an abstract number field::

        sage: K.<a> = NumberField(x^3-2, embedding=1.2)
        sage: RR.coerce_map_from(K)
        Composite map:
          From: Number Field in a with defining polynomial x^3 - 2
          To:   Real Field with 53 bits of precision
          Defn:   Generic morphism:
                  From: Number Field in a with defining polynomial x^3 - 2
                  To:   Real Lazy Field
                  Defn: a -> 1.259921049894873?
                then
                  Conversion via _mpfr_ method map:
                  From: Real Lazy Field
                  To:   Real Field with 53 bits of precision
        sage: RR(a)
        1.25992104989487
        sage: 1.1 + a
        2.35992104989487
        sage: b = 1/(a+1); b
        1/3*a^2 - 1/3*a + 1/3
        sage: RR(b)
        0.442493334024442
        sage: L.<b> = NumberField(x^6-2, embedding=1.1)
        sage: L(a)
        b^2
        sage: a + b
        b^2 + b

    Note that the image only needs to be specified to enough precision
    to distinguish roots, and is exactly computed to any needed
    precision::

        sage: RealField(200)(a)
        1.2599210498948731647672106072782283505702514647015079800820

    One can embed into any other field::

        sage: K.<a> = NumberField(x^3-2, embedding=CC.gen()-0.6)
        sage: CC(a)
        -0.629960524947436 + 1.09112363597172*I
        sage: L = Qp(5)
        sage: f = polygen(L)^3 - 2
        sage: K.<a> = NumberField(x^3-2, embedding=f.roots()[0][0])
        sage: a + L(1)
        4 + 2*5^2 + 2*5^3 + 3*5^4 + 5^5 + 4*5^6 + 2*5^8 + 3*5^9 + 4*5^12 + 4*5^14 + 4*5^15 + 3*5^16 + 5^17 + 5^18 + 2*5^19 + O(5^20)
        sage: L.<b> = NumberField(x^6-x^2+1/10, embedding=1)
        sage: K.<a> = NumberField(x^3-x+1/10, embedding=b^2)
        sage: a+b
        b^2 + b
        sage: CC(a) == CC(b)^2
        True
        sage: K.coerce_embedding()
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Number Field in b with defining polynomial x^6 - x^2 + 1/10
          Defn: a -> b^2

    The ``QuadraticField`` and ``CyclotomicField`` constructors
    create an embedding by default unless otherwise specified.

    ::

        sage: K.<zeta> = CyclotomicField(15)
        sage: CC(zeta)
        0.913545457642601 + 0.406736643075800*I
        sage: L.<sqrtn3> = QuadraticField(-3)
        sage: K(sqrtn3)
        2*zeta^5 + 1
        sage: sqrtn3 + zeta
        2*zeta^5 + zeta + 1

    An example involving a variable name that defines a function in
    PARI::

        sage: theta = polygen(QQ, 'theta')
        sage: M.<z> = NumberField([theta^3 + 4, theta^2 + 3]); M
        Number Field in z0 with defining polynomial theta^3 + 4 over its base field

    TESTS::

        sage: x = QQ['x'].gen()
        sage: y = ZZ['y'].gen()
        sage: K = NumberField(x^3 + x + 3, 'a'); K
        Number Field in a with defining polynomial x^3 + x + 3
        sage: K.defining_polynomial().parent()
        Univariate Polynomial Ring in x over Rational Field

    ::

        sage: L = NumberField(y^3 + y + 3, 'a'); L
        Number Field in a with defining polynomial y^3 + y + 3
        sage: L.defining_polynomial().parent()
        Univariate Polynomial Ring in y over Rational Field

    ::

        sage: sage.rings.number_field.number_field._nf_cache = {}
        sage: K.<x> = CyclotomicField(5)[]
        sage: W.<a> = NumberField(x^2 + 1); W
        Number Field in a with defining polynomial x^2 + 1 over its base field
        sage: sage.rings.number_field.number_field._nf_cache = {}
        sage: W1 = NumberField(x^2+1,'a')
        sage: K.<x> = CyclotomicField(5)[]
        sage: W.<a> = NumberField(x^2 + 1); W
        Number Field in a with defining polynomial x^2 + 1 over its base field

    The following has been fixed in :trac:`8800`::

        sage: P.<x> = QQ[]
        sage: K.<a> = NumberField(x^3-5,embedding=0)
        sage: L.<b> = K.extension(x^2+a)
        sage: F, R = L.construction()
        sage: F(R) == L    # indirect doctest
        True

    """
    if name is None and names is None:
        raise TypeError, "You must specify the name of the generator."
    if not names is None:
        name = names

    if isinstance(polynomial, (list, tuple)):
        return NumberFieldTower(polynomial, name, embeddings=embedding)

    name = sage.structure.parent_gens.normalize_names(1, name)

    if not isinstance(polynomial, polynomial_element.Polynomial):
        try:
            polynomial = polynomial.polynomial(QQ)
        except (AttributeError, TypeError):
            raise TypeError, "polynomial (=%s) must be a polynomial."%repr(polynomial)

    # convert ZZ to QQ
    R = polynomial.base_ring()
    Q = polynomial.parent().base_extend(R.fraction_field())
    polynomial = Q(polynomial)

    if cache:
        key = (polynomial, polynomial.base_ring(), name, latex_name,
               embedding, embedding.parent() if embedding is not None else None,
               assume_disc_small, None if maximize_at_primes is None else tuple(maximize_at_primes))
        if _nf_cache.has_key(key):
            K = _nf_cache[key]()
            if not K is None: return K

    if isinstance(R, NumberField_generic):
        S = R.extension(polynomial, name, check=check)
        if cache:
            _nf_cache[key] = weakref.ref(S)
        return S

    if polynomial.degree() == 2:
        K = NumberField_quadratic(polynomial, name, latex_name, check, embedding,
             assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
    else:
        K = NumberField_absolute(polynomial, name, latex_name, check, embedding,
             assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)

    if cache:
        _nf_cache[key] = weakref.ref(K)
    return K


def NumberFieldTower(v, names, check=True, embeddings=None):
    """
    Return the tower of number fields defined by the polynomials or
    number fields in the list v.

    This is the field `K_0` generated by a root of ``v[0]`` over its base field
    `K_1`, which is in turn generated by a root of ``v[1]`` over its base field
    `K_2`, etc, until `K_d = \QQ` where `d` is the size of the list `v`.  If
    all is False, then each ``v[i]`` must be irreducible over the previous
    fields.  Otherwise a list of all possible fields defined by all those
    polynomials is output.

    If names defines a variable name a, say, then the generators of the
    intermediate number fields are a0, a1, a2, ...

    INPUT:

    -  ``v`` - a list of polynomials or number fields
    -  ``names`` - variable names
    -  ``check`` - bool (default: True) only relevant if all is False. Then
       check irreducibility of each input polynomial.

    OUTPUT: a single number field or a list of number fields

    EXAMPLES::

        sage: k.<a,b,c> = NumberField([x^2 + 1, x^2 + 3, x^2 + 5]); k # indirect doctest
        Number Field in a with defining polynomial x^2 + 1 over its base field
        sage: a^2
        -1
        sage: b^2
        -3
        sage: c^2
        -5
        sage: (a+b+c)^2
        (2*b + 2*c)*a + 2*c*b - 9

    The Galois group is a product of 3 groups of order 2::

        sage: k.galois_group(type="pari")
        Galois group PARI group [8, 1, 3, "E(8)=2[x]2[x]2"] of degree 8 of the Number Field in a with defining polynomial x^2 + 1 over its base field

    Repeatedly calling base_field allows us to descend the internally
    constructed tower of fields::

        sage: k.base_field()
        Number Field in b with defining polynomial x^2 + 3 over its base field
        sage: k.base_field().base_field()
        Number Field in c with defining polynomial x^2 + 5
        sage: k.base_field().base_field().base_field()
        Rational Field

    In the following example the second polynomial is reducible over
    the first, so we get an error::

        sage: v = NumberField([x^3 - 2, x^3 - 2], names='a')
        Traceback (most recent call last):
        ...
        ValueError: defining polynomial (x^3 - 2) must be irreducible

    We mix polynomial parent rings::

        sage: k.<y> = QQ[]
        sage: m = NumberField([y^3 - 3, x^2 + x + 1, y^3 + 2], 'beta')
        sage: m
        Number Field in beta0 with defining polynomial y^3 - 3 over its base field
        sage: m.base_field ()
        Number Field in beta1 with defining polynomial x^2 + x + 1 over its base field

    A tower of quadratic fields::

        sage: K.<a> = NumberField([x^2 + 3, x^2 + 2, x^2 + 1])
        sage: K
        Number Field in a0 with defining polynomial x^2 + 3 over its base field
        sage: K.base_field()
        Number Field in a1 with defining polynomial x^2 + 2 over its base field
        sage: K.base_field().base_field()
        Number Field in a2 with defining polynomial x^2 + 1

    A bigger tower of quadratic fields::

        sage: K.<a2,a3,a5,a7> = NumberField([x^2 + p for p in [2,3,5,7]]); K
        Number Field in a2 with defining polynomial x^2 + 2 over its base field
        sage: a2^2
        -2
        sage: a3^2
        -3
        sage: (a2+a3+a5+a7)^3
        ((6*a5 + 6*a7)*a3 + 6*a7*a5 - 47)*a2 + (6*a7*a5 - 45)*a3 - 41*a5 - 37*a7

    The function can also be called by name::

        sage: NumberFieldTower([x^2 + 1, x^2 + 2], ['a','b'])
        Number Field in a with defining polynomial x^2 + 1 over its base field
    """
    # there is an "all" option below -- we do not use it, since
    # I couldn't get it to work with PARI reliably, and it isn't
    # very useful in practice.
    try:
        names = sage.structure.parent_gens.normalize_names(len(v), names)
    except IndexError:
        names = sage.structure.parent_gens.normalize_names(1, names)
        if len(v) > 1:
            names = ['%s%s'%(names[0], i) for i in range(len(v))]

    if embeddings is None:
        embeddings = [None] * len(v)

    if not isinstance(v, (list, tuple)):
        raise TypeError, "v must be a list or tuple"
    if len(v) == 0:
        return QQ
    if len(v) == 1:
        #if all:
        #    f = v[0]
        #    if not isinstance(f, polynomial_element.Polynomial):
        #        f = QQ['x'](v[0])
        #    F = f.factor()
        #    return [NumberField(F[i][0], names='%s%s'%(names[0],str(i))) for i in range(len(F))]
        #else:
        return NumberField(v[0], names=names, embedding=embeddings[0])
    f = v[0]
    w = NumberFieldTower(v[1:], names=names[1:], embeddings=embeddings[1:])
    if isinstance(f, polynomial_element.Polynomial):
        var = f.variable_name()
    else:
        var = 'x'

    name = names[0]
    #if all:
    #    R = w[0][var]  # polynomial ring
    #else:
    R = w[var]  # polynomial ring

    f = R(f)
    i = 0

    sep = chr(ord(name[0]) + 1)
    #if all:
    #    z = []
    #    for k in w:
    #        for g, _ in f.factor():
    #            z.append(k.extension(g, names='%s%s%s'%(name, sep, str(i)), check=False))
    #            i += 1
    #    return z
    #else:
    return w.extension(f, name, check=check, embedding=embeddings[0])


def QuadraticField(D, name='a', check=True, embedding=True, latex_name='sqrt', **args):
    r"""
    Return a quadratic field obtained by adjoining a square root of
    `D` to the rational numbers, where `D` is not a
    perfect square.

    INPUT:


    -  ``D`` - a rational number

    -  ``name`` - variable name (default: 'a')

    -  ``check`` - bool (default: True)

    -  ``embedding`` - bool or square root of D in an
       ambient field (default: True)

    - ``latex_name`` - latex variable name (default: \sqrt{D})


    OUTPUT: A number field defined by a quadratic polynomial. Unless
    otherwise specified, it has an embedding into `\RR` or
    `\CC` by sending the generator to the positive
    or upper-half-plane root.

    EXAMPLES::

        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3
        sage: K.<theta> = QuadraticField(3); K
        Number Field in theta with defining polynomial x^2 - 3
        sage: RR(theta)
        1.73205080756888
        sage: QuadraticField(9, 'a')
        Traceback (most recent call last):
        ...
        ValueError: D must not be a perfect square.
        sage: QuadraticField(9, 'a', check=False)
        Number Field in a with defining polynomial x^2 - 9

    Quadratic number fields derive from general number fields.

    ::

        sage: from sage.rings.number_field.number_field import is_NumberField
        sage: type(K)
        <class 'sage.rings.number_field.number_field.NumberField_quadratic_with_category'>
        sage: is_NumberField(K)
        True

    Quadratic number fields are cached::

        sage: QuadraticField(-11, 'a') is QuadraticField(-11, 'a')
        True

    By default, quadratic fields come with a nice latex representation::

        sage: K.<a> = QuadraticField(-7)
        sage: latex(K)
        \Bold{Q}(\sqrt{-7})
        sage: latex(a)
        \sqrt{-7}
        sage: latex(1/(1+a))
        -\frac{1}{8} \sqrt{-7} + \frac{1}{8}
        sage: K.latex_variable_name()
        '\\sqrt{-7}'

    We can provide our own name as well::

        sage: K.<a> = QuadraticField(next_prime(10^10), latex_name=r'\sqrt{D}')
        sage: 1+a
        a + 1
        sage: latex(1+a)
        \sqrt{D} + 1
        sage: latex(QuadraticField(-1, 'a', latex_name=None).gen())
        a

    The name of the generator does not interfere with Sage preparser, see #1135::

        sage: K1 = QuadraticField(5, 'x')
        sage: K2.<x> = QuadraticField(5)
        sage: K3.<x> = QuadraticField(5, 'x')
        sage: K1 is K2
        True
        sage: K1 is K3
        True
        sage: K1
        Number Field in x with defining polynomial x^2 - 5


    Note that, in presence of two different names for the generator,
    the name given by the preparser takes precedence::

        sage: K4.<y> = QuadraticField(5, 'x'); K4
        Number Field in y with defining polynomial x^2 - 5
        sage: K1 == K4
        False

    TESTS::

        sage: QuadraticField(-11, 'a') is QuadraticField(-11, 'a', latex_name='Z')
        False
        sage: QuadraticField(-11, 'a') is QuadraticField(-11, 'a', latex_name=None)
        False
    """
    D = QQ(D)
    if check:
        if D.is_square():
            raise ValueError, "D must not be a perfect square."
    R = QQ['x']
    f = R([-D, 0, 1])
    if embedding is True:
        if D > 0:
            embedding = RLF(D).sqrt()
        else:
            embedding = CLF(D).sqrt()
    if latex_name == 'sqrt':
        latex_name = r'\sqrt{%s}' % D
    return NumberField(f, name, check=False, embedding=embedding, latex_name=latex_name, **args)

def is_AbsoluteNumberField(x):
    """
    Return True if x is an absolute number field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import is_AbsoluteNumberField
        sage: is_AbsoluteNumberField(NumberField(x^2+1,'a'))
        True
        sage: is_AbsoluteNumberField(NumberField([x^3 + 17, x^2+1],'a'))
        False

    The rationals are a number field, but they're not of the absolute
    number field class.

    ::

        sage: is_AbsoluteNumberField(QQ)
        False
    """
    return isinstance(x, NumberField_absolute)

def is_QuadraticField(x):
    r"""
    Return True if x is of the quadratic *number* field type.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import is_QuadraticField
        sage: is_QuadraticField(QuadraticField(5,'a'))
        True
        sage: is_QuadraticField(NumberField(x^2 - 5, 'b'))
        True
        sage: is_QuadraticField(NumberField(x^3 - 5, 'b'))
        False

    A quadratic field specially refers to a number field, not a finite
    field::

        sage: is_QuadraticField(GF(9,'a'))
        False
    """
    return isinstance(x, NumberField_quadratic)


_cyclo_cache = {}
def CyclotomicField(n=0, names=None, bracket="()", embedding=True):
    r"""
    Return the `n`-th cyclotomic field, where n is a positive integer,
    or the universal cyclotomic field if ``n==0``.

    For the documentation of the universal cyclotomic field, see :class:`~sage.rings.universal_cyclotomic_field.universal_cyclotomic_field.UniversalCyclotomicField`.

    INPUT:

    -  ``n`` - a nonnegative integer, default:``0``

    -  ``names`` - name of generator (optional - defaults to zetan)

    - ``bracket`` - Defines the brackets in the case of ``n==0``, and
      is ignored otherwise. Can be any even length string, with ``"()"`` being the default.

    -  ``embedding`` - bool or n-th root of unity in an
       ambient field (default True)

    EXAMPLES:

    If called without a parameter, we get the :class:`universal cyclotomic field<sage.rings.universal_cyclotomic_field.universal_cyclotomic_field.UniversalCyclotomicField>`::

        sage: CyclotomicField()
        Universal Cyclotomic Field

    We create the `7`\th cyclotomic field
    `\QQ(\zeta_7)` with the default generator name.

    ::

        sage: k = CyclotomicField(7); k
        Cyclotomic Field of order 7 and degree 6
        sage: k.gen()
        zeta7

    The default embedding sends the generator to the complex primitive
    `n^{th}` root of unity of least argument.

    ::

        sage: CC(k.gen())
        0.623489801858734 + 0.781831482468030*I

    Cyclotomic fields are of a special type.

    ::

        sage: type(k)
        <class 'sage.rings.number_field.number_field.NumberField_cyclotomic_with_category'>

    We can specify a different generator name as follows.

    ::

        sage: k.<z7> = CyclotomicField(7); k
        Cyclotomic Field of order 7 and degree 6
        sage: k.gen()
        z7

    The `n` must be an integer.

    ::

        sage: CyclotomicField(3/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    The degree must be nonnegative.

    ::

        sage: CyclotomicField(-1)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be a positive integer

    The special case `n=1` does *not* return the rational
    numbers::

        sage: CyclotomicField(1)
        Cyclotomic Field of order 1 and degree 1

    Due to their default embedding into `\CC`,
    cyclotomic number fields are all compatible.

    ::

        sage: cf30 = CyclotomicField(30)
        sage: cf5 = CyclotomicField(5)
        sage: cf3 = CyclotomicField(3)
        sage: cf30.gen() + cf5.gen() + cf3.gen()
        zeta30^6 + zeta30^5 + zeta30 - 1
        sage: cf6 = CyclotomicField(6) ; z6 = cf6.0
        sage: cf3 = CyclotomicField(3) ; z3 = cf3.0
        sage: cf3(z6)
        zeta3 + 1
        sage: cf6(z3)
        zeta6 - 1
        sage: cf9 = CyclotomicField(9) ; z9 = cf9.0
        sage: cf18 = CyclotomicField(18) ; z18 = cf18.0
        sage: cf18(z9)
        zeta18^2
        sage: cf9(z18)
        -zeta9^5
        sage: cf18(z3)
        zeta18^3 - 1
        sage: cf18(z6)
        zeta18^3
        sage: cf18(z6)**2
        zeta18^3 - 1
        sage: cf9(z3)
        zeta9^3
    """
    n = ZZ(n)
    if n == 0:
        from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import UniversalCyclotomicField
        return UniversalCyclotomicField(names=names, bracket=bracket, embedding=embedding)
    if n < 0:
        raise ValueError, "n (=%s) must be a positive integer"%n

    if names is None:
        names = "zeta%s"%n
    names = sage.structure.parent_gens.normalize_names(1, names)
    if embedding is True:
        embedding = (2 * CLF.pi() * CLF.gen() / n).exp()
    key = (n, names, embedding)
    if _cyclo_cache.has_key(key):
        K = _cyclo_cache[key]()
        if not K is None: return K
    K = NumberField_cyclotomic(n, names, embedding=embedding)
    _cyclo_cache[key] = weakref.ref(K)
    return K

def is_CyclotomicField(x):
    """
    Return True if x is a cyclotomic field, i.e., of the special
    cyclotomic field class. This function does not return True for a
    number field that just happens to be isomorphic to a cyclotomic
    field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import is_CyclotomicField
        sage: is_CyclotomicField(NumberField(x^2 + 1,'zeta4'))
        False
        sage: is_CyclotomicField(CyclotomicField(4))
        True
        sage: is_CyclotomicField(CyclotomicField(1))
        True
        sage: is_CyclotomicField(QQ)
        False
        sage: is_CyclotomicField(7)
        False
    """
    return isinstance(x, NumberField_cyclotomic)

import number_field_base

is_NumberField = number_field_base.is_NumberField

class NumberField_generic(number_field_base.NumberField):
    """
    EXAMPLES::

        sage: K.<a> = NumberField(x^3 - 2); K
        Number Field in a with defining polynomial x^3 - 2
        sage: TestSuite(K).run()
    """
    def __init__(self, polynomial, name,
                 latex_name=None, check=True, embedding=None,
                 category = None,
                 assume_disc_small=False, maximize_at_primes=None):
        """
        Create a number field.

        EXAMPLES::

            sage: NumberField(x^97 - 19, 'a')
            Number Field in a with defining polynomial x^97 - 19

        The defining polynomial must be irreducible::

            sage: K.<a> = NumberField(x^2 - 1)
            Traceback (most recent call last):
            ...
            ValueError: defining polynomial (x^2 - 1) must be irreducible

        If you use check=False, you avoid checking irreducibility of the
        defining polynomial, which can save time.

        ::

            sage: K.<a> = NumberField(x^2 - 1, check=False)

        It can also be dangerous::

            sage: (a-1)*(a+1)
            0

        The constructed object is in the category of number fields::

            sage: NumberField(x^2 + 3, 'a').category()
            Category of number fields
            sage: category(NumberField(x^2 + 3, 'a'))
            Category of number fields

        The special types of number fields, e.g., quadratic fields, do
        not have (yet?) their own category::

            sage: QuadraticField(2,'d').category()
            Category of number fields

        TESTS::

            sage: NumberField(ZZ['x'].0^4 + 23, 'a')
            Number Field in a with defining polynomial x^4 + 23
            sage: NumberField(QQ['x'].0^4 + 23, 'a')
            Number Field in a with defining polynomial x^4 + 23
            sage: NumberField(GF(7)['x'].0^4 + 23, 'a')
            Traceback (most recent call last):
            ...
            TypeError: polynomial must be defined over rational field
        """
        self._assume_disc_small = assume_disc_small
        self._maximize_at_primes = maximize_at_primes
        from sage.categories.number_fields import NumberFields
        default_category = NumberFields()
        if category is None:
            category = default_category
        else:
            assert category.is_subcategory(default_category), "%s is not a subcategory of %s"%(category, default_category)

        ParentWithGens.__init__(self, QQ, name, category=category)
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError, "polynomial (=%s) must be a polynomial"%repr(polynomial)

        if check:
            if not polynomial.parent().base_ring() == QQ:
                raise TypeError, "polynomial must be defined over rational field"
            if not polynomial.is_monic():
                raise NotImplementedError, "number fields for non-monic polynomials not yet implemented."
            if not polynomial.is_irreducible():
                raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial

        self._assign_names(name)
        if latex_name is None:
            self.__latex_variable_name = latex_variable_name(self.variable_name())
        else:
            self.__latex_variable_name = latex_name
        self.__polynomial = polynomial
        self.__pari_bnf_certified = False
        self._integral_basis_dict = {}
        embedding = number_field_morphisms.create_embedding_from_approx(self, embedding)
        self._populate_coercion_lists_(embedding=embedding)

    @cached_method
    def _magma_polynomial_(self, magma):
        """
        Return Magma version of the defining polynomial of this number field.

        EXAMPLES::

            sage: R.<x> = QQ[]                                                   # optional - magma
            sage: K.<a> = NumberField(x^3+2)                                     # optional - magma
            sage: K._magma_polynomial_(magma)                                    # optional - magma
            x^3 + 2
            sage: magma2=Magma()                                                 # optional - magma
            sage: K._magma_polynomial_(magma2)                                   # optional - magma
            x^3 + 2
            sage: K._magma_polynomial_(magma) is K._magma_polynomial_(magma)     # optional - magma
            True
            sage: K._magma_polynomial_(magma) is K._magma_polynomial_(magma2)    # optional - magma
            False
        """
        # NB f must not be garbage-collected, otherwise the
        # return value of this function is invalid
        return magma(self.defining_polynomial())

    def _magma_init_(self, magma):
        """
        Return a Magma version of this number field.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: K.<a> = NumberField(t^2 + 1)
            sage: K._magma_init_(magma)                            # optional - magma
            'SageCreateWithNames(NumberField(_sage_[...]),["a"])'
            sage: L = magma(K)    # optional - magma
            sage: L               # optional - magma
            Number Field with defining polynomial t^2 + 1 over the Rational Field
            sage: L.sage()        # optional - magma
            Number Field in a with defining polynomial t^2 + 1
            sage: L.sage() is K   # optional - magma
            True
            sage: L.1             # optional - magma
            a
            sage: L.1^2           # optional - magma
            -1
            sage: m = magma(a+1/2); m    # optional - magma
            1/2*(2*a + 1)
            sage: m.sage()        # optional - magma
            a + 1/2

        A relative number field::

            sage: S.<u> = K[]
            sage: M.<b> = NumberField(u^3+u+a)
            sage: L = magma(M)    # optional - magma
            sage: L               # optional - magma
            Number Field with defining polynomial u^3 + u + a over its ground field
            sage: L.sage() is M   # optional - magma
            True
        """
        # Get magma version of defining polynomial of this number field
        f = self._magma_polynomial_(magma)
        s = 'NumberField(%s)'%f.name()
        return magma._with_names(s, self.variable_names())

    def construction(self):
        r"""
        Construction of self

        EXAMPLE::

            sage: K.<a>=NumberField(x^3+x^2+1,embedding=CC.gen())
            sage: F,R = K.construction()
            sage: F
            AlgebraicExtensionFunctor
            sage: R
            Rational Field
            sage: F(R) == K
            True

        Note that, if a number field is provided with an embedding,
        the construction functor applied to the rationals is not
        necessarily identic with the number field. The reason is
        that the construction functor uses a value for the embedding
        that is equivalent, but not necessarily equal, to the one
        provided in the definition of the number field::

            sage: F(R) is K
            False
            sage: F.embeddings
            [0.2327856159383841? + 0.7925519925154479?*I]

        TEST::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: R.<t> = ZZ[]
            sage: a+t     # indirect doctest
            t + a
            sage: (a+t).parent()
            Univariate Polynomial Ring in t over Number Field in a with defining polynomial x^3 + x + 1

        The construction works for non-absolute number fields as well::

            sage: K.<a,b,c>=NumberField([x^3+x^2+1,x^2+1,x^7+x+1])
            sage: F,R = K.construction()
            sage: F(R) == K
            True

        ::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-5,embedding=0)
            sage: L.<b> = K.extension(x^2+a)
            sage: a*b
            a*b

        """
        from sage.categories.pushout import AlgebraicExtensionFunctor
        from sage.all import QQ
        if self.is_absolute():
            return (AlgebraicExtensionFunctor([self.polynomial()], [self.variable_name()], [None if self.coerce_embedding() is None else self.coerce_embedding()(self.gen())]), QQ)
        names = self.variable_names()
        polys = []
        embeddings = []
        K = self
        while (1):
            if K.is_absolute():
                break
            polys.append(K.relative_polynomial())
            embeddings.append(None if K.coerce_embedding() is None else K.coerce_embedding()(self.gen()))
            K = K.base_field()
        polys.append(K.relative_polynomial())
        embeddings.append(None if K.coerce_embedding() is None else K.coerce_embedding()(K.gen()))
        return (AlgebraicExtensionFunctor(polys, names, embeddings), QQ)

    def _element_constructor_(self, x):
        r"""
        Make x into an element of this number field, possibly not canonically.

        INPUT:

            - ``x`` - the element

        OUTPUT:

            ``x``, as an element of this number field

        TESTS::

            sage: K.<a> = NumberField(x^3 + 17)
            sage: K(a) is a # indirect doctest
            True
            sage: K('a^2 + 2/3*a + 5')
            a^2 + 2/3*a + 5
            sage: K('1').parent()
            Number Field in a with defining polynomial x^3 + 17
            sage: K(3/5).parent()
            Number Field in a with defining polynomial x^3 + 17
            sage: K.<a> = NumberField(polygen(QQ)^2 - 5)
            sage: F.<b> = K.extension(polygen(K))
            sage: F([a])
            a

        Using a GAP element may be tricky, as it may contain
        an exclamation mark::

            sage: L.<tau> = NumberField(x^3-2)
            sage: gap(tau^3)
            2
            sage: gap(tau)^3
            !2
            sage: L(gap(tau)^3)     # indirect doctest
            2

        """
        if isinstance(x, number_field_element.NumberFieldElement):
            K = x.parent()
            if K is self:
                return x
            elif K == self:
                return self._element_class(self, x.polynomial())
            elif isinstance(x, (number_field_element.OrderElement_absolute,
                                number_field_element.OrderElement_relative,
                                number_field_element_quadratic.OrderElement_quadratic)):
                L = K.number_field()
                if L == self:
                    return self._element_class(self, x)
                x = L(x)
            return self._coerce_from_other_number_field(x)
        elif sage.interfaces.gap.is_GapElement(x):
            s = repr(x)
            if self.variable_name() in s:
                return self._coerce_from_str(s)
            return self._coerce_from_str(s.replace('!',''))
        elif isinstance(x,str):
            return self._coerce_from_str(x)
        elif isinstance(x, (tuple, list)) or \
                (isinstance(x, sage.modules.free_module_element.FreeModuleElement) and
                 self.base_ring().has_coerce_map_from(x.parent().base_ring())):
            if len(x) != self.relative_degree():
                raise ValueError, "Length must be equal to the degree of this number field"
            result = x[0]
            for i in xrange(1,self.relative_degree()):
                result += x[i]*self.gen(0)**i
            return result
        return self._coerce_non_number_field_element_in(x)

    def _coerce_from_str(self, x):
        """
        Coerce a string representation of an element of this
        number field into this number field.

        INPUT:
            x -- string

        EXAMPLES:
            sage: k.<theta25> = NumberField(x^3+(2/3)*x+1)
            sage: k._coerce_from_str('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1

        This function is called by the coerce method when it gets a string
        as input:
            sage: k('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1
        """
        # provide string coercion, as
        # for finite fields
        w = sage.misc.all.sage_eval(x,locals=self.gens_dict())
        if not (is_Element(w) and w.parent() is self):
            return self(w)
        else:
            return w

    def _coerce_from_universal_cyclotomic_field(self, x):
        """
        Coerce an element of the universal cyclotomic field
        into this cyclotomic field.

        EXAMPLES::

            sage: CF = CyclotomicField(5)
            sage: UCF.<E> = UniversalCyclotomicField()

            sage: CF(E(5)) # indirect doctest
            zeta5

            sage: CF = CyclotomicField(7)
            sage: CF(E(5)) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: The element E(5) cannot be converted to Cyclotomic Field of order 7 and degree 6

            sage: CF = CyclotomicField(10)
            sage: CF(E(5)) # indirect doctest
            zeta10^2

        Matrices are correctly dealt with::

            sage: M = Matrix(UCF,2,[E(3),E(4),E(5),E(6)]); M
            [   E(3)    E(4)]
            [   E(5) -E(3)^2]

            sage: Matrix(CyclotomicField(60),M) # indirect doctest
            [zeta60^10 - 1     zeta60^15]
            [    zeta60^12     zeta60^10]

        Using a non-standard embedding::

            sage: CF = CyclotomicField(5,embedding=CC(exp(4*pi*i/5)))
            sage: x = E(5)
            sage: CC(x)
            0.309016994374947 + 0.951056516295154*I
            sage: CC(CF(x))
            0.309016994374947 + 0.951056516295154*I

            sage: UCF.<E> = UniversalCyclotomicField(embedding=None)
            sage: CF(E(5))
            Traceback (most recent call last):
            ...
            TypeError: The element E(5) cannot be converted to Cyclotomic Field of order 5 and degree 4
        """
        from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import UniversalCyclotomicField
        assert isinstance(x,UniversalCyclotomicField.Element)
        if not x.parent()._has_standard_embedding or self._n() % x.field_order():
            raise TypeError("The element %s cannot be converted to %s"%(x,self))
        return self(x.to_cyclotomic_field())

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from self to the number field
        codomain.

        The cat option is currently ignored.

        EXAMPLES: This function is implicitly called by the Hom method or
        function.

        ::

            sage: K.<i> = NumberField(x^2 + 1); K
            Number Field in i with defining polynomial x^2 + 1
            sage: K.Hom(K) # indirect doctest
            Automorphism group of Number Field in i with defining polynomial x^2 + 1
            sage: Hom(K, QuadraticField(-1, 'b'))
            Set of field embeddings from Number Field in i with defining polynomial x^2 + 1 to Number Field in b with defining polynomial x^2 + 1

        CHECKME: handling of the case where codomain is not a number field?

           sage: Hom(K, VectorSpace(QQ,3))
           Set of Morphisms from Number Field in i with defining polynomial x^2 + 1 to Vector space of dimension 3 over Rational Field in Category of commutative additive groups
        """
        if is_NumberFieldHomsetCodomain(codomain):
            import morphism
            return morphism.NumberFieldHomset(self, codomain)
        else:
            raise TypeError

    def _set_structure(self, from_self, to_self, unsafe_force_change=False):
        """
        Internal function to set the structure fields of this number
        field.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-23)
            sage: L.<b> = K.change_names()
            sage: L.structure() # indirect doctest
            (Isomorphism given by variable name change map:
              From: Number Field in b with defining polynomial x^2 + 23
              To:   Number Field in a with defining polynomial x^2 + 23, Isomorphism given by variable name change map:
              From: Number Field in a with defining polynomial x^2 + 23
              To:   Number Field in b with defining polynomial x^2 + 23)
        """
        # Note -- never call this on a cached number field, since
        # that could eventually lead to problems.
        if unsafe_force_change:
            self.__from_self = from_self
            self.__to_self = to_self
            return
        try:
            self.__from_self
        except AttributeError:
            self.__from_self = from_self
            self.__to_self = to_self
        else:
            raise ValueError, "number field structure is immutable."

    def structure(self):
        """
        Return fixed isomorphism or embedding structure on self.

        This is used to record various isomorphisms or embeddings that
        arise naturally in other constructions.

        EXAMPLES::

            sage: K.<z> = NumberField(x^2 + 3)
            sage: L.<a> = K.absolute_field(); L
            Number Field in a with defining polynomial x^2 + 3
            sage: L.structure()
            (Isomorphism given by variable name change map:
              From: Number Field in a with defining polynomial x^2 + 3
              To:   Number Field in z with defining polynomial x^2 + 3,
             Isomorphism given by variable name change map:
              From: Number Field in z with defining polynomial x^2 + 3
              To:   Number Field in a with defining polynomial x^2 + 3)

            sage: K.<a> = QuadraticField(-3)
            sage: R.<y> = K[]
            sage: D.<x0> = K.extension(y)
            sage: D_abs.<y0> = D.absolute_field()
            sage: D_abs.structure()[0](y0)
            -a
        """
        try:
            if self.__to_self is not None:
                return self.__from_self, self.__to_self
            else:
                return self.__from_self, self.__to_self
        except AttributeError:
            f = self.hom(self)
            self._set_structure(f,f)
            return f, f

    def completion(self, p, prec, extras={}):
        """
        Returns the completion of self at `p` to the specified
        precision. Only implemented at archimedean places, and then only if
        an embedding has been fixed.

        EXAMPLES::

            sage: K.<a> = QuadraticField(2)
            sage: K.completion(infinity, 100)
            Real Field with 100 bits of precision
            sage: K.<zeta> = CyclotomicField(12)
            sage: K.completion(infinity, 53, extras={'type': 'RDF'})
            Complex Double Field
            sage: zeta + 1.5                            # implicit test
            2.36602540378444 + 0.500000000000000*I
        """
        if p == infinity.infinity:
            gen_image = self.gen_embedding()
            if gen_image is not None:
                if gen_image in RDF:
                    return QQ.completion(p, prec, extras)
                elif gen_image in CDF:
                    return QQ.completion(p, prec, extras).algebraic_closure()
            raise ValueError, "No embedding into the complex numbers has been specified."
        else:
            raise NotImplementedError

    def primitive_element(self):
        r"""
        Return a primitive element for this field, i.e., an element that
        generates it over `\QQ`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.primitive_element()
            a
            sage: K.<a,b,c> = NumberField([x^2-2,x^2-3,x^2-5])
            sage: K.primitive_element()
            a - b + c
            sage: alpha = K.primitive_element(); alpha
            a - b + c
            sage: alpha.minpoly()
            x^2 + (2*b - 2*c)*x - 2*c*b + 6
            sage: alpha.absolute_minpoly()
            x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
        """
        try:
            return self.__primitive_element
        except AttributeError:
            pass
        K = self.absolute_field('a')
        from_K, to_K = K.structure()
        self.__primitive_element = from_K(K.gen())
        return self.__primitive_element

    def random_element(self, num_bound=None, den_bound=None,
                       integral_coefficients=False, distribution=None):
        r"""
        Return a random element of this number field.

        INPUT:

        - ``num_bound`` - Bound on numerator of the coefficients of
                          the resulting element

        - ``den_bound`` - Bound on denominators of the coefficients
                          of the resulting element

        - ``integral_coefficients`` (default: False) - If True, then
                          the resulting element will have integral
                          coefficients. This option overrides any
                          value of `den_bound`.

        - ``distribution`` - Distribution to use for the coefficients
                          of the resulting element

        OUTPUT:

        - Element of this number field

        EXAMPLES::

            sage: K.<j> = NumberField(x^8+1)
            sage: K.random_element()
            1/2*j^7 - j^6 - 12*j^5 + 1/2*j^4 - 1/95*j^3 - 1/2*j^2 - 4

            sage: K.<a,b,c> = NumberField([x^2-2,x^2-3,x^2-5])
            sage: K.random_element()
            ((6136*c - 7489/3)*b + 5825/3*c - 71422/3)*a + (-4849/3*c + 58918/3)*b - 45718/3*c + 75409/12

            sage: K.<a> = NumberField(x^5-2)
            sage: K.random_element(integral_coefficients=True)
            a^3 + a^2 - 3*a - 1

        TESTS::

            sage: K.<a> = NumberField(x^5-2)
            sage: K.random_element(-1)
            Traceback (most recent call last):
            ...
            TypeError: x must be < y
            sage: K.random_element(5,0)
            Traceback (most recent call last):
            ...
            TypeError: x must be < y
            sage: QQ[I].random_element(0)
            Traceback (most recent call last):
            ...
            TypeError: x must be > 0
        """
        if integral_coefficients:
            den_bound = 1

        return self._zero_element._random_element(num_bound=num_bound,
                                                  den_bound=den_bound,
                                                  distribution=distribution)

    def subfield(self, alpha, name=None, names=None):
        r"""
        Return a number field `K` isomorphic to `\QQ(\alpha)`
        (if this is an absolute number field) or `L(\alpha)` (if this
        is a relative extension `M/L`) and a map from K to self that
        sends the generator of K to alpha.

        INPUT:

        -  ``alpha`` - an element of self, or something that
           coerces to an element of self.

        OUTPUT:

        - ``K`` - a number field
        - ``from_K`` - a homomorphism from K to self that
          sends the generator of K to alpha.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 3); K
            Number Field in a with defining polynomial x^4 - 3
            sage: H.<b>, from_H = K.subfield(a^2)
            sage: H
            Number Field in b with defining polynomial x^2 - 3
            sage: from_H(b)
            a^2
            sage: from_H
            Ring morphism:
              From: Number Field in b with defining polynomial x^2 - 3
              To:   Number Field in a with defining polynomial x^4 - 3
              Defn: b |--> a^2

        A relative example. Note that the result returned is the subfield generated
        by `\alpha` over ``self.base_field()``, not over `\QQ` (see :trac:`5392`)::

            sage: L.<a> = NumberField(x^2 - 3)
            sage: M.<b> = L.extension(x^4 + 1)
            sage: K, phi = M.subfield(b^2)
            sage: K.base_field() is L
            True

        Subfields inherit embeddings::

            sage: K.<z> = CyclotomicField(5)
            sage: L, K_from_L = K.subfield(z-z^2-z^3+z^4)
            sage: L
            Number Field in z0 with defining polynomial x^2 - 5
            sage: CLF_from_K = K.coerce_embedding(); CLF_from_K
            Generic morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Lazy Field
              Defn: z -> 0.309016994374948? + 0.951056516295154?*I
            sage: CLF_from_L = L.coerce_embedding(); CLF_from_L
            Generic morphism:
              From: Number Field in z0 with defining polynomial x^2 - 5
              To:   Complex Lazy Field
              Defn: z0 -> 2.23606797749979? + 0.?e-14*I

        Check transitivity::

            sage: CLF_from_L(L.gen())
            2.23606797749979? + 0.?e-14*I
            sage: CLF_from_K(K_from_L(L.gen()))
            2.23606797749979? + 0.?e-14*I

        If `self` has no specified embedding, then `K` comes with an
        embedding in `self`::

            sage: K.<a> = NumberField(x^6 - 6*x^4 + 8*x^2 - 1)
            sage: L.<b>, from_L = K.subfield(a^2)
            sage: L
            Number Field in b with defining polynomial x^3 - 6*x^2 + 8*x - 1
            sage: L.gen_embedding()
            a^2

        You can also view a number field as having a different generator by
        just choosing the input to generate the whole field; for that it is
        better to use ``self.change_generator``, which gives
        isomorphisms in both directions.
        """
        if not names is None:
            name = names
        if name is None:
            name = self.variable_name() + '0'
        beta = self(alpha)
        f = beta.minpoly()
        # If self has a specified embedding, K should inherit it
        if self.coerce_embedding() is not None:
            emb = self.coerce_embedding()(beta)
        else:
            # Otherwise K should at least come with an embedding in self
            emb = beta
        K = NumberField(f, names=name, embedding=emb)
        from_K = K.hom([beta])
        return K, from_K

    def change_generator(self, alpha, name=None, names=None):
        r"""
        Given the number field self, construct another isomorphic number
        field `K` generated by the element alpha of self, along
        with isomorphisms from `K` to self and from self to
        `K`.

        EXAMPLES::

            sage: L.<i> = NumberField(x^2 + 1); L
            Number Field in i with defining polynomial x^2 + 1
            sage: K, from_K, to_K = L.change_generator(i/2 + 3)
            sage: K
            Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
            sage: from_K
            Ring morphism:
              From: Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
              To:   Number Field in i with defining polynomial x^2 + 1
              Defn: i0 |--> 1/2*i + 3
            sage: to_K
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
              Defn: i |--> 2*i0 - 6

        We can also do

        ::

            sage: K.<c>, from_K, to_K = L.change_generator(i/2 + 3); K
            Number Field in c with defining polynomial x^2 - 6*x + 37/4


        We compute the image of the generator `\sqrt{-1}` of `L`.

        ::

            sage: to_K(i)
            2*c - 6

        Note that the image is indeed a square root of -1.

        ::

            sage: to_K(i)^2
            -1
            sage: from_K(to_K(i))
            i
            sage: to_K(from_K(c))
            c
        """
        if not names is None:
            name = names
        alpha = self(alpha)
        K, from_K = self.subfield(alpha, name=name)
        if K.degree() != self.degree():
            raise ValueError, "alpha must generate a field of degree %s, but alpha generates a subfield of degree %s"%(self.degree(), K.degree())
        # Now compute to_K, which is an isomorphism
        # from self to K such that from_K(to_K(x)) == x for all x,
        # and to_K(from_K(y)) == y.
        # To do this, we must compute the image of self.gen()
        # under to_K.   This means writing self.gen() as a
        # polynomial in alpha, which is possible by the degree
        # check above.  This latter we do by linear algebra.
        phi = alpha.coordinates_in_terms_of_powers()
        c = phi(self.gen())
        to_K = self.hom([K(c)])
        return K, from_K, to_K

    def is_absolute(self):
        """
        Returns True if self is an absolute field.

        This function will be implemented in the derived classes.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.is_absolute()
            True
        """
        raise NotImplementedError

    def is_relative(self):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^10 - 2)
            sage: K.is_absolute()
            True
            sage: K.is_relative()
            False
        """
        return not self.is_absolute()

    def absolute_field(self, names):
        """
        Returns self as an absolute extension over QQ.

        OUTPUT:


        -  ``K`` - this number field (since it is already
           absolute)


        Also, ``K.structure()`` returns from_K and to_K,
        where from_K is an isomorphism from K to self and to_K is an
        isomorphism from self to K.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.absolute_field('a')
            Number Field in a with defining polynomial x^4 + x^3 + x^2 + x + 1
        """
        try:
            return self.__absolute_field[names]
        except KeyError:
            pass
        except AttributeError:
            self.__absolute_field = {}
        K = NumberField(self.defining_polynomial(), names, cache=False)
        K._set_structure(maps.NameChangeMap(K, self), maps.NameChangeMap(self, K))
        self.__absolute_field[names] = K
        return K

    def is_isomorphic(self, other):
        """
        Return True if self is isomorphic as a number field to other.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: m.<b> = NumberField(x^2 + 4)
            sage: k.is_isomorphic(m)
            True
            sage: m.<b> = NumberField(x^2 + 5)
            sage: k.is_isomorphic (m)
            False

        ::

            sage: k = NumberField(x^3 + 2, 'a')
            sage: k.is_isomorphic(NumberField((x+1/3)^3 + 2, 'b'))
            True
            sage: k.is_isomorphic(NumberField(x^3 + 4, 'b'))
            True
            sage: k.is_isomorphic(NumberField(x^3 + 5, 'b'))
            False
        """
        if not isinstance(other, NumberField_generic):
            raise ValueError, "other must be a generic number field."
        t = self.pari_polynomial().nfisisom(other.pari_polynomial())
        return t != 0

    def is_totally_real(self):
        """
        Return True if self is totally real, and False otherwise.

        Totally real means that every isomorphic embedding of self into the
        complex numbers has image contained in the real numbers.

        EXAMPLES::

            sage: NumberField(x^2+2, 'alpha').is_totally_real()
            False
            sage: NumberField(x^2-2, 'alpha').is_totally_real()
            True
            sage: NumberField(x^4-2, 'alpha').is_totally_real()
            False
        """
        return self.signature()[1] == 0

    def is_totally_imaginary(self):
        """
        Return True if self is totally imaginary, and False otherwise.

        Totally imaginary means that no isomorphic embedding of self into
        the complex numbers has image contained in the real numbers.

        EXAMPLES::

            sage: NumberField(x^2+2, 'alpha').is_totally_imaginary()
            True
            sage: NumberField(x^2-2, 'alpha').is_totally_imaginary()
            False
            sage: NumberField(x^4-2, 'alpha').is_totally_imaginary()
            False
        """
        return self.signature()[0] == 0

    def is_CM(self):
        r"""
        Return True if self is a CM field (i.e. a totally imaginary
        quadratic extension of a totally real field).

        EXAMPLES::

            sage: Q.<a> = NumberField(x - 1)
            sage: Q.is_CM()
            False
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.is_CM()
            True
            sage: L.<zeta20> = CyclotomicField(20)
            sage: L.is_CM()
            True
            sage: K.<omega> = QuadraticField(-3)
            sage: K.is_CM()
            True
            sage: L.<sqrt5> = QuadraticField(5)
            sage: L.is_CM()
            False
            sage: F.<a> = NumberField(x^3 - 2)
            sage: F.is_CM()
            False
            sage: F.<a> = NumberField(x^4-x^3-3*x^2+x+1)
            sage: F.is_CM()
            False

        The following are non-CM totally imaginary fields.

        ::

            sage: F.<a> = NumberField(x^4 + x^3 - x^2 - x + 1)
            sage: F.is_totally_imaginary()
            True
            sage: F.is_CM()
            False
            sage: F2.<a> = NumberField(x^12 - 5*x^11 + 8*x^10 - 5*x^9 - \
                                       x^8 + 9*x^7 + 7*x^6 - 3*x^5 + 5*x^4 + \
                                       7*x^3 - 4*x^2 - 7*x + 7)
            sage: F2.is_totally_imaginary()
            True
            sage: F2.is_CM()
            False

        The following is a non-cyclotomic CM field.

        ::

            sage: M.<a> = NumberField(x^4 - x^3 - x^2 - 2*x + 4)
            sage: M.is_CM()
            True

        Now, we construct a totally imaginary quadratic extension of a
        totally real field (which is not cyclotomic).

        ::

            sage: E_0.<a> = NumberField(x^7 - 4*x^6 - 4*x^5 + 10*x^4 + 4*x^3 - \
                                        6*x^2 - x + 1)
            sage: E_0.is_totally_real()
            True
            sage: E.<b> = E_0.extension(x^2 + 1)
            sage: E.is_CM()
            True

        Finally, a CM field that is given as an extension that is not CM.

        ::

            sage: E_0.<a> = NumberField(x^2 - 4*x + 16)
            sage: y = polygen(E_0)
            sage: E.<z> = E_0.extension(y^2 - E_0.gen() / 2)
            sage: E.is_CM()
            True
            sage: E.is_CM_extension()
            False

        """

        #Return cached answer if available
        try:
            return self.__is_CM
        except(AttributeError):
            pass

        #Then, deal with simple cases
        if is_odd(self.absolute_degree()):
            self.__is_CM = False
            return False
        if isinstance(
           self, sage.rings.number_field.number_field.NumberField_quadratic):
            self.__is_CM = (self.discriminant() < 0)
            return self.__is_CM
        if isinstance(
           self, sage.rings.number_field.number_field.NumberField_cyclotomic):
            self.__is_CM = True
            return True
        if not self.is_totally_imaginary():
            self.__is_CM = False
            return False
        if self.is_absolute():
            K = self
        else:
            F = self.base_field()
            if F.absolute_degree() == self.absolute_degree() / 2:
                if F.is_totally_real():
                    self.__is_CM = True
                    self.__max_tot_real_sub = [F, self.coerce_map_from(F)]
                    return True
            K = self.absolute_field('z')

        #Check for index 2 subextensions that are totally real
        possibilities = K.subfields(K.absolute_degree()/2)
        for F, phi, _ in possibilities:
            if F.is_totally_real():
                self.__is_CM = True
                if self.is_relative():
                    phi = phi.post_compose(K.structure()[0])
                self.__max_tot_real_sub = [F, phi]
                return True
        self.__is_CM = False
        return False

    def complex_conjugation(self):
        """
        Return the complex conjugation of self.

        This is only well-defined for fields contained in CM fields
        (i.e. for totally real fields and CM fields). Recall that a CM
        field is a totally imaginary quadratic extension of a totally
        real field. For other fields, a ValueError is raised.

        EXAMPLES::

            sage: QuadraticField(-1, 'I').complex_conjugation()
            Ring endomorphism of Number Field in I with defining polynomial x^2 + 1
              Defn: I |--> -I
            sage: CyclotomicField(8).complex_conjugation()
            Ring endomorphism of Cyclotomic Field of order 8 and degree 4
              Defn: zeta8 |--> -zeta8^3
            sage: QuadraticField(5, 'a').complex_conjugation()
            Identity endomorphism of Number Field in a with defining polynomial x^2 - 5
            sage: F = NumberField(x^4 + x^3 - 3*x^2 - x + 1, 'a')
            sage: F.is_totally_real()
            True
            sage: F.complex_conjugation()
            Identity endomorphism of Number Field in a with defining polynomial x^4 + x^3 - 3*x^2 - x + 1
            sage: F.<b> = NumberField(x^2 - 2)
            sage: F.extension(x^2 + 1, 'a').complex_conjugation()
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + 1 over its base field
              Defn: a |--> -a
                    b |--> b
            sage: F2.<b> = NumberField(x^2 + 2)
            sage: K2.<a> = F2.extension(x^2 + 1)
            sage: cc = K2.complex_conjugation()
            sage: cc(a)
            -a
            sage: cc(b)
            -b

        """

        #Return cached answer if available
        try:
            return self.__complex_conjugation
        except(AttributeError):
            pass

        #Then, deal with simple cases
        if isinstance(
           self, sage.rings.number_field.number_field.NumberField_quadratic):
            disc = self.discriminant()
            if disc > 0:
                self.__complex_conjugation = self.coerce_map_from(self)
                return self.__complex_conjugation
            else:
                a = self.gen()
                r = a.trace()
                iy = a - r / 2
                self.__complex_conjugation = self.hom([a - 2 * iy], check=False)
            return self.__complex_conjugation
        if isinstance(
           self, sage.rings.number_field.number_field.NumberField_cyclotomic):
            zeta = self.gen()
            self.__complex_conjugation = self.hom([zeta ** (-1)], check=False)
            return self.__complex_conjugation
        if self.is_totally_real():
            self.__complex_conjugation = self.coerce_map_from(self)
            return self.__complex_conjugation

        if not self.is_CM():
            raise ValueError('Complex conjugation is only well-defined for fields contained in CM fields.')

        #In the remaining case, self.is_CM() should have cached __max_tot_real_sub
        try:
            F, phi = self.__max_tot_real_sub
        except(AttributeError):
            F, phi = self.maximal_totally_real_subfield()
        if self.is_absolute():
            K_rel = self.relativize(phi, self.variable_name() * 2)
            to_abs, from_abs = K_rel.structure()
            self.__complex_conjugation = K_rel.automorphisms()[1].pre_compose( \
               from_abs).post_compose(to_abs)
            self.__complex_conjugation = self.hom([self.__complex_conjugation(self.gen())], check=False)
            return self.__complex_conjugation
        else:
            if self.is_CM_extension():
                return self.automorphisms()[1]
            K_abs = self.absolute_field(self.variable_name() * 2)
            to_self, from_self = K_abs.structure()
            K_rel = K_abs.relativize(phi.post_compose(from_self), self.variable_name() * 3)
            to_abs, from_abs = K_rel.structure()
            self.__complex_conjugation = K_rel.automorphisms()[1].pre_compose(from_abs).post_compose(to_abs)
            self.__complex_conjugation = K_abs.hom([self.__complex_conjugation(K_abs.gen())], check=False)
            self.__complex_conjugation = self.__complex_conjugation.pre_compose(from_self).post_compose(to_self)
            return self.__complex_conjugation

    def maximal_totally_real_subfield(self):
        """
        Return the maximal totally real subfield of self together with an embedding of it into self.

        EXAMPLES::

            sage: F.<a> = QuadraticField(11)
            sage: F.maximal_totally_real_subfield()
            [Number Field in a with defining polynomial x^2 - 11, Identity endomorphism of Number Field in a with defining polynomial x^2 - 11]
            sage: F.<a> = QuadraticField(-15)
            sage: F.maximal_totally_real_subfield()
            [Rational Field, Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 + 15]
            sage: F.<a> = CyclotomicField(29)
            sage: F.maximal_totally_real_subfield()
            (Number Field in a0 with defining polynomial x^14 + x^13 - 13*x^12 - 12*x^11 + 66*x^10 + 55*x^9 - 165*x^8 - 120*x^7 + 210*x^6 + 126*x^5 - 126*x^4 - 56*x^3 + 28*x^2 + 7*x - 1, Ring morphism:
              From: Number Field in a0 with defining polynomial x^14 + x^13 - 13*x^12 - 12*x^11 + 66*x^10 + 55*x^9 - 165*x^8 - 120*x^7 + 210*x^6 + 126*x^5 - 126*x^4 - 56*x^3 + 28*x^2 + 7*x - 1
              To:   Cyclotomic Field of order 29 and degree 28
              Defn: a0 |--> -a^27 - a^26 - a^25 - a^24 - a^23 - a^22 - a^21 - a^20 - a^19 - a^18 - a^17 - a^16 - a^15 - a^14 - a^13 - a^12 - a^11 - a^10 - a^9 - a^8 - a^7 - a^6 - a^5 - a^4 - a^3 - a^2 - 1)
            sage: F.<a> = NumberField(x^3 - 2)
            sage: F.maximal_totally_real_subfield()
            [Rational Field, Conversion map:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^3 - 2]
            sage: F.<a> = NumberField(x^4 - x^3 - x^2 + x + 1)
            sage: F.maximal_totally_real_subfield()
            [Rational Field, Conversion map:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^4 - x^3 - x^2 + x + 1]
            sage: F.<a> = NumberField(x^4 - x^3 + 2*x^2 + x + 1)
            sage: F.maximal_totally_real_subfield()
            [Number Field in a1 with defining polynomial x^2 - x - 1, Ring morphism:
              From: Number Field in a1 with defining polynomial x^2 - x - 1
              To:   Number Field in a with defining polynomial x^4 - x^3 + 2*x^2 + x + 1
              Defn: a1 |--> -1/2*a^3 - 1/2]
            sage: F.<a> = NumberField(x^4-4*x^2-x+1)
            sage: F.maximal_totally_real_subfield()
            [Number Field in a with defining polynomial x^4 - 4*x^2 - x + 1, Identity endomorphism of Number Field in a with defining polynomial x^4 - 4*x^2 - x + 1]

        An example of a relative extension where the base field is not the maximal totally real subfield.

        ::

            sage: E_0.<a> = NumberField(x^2 - 4*x + 16)
            sage: y = polygen(E_0)
            sage: E.<z> = E_0.extension(y^2 - E_0.gen() / 2)
            sage: E.maximal_totally_real_subfield()
            [Number Field in z2 with defining polynomial x^2 - 6, Composite map:
              From: Number Field in z2 with defining polynomial x^2 - 6
              To:   Number Field in z with defining polynomial x^2 - 1/2*a over its base field
              Defn:   Ring morphism:
                      From: Number Field in z2 with defining polynomial x^2 - 6
                      To:   Number Field in z with defining polynomial x^4 - 2*x^2 + 4
                      Defn: z2 |--> -1/2*z^3 + 2*z
                    then
                      Isomorphism map:
                      From: Number Field in z with defining polynomial x^4 - 2*x^2 + 4
                      To:   Number Field in z with defining polynomial x^2 - 1/2*a over its base field]

        """

        try:
            return self.__max_tot_real_sub
        except(AttributeError):
            pass

        if isinstance(
           self, sage.rings.number_field.number_field.NumberField_quadratic):
            if self.discriminant() > 0:
                self.__max_tot_real_sub = [self, self.coerce_map_from(self)]
                return self.__max_tot_real_sub
            else:
                self.__max_tot_real_sub = [QQ, self.coerce_map_from(QQ)]
            return self.__max_tot_real_sub
        if isinstance(
           self, sage.rings.number_field.number_field.NumberField_cyclotomic):
            zeta = self.gen()
            self.__max_tot_real_sub = self.subfield(zeta + zeta ** (-1))
            return self.__max_tot_real_sub
        if self.is_totally_real():
            self.__max_tot_real_sub = [self, self.coerce_map_from(self)]
            return self.__max_tot_real_sub
        if self.is_absolute():
            K = self
        else:
            if self.is_CM_extension():
                self.__max_tot_real_sub = [self.base_field(), self.coerce_map_from(self.base_field())]
                return self.__max_tot_real_sub
            K = self.absolute_field('z')

        d = K.absolute_degree()
        divs = d.divisors()[1:-1]
        divs.reverse()
        for i in divs:
            possibilities = K.subfields(i)
            for F, phi, _ in possibilities:
                if F.is_totally_real():
                    if self.is_relative():
                        phi = phi.post_compose(K.structure()[0])
                    self.__max_tot_real_sub = [F, phi]
                    return self.__max_tot_real_sub
        self.__max_tot_real_sub = [QQ, self.coerce_map_from(QQ)]
        return self.__max_tot_real_sub

    def complex_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this number field into the approximate
        complex field with precision prec.

        This always embeds into an MPFR based complex field.  If you
        want embeddings into the 53-bit double precision, which is
        faster, use ``self.embeddings(CDF)``.

        EXAMPLES::

            sage: k.<a> = NumberField(x^5 + x + 17)
            sage: v = k.complex_embeddings()
            sage: ls = [phi(k.0^2) for phi in v] ; ls # random order
            [2.97572074038...,
             -2.40889943716 + 1.90254105304*I,
             -2.40889943716 - 1.90254105304*I,
             0.921039066973 + 3.07553311885*I,
             0.921039066973 - 3.07553311885*I]
            sage: K.<a> = NumberField(x^3 + 2)
            sage: ls = K.complex_embeddings() ; ls # random order
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Complex Double Field
              Defn: a |--> -1.25992104989...,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Complex Double Field
              Defn: a |--> 0.629960524947 - 1.09112363597*I,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Complex Double Field
              Defn: a |--> 0.629960524947 + 1.09112363597*I
            ]
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        return self.embeddings(CC)

    def real_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this number field into the approximate
        real field with precision prec.

        If prec is 53 (the default), then the real double field is
        used; otherwise the arbitrary precision (but slow) real field
        is used.  If you want embeddings into the 53-bit double
        precision, which is faster, use ``self.embeddings(RDF)``.

        .. NOTE::

            This function uses finite precision real numbers.
            In functions that should output proven results, one
            could use ``self.embeddings(AA)`` instead.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.real_embeddings()
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Real Field with 53 bits of precision
              Defn: a |--> -1.25992104989487
            ]
            sage: K.real_embeddings(16)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Real Field with 16 bits of precision
              Defn: a |--> -1.260
            ]
            sage: K.real_embeddings(100)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Real Field with 100 bits of precision
              Defn: a |--> -1.2599210498948731647672106073
            ]

        As this is a numerical function, the number of embeddings
        may be incorrect if the precision is too low::

            sage: K = NumberField(x^2+2*10^1000*x + 10^2000+1, 'a')
            sage: len(K.real_embeddings())
            2
            sage: len(K.real_embeddings(100))
            2
            sage: len(K.real_embeddings(10000))
            0
            sage: len(K.embeddings(AA))
            0

        """
        K = sage.rings.real_mpfr.RealField(prec)
        return self.embeddings(K)

    def specified_complex_embedding(self):
        r"""
        Returns the embedding of this field into the complex numbers which has
        been specified.

        Fields created with the ``QuadraticField`` or
        ``CyclotomicField`` constructors come with an implicit
        embedding. To get one of these fields without the embedding, use
        the generic ``NumberField`` constructor.

        EXAMPLES::

            sage: QuadraticField(-1, 'I').specified_complex_embedding()
            Generic morphism:
              From: Number Field in I with defining polynomial x^2 + 1
              To:   Complex Lazy Field
              Defn: I -> 1*I

        ::

            sage: QuadraticField(3, 'a').specified_complex_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^2 - 3
              To:   Real Lazy Field
              Defn: a -> 1.732050807568878?

        ::

            sage: CyclotomicField(13).specified_complex_embedding()
            Generic morphism:
              From: Cyclotomic Field of order 13 and degree 12
              To:   Complex Lazy Field
              Defn: zeta13 -> 0.885456025653210? + 0.464723172043769?*I

        Most fields don't implicitly have embeddings unless explicitly
        specified::

            sage: NumberField(x^2-2, 'a').specified_complex_embedding() is None
            True
            sage: NumberField(x^3-x+5, 'a').specified_complex_embedding() is None
            True
            sage: NumberField(x^3-x+5, 'a', embedding=2).specified_complex_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - x + 5
              To:   Real Lazy Field
              Defn: a -> -1.904160859134921?
            sage: NumberField(x^3-x+5, 'a', embedding=CDF.0).specified_complex_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - x + 5
              To:   Complex Lazy Field
              Defn: a -> 0.952080429567461? + 1.311248044077123?*I

        This function only returns complex embeddings::

            sage: K.<a> = NumberField(x^2-2, embedding=Qp(7)(2).sqrt())
            sage: K.specified_complex_embedding() is None
            True
            sage: K.gen_embedding()
            3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8 + 6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16 + 7^17 + 4*7^18 + 6*7^19 + O(7^20)
            sage: K.coerce_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^2 - 2
              To:   7-adic Field with capped relative precision 20
              Defn: a -> 3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8 + 6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16 + 7^17 + 4*7^18 + 6*7^19 + O(7^20)
        """
        embedding = self.coerce_embedding()
        if embedding is not None:
            from sage.rings.real_mpfr import mpfr_prec_min
            from sage.rings.complex_field import ComplexField
            if ComplexField(mpfr_prec_min()).has_coerce_map_from(embedding.codomain()):
                return embedding

    def gen_embedding(self):
        """
        If an embedding has been specified, return the image of the
        generator under that embedding. Otherwise return None.

        EXAMPLES::

            sage: QuadraticField(-7, 'a').gen_embedding()
            2.645751311064591?*I
            sage: NumberField(x^2+7, 'a').gen_embedding() # None
        """
        embedding = self.coerce_embedding()
        if embedding is None:
            return None
        else:
            return embedding(self.gen())

    def algebraic_closure(self):
        """
        Return the algebraic closure of self (which is QQbar).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: K.algebraic_closure()
            Algebraic Field
            sage: K.<a> = NumberField(x^3-2)
            sage: K.algebraic_closure()
            Algebraic Field
            sage: K = CyclotomicField(23)
            sage: K.algebraic_closure()
            Algebraic Field
        """
        return sage.rings.all.QQbar

    def latex_variable_name(self, name=None):
        """
        Return the latex representation of the variable name for this
        number field.

        EXAMPLES::

            sage: NumberField(x^2 + 3, 'a').latex_variable_name()
            'a'
            sage: NumberField(x^3 + 3, 'theta3').latex_variable_name()
            '\\theta_{3}'
            sage: CyclotomicField(5).latex_variable_name()
            '\\zeta_{5}'
        """
        if name is None:
            return self.__latex_variable_name
        else:
            self.__latex_variable_name = name

    def _repr_(self):
        """
        Return string representation of this number field.

        EXAMPLES::

            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._repr_()
            'Number Field in a with defining polynomial x^13 - 2/3*x + 3'
        """
        return "Number Field in %s with defining polynomial %s"%(
                   self.variable_name(), self.polynomial())

    def _latex_(self):
        r"""
        Return latex representation of this number field. This is viewed as
        a polynomial quotient ring over a field.

        EXAMPLES::

            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._latex_()
            '\\Bold{Q}[a]/(a^{13} - \\frac{2}{3} a + 3)'
            sage: latex(k)
            \Bold{Q}[a]/(a^{13} - \frac{2}{3} a + 3)

        Numbered variables are often correctly typeset::

            sage: k.<theta25> = NumberField(x^25+x+1)
            sage: print k._latex_()
            \Bold{Q}[\theta_{25}]/(\theta_{25}^{25} + \theta_{25} + 1)
        """
        return "%s[%s]/(%s)"%(latex(QQ), self.latex_variable_name(),
                              self.polynomial()._latex_(self.latex_variable_name()))

    def __cmp__(self, other):
        """
        Compare a number field with something else.

        INPUT:


        -  ``other`` - arbitrary Python object.


        If other is not a number field, then the types of self and other
        are compared. If both are number fields, then the variable names
        are compared. If those are the same, then the underlying defining
        polynomials are compared. If the polynomials are the same, the
        number fields are considered "equal", but need not be identical.
        Coercion between equal number fields is allowed.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 + 2); m.<b> = NumberField(x^3 + 2)
            sage: cmp(k,m)
            -1
            sage: cmp(m,k)
            1
            sage: k == QQ
            False
            sage: k.<a> = NumberField(x^3 + 2); m.<a> = NumberField(x^3 + 2)
            sage: k is m
            True
            sage: m = loads(dumps(k))
            sage: k is m
            False
            sage: k == m
            True

        TESTS::

            sage: x = QQ['x'].gen()
            sage: y = ZZ['y'].gen()
            sage: K = NumberField(x^3 + x + 3, 'a'); K
            Number Field in a with defining polynomial x^3 + x + 3
            sage: K.defining_polynomial().parent()
            Univariate Polynomial Ring in x over Rational Field

        ::

            sage: L = NumberField(y^3 + y + 3, 'a'); L
            Number Field in a with defining polynomial y^3 + y + 3
            sage: L.defining_polynomial().parent()
            Univariate Polynomial Ring in y over Rational Field
            sage: L == K
            True

        ::

            sage: NumberField(ZZ['x'].0^4 + 23, 'a') == NumberField(ZZ['y'].0^4 + 23, 'a')
            True
            sage: NumberField(ZZ['x'].0^4 + 23, 'a') == NumberField(QQ['y'].0^4 + 23, 'a')
            True
            sage: NumberField(QQ['x'].0^4 + 23, 'a') == NumberField(QQ['y'].0^4 + 23, 'a')
            True

        ::

            sage: x = var('x'); y = ZZ['y'].gen()
            sage: NumberField(x^3 + x + 5, 'a') == NumberField(y^3 + y + 5, 'a')
            True
            sage: NumberField(x^3 + x + 5, 'a') == NumberField(y^4 + y + 5, 'a')
            False
            sage: NumberField(x^3 + x + 5, 'a') == NumberField(x^3 + x + 5, 'b')
            False
            sage: QuadraticField(2, 'a', embedding=2) == QuadraticField(2, 'a', embedding=-2)
            False
        """
        if not isinstance(other, NumberField_generic):
            return cmp(type(self), type(other))
        c = cmp(self.variable_name(), other.variable_name())
        if c: return c
        # compare coefficients so that the polynomial variable does not count
        c = cmp(list(self.__polynomial), list(other.__polynomial))
        if c: return c
        # Now we compare the embeddings (if any).
        f, g = self.coerce_embedding(), other.coerce_embedding()
        if f is None and g is None:
            return 0
        elif f is None:
            return -1
        elif g is None:
            return 1
        else:
            return cmp(self.coerce_embedding()(self.gen()),
                       other.coerce_embedding()(other.gen()))

    def _ideal_class_(self, n=0):
        """
        Return the Python class used in defining the zero ideal of the ring
        of integers of this number field.

        This function is required by the general ring/ideal machinery. The
        value defined here is the default value for all number fields.

        EXAMPLES::

            sage: NumberField(x^2 + 2, 'c')._ideal_class_()
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldIdeal'>
        """
        return sage.rings.number_field.number_field_ideal.NumberFieldIdeal

    def _fractional_ideal_class_(self):
        """
        Return the Python class used in defining fractional ideals of the
        ring of integers of this number field.

        This function is required by the general ring/ideal machinery. The
        value defined here is the default value for all number fields
        *except* relative number fields; this function is overridden by
        one of the same name on class NumberField_relative.

        EXAMPLES::

            sage: NumberField(x^2 + 2, 'c')._fractional_ideal_class_()
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
        """
        return sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal

    def ideal(self, *gens, **kwds):
        """
        K.ideal() returns a fractional ideal of the field, except for the
        zero ideal which is not a fractional ideal.

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: K.ideal(2)
            Fractional ideal (2)
            sage: K.ideal(2+i)
            Fractional ideal (i + 2)
            sage: K.ideal(0)
            Ideal (0) of Number Field in i with defining polynomial x^2 + 1
        """
        try:
            return self.fractional_ideal(*gens, **kwds)
        except ValueError:
            return sage.rings.ring.Ring.ideal(self, gens, **kwds)

    def fractional_ideal(self, *gens, **kwds):
        r"""
        Return the ideal in `\mathcal{O}_K` generated by gens.
        This overrides the ``sage.rings.ring.Field`` method to
        use the ``sage.rings.ring.Ring`` one instead, since
        we're not really concerned with ideals in a field but in its ring
        of integers.

        INPUT:


        -  ``gens`` - a list of generators, or a number field
           ideal.


        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: K.fractional_ideal([1/a])
            Fractional ideal (1/2*a^2)

        One can also input a number field ideal itself,
        or, more usefully, for a tower of number fields an ideal
        in one of the fields lower down the tower.

        ::

            sage: K.fractional_ideal(K.ideal(a))
            Fractional ideal (a)
            sage: L.<b> = K.extension(x^2 - 3, x^2 + 1)
            sage: M.<c> = L.extension(x^2 + 1)
            sage: L.ideal(K.ideal(2, a))
            Fractional ideal (a)
            sage: M.ideal(K.ideal(2, a)) == M.ideal(a*(b - c)/2)
            True

        The zero ideal is not a fractional ideal!

        ::

            sage: K.fractional_ideal(0)
            Traceback (most recent call last):
            ...
            ValueError: gens must have a nonzero element (zero ideal is not a fractional ideal)
        """
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        if len(gens) == 1 and isinstance(gens[0], NumberFieldFractionalIdeal):
            I = gens[0]
            if I.number_field() is self:
                return I
            else:
                gens = I.gens()
        return self._fractional_ideal_class_()(self, gens, **kwds)

    def ideals_of_bdd_norm(self, bound):
        """
        All integral ideals of bounded norm.

        INPUT:


        -  ``bound`` - a positive integer


        OUTPUT: A dict of all integral ideals I such that Norm(I) <= bound,
        keyed by norm.

        EXAMPLE::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: d = K.ideals_of_bdd_norm(10)
            sage: for n in d:
            ...       print n
            ...       for I in d[n]:
            ...           print I
            1
            Fractional ideal (1)
            2
            Fractional ideal (2, 1/2*a - 1/2)
            Fractional ideal (2, 1/2*a + 1/2)
            3
            Fractional ideal (3, 1/2*a - 1/2)
            Fractional ideal (3, 1/2*a + 1/2)
            4
            Fractional ideal (4, 1/2*a + 3/2)
            Fractional ideal (2)
            Fractional ideal (4, 1/2*a + 5/2)
            5
            6
            Fractional ideal (-1/2*a + 1/2)
            Fractional ideal (6, 1/2*a + 5/2)
            Fractional ideal (6, 1/2*a + 7/2)
            Fractional ideal (1/2*a + 1/2)
            7
            8
            Fractional ideal (-1/2*a - 3/2)
            Fractional ideal (4, a - 1)
            Fractional ideal (4, a + 1)
            Fractional ideal (1/2*a - 3/2)
            9
            Fractional ideal (9, 1/2*a + 11/2)
            Fractional ideal (3)
            Fractional ideal (9, 1/2*a + 7/2)
            10
        """
        hnf_ideals = self.pari_nf().ideallist(bound)
        d = {}
        for i in xrange(bound):
            d[i+1] = [self.ideal(hnf) for hnf in hnf_ideals[i]]
        return d

    def primes_above(self, x, degree=None):
        r"""
        Return prime ideals of self lying over x.

        INPUT:


        -  ``x``: usually an element or ideal of self. It
           should be such that self.ideal(x) is sensible. This excludes x=0.

        -  ``degree`` (default: None): None or an integer.
           If None, find all primes above x of any degree. If an integer, find
           all primes above x such that the resulting residue field has
           exactly this degree.


        OUTPUT: A list of prime ideals of self lying over x. If degree
        is specified and no such ideal exists, returns the empty list.
        The output is sorted by residue degree first, then by
        underlying prime (or equivalently, by norm).

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

        ::

            sage: P2s = F.primes_above(2)
            sage: P2s # random
            [Fractional ideal (-t)]
            sage: all(2 in P2 for P2 in P2s)
            True
            sage: all(P2.is_prime() for P2 in P2s)
            True
            sage: [ P2.norm() for P2 in P2s ]
            [2]

        ::

            sage: P3s = F.primes_above(3)
            sage: P3s # random
            [Fractional ideal (t + 1)]
            sage: all(3 in P3 for P3 in P3s)
            True
            sage: all(P3.is_prime() for P3 in P3s)
            True
            sage: [ P3.norm() for P3 in P3s ]
            [3]

        The ideal (3) is totally ramified in F, so there is no degree 2
        prime above 3::

            sage: F.primes_above(3, degree=2)
            []
            sage: [ id.residue_class_degree() for id, _ in F.ideal(3).factor() ]
            [1]

        Asking for a specific degree works::

            sage: P5_1s = F.primes_above(5, degree=1)
            sage: P5_1s # random
            [Fractional ideal (-t^2 - 1)]
            sage: P5_1 = P5_1s[0]; P5_1.residue_class_degree()
            1

        ::

            sage: P5_2s = F.primes_above(5, degree=2)
            sage: P5_2s # random
            [Fractional ideal (t^2 - 2*t - 1)]
            sage: P5_2 = P5_2s[0]; P5_2.residue_class_degree()
            2

        Works in relative extensions too::

            sage: PQ.<X> = QQ[]
            sage: F.<a, b> = NumberField([X^2 - 2, X^2 - 3])
            sage: PF.<Y> = F[]
            sage: K.<c> = F.extension(Y^2 - (1 + a)*(a + b)*a*b)
            sage: I = F.ideal(a + 2*b)
            sage: P, Q = K.primes_above(I)
            sage: K.ideal(I) == P^4*Q
            True
            sage: K.primes_above(I, degree=1) == [P]
            True
            sage: K.primes_above(I, degree=4) == [Q]
            True

        It doesn't make sense to factor the ideal (0), so this raises an error::

            sage: F.prime_above(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'prime_factors'
        """
        if degree is not None:
            degree = ZZ(degree)
        facs = [ (id.residue_class_degree(), id.absolute_norm(), id) for id in self.prime_factors(x) ]
        facs.sort() # sorts on residue_class_degree(), lowest first
        if degree is None:
            return [ id for d, n, id in facs ]
        else:
            return [ id for d, n, id in facs if d == degree ]

    def prime_above(self, x, degree=None):
        r"""
        Return a prime ideal of self lying over x.

        INPUT:


        -  ``x``: usually an element or ideal of self. It
           should be such that self.ideal(x) is sensible. This excludes x=0.

        -  ``degree`` (default: None): None or an integer.
           If one, find a prime above x of any degree. If an integer, find a
           prime above x such that the resulting residue field has exactly
           this degree.


        OUTPUT: A prime ideal of self lying over x. If degree is specified
        and no such ideal exists, raises a ValueError.

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

        ::

            sage: P2 = F.prime_above(2)
            sage: P2 # random
            Fractional ideal (-t)
            sage: 2 in P2
            True
            sage: P2.is_prime()
            True
            sage: P2.norm()
            2

        ::

            sage: P3 = F.prime_above(3)
            sage: P3 # random
            Fractional ideal (t + 1)
            sage: 3 in P3
            True
            sage: P3.is_prime()
            True
            sage: P3.norm()
            3

        The ideal (3) is totally ramified in F, so there is no degree 2
        prime above 3::

            sage: F.prime_above(3, degree=2)
            Traceback (most recent call last):
            ...
            ValueError: No prime of degree 2 above Fractional ideal (3)
            sage: [ id.residue_class_degree() for id, _ in F.ideal(3).factor() ]
            [1]

        Asking for a specific degree works::

            sage: P5_1 = F.prime_above(5, degree=1)
            sage: P5_1 # random
            Fractional ideal (-t^2 - 1)
            sage: P5_1.residue_class_degree()
            1

        ::

            sage: P5_2 = F.prime_above(5, degree=2)
            sage: P5_2 # random
            Fractional ideal (t^2 - 2*t - 1)
            sage: P5_2.residue_class_degree()
            2

        Relative number fields are ok::

            sage: G = F.extension(x^2 - 11, 'b')
            sage: G.prime_above(7)
            Fractional ideal (b + 2)

        It doesn't make sense to factor the ideal (0)::

            sage: F.prime_above(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'prime_factors'

        """
        ids = self.primes_above(x, degree)
        if not ids:
            raise ValueError, "No prime of degree %s above %s" % (degree, self.ideal(x))
        return ids[0]

    def primes_of_degree_one_iter(self, num_integer_primes=10000, max_iterations=100):
        r"""
        Return an iterator yielding prime ideals of absolute degree one and
        small norm.

        .. warning::

           It is possible that there are no primes of `K` of
           absolute degree one of small prime norm, and it possible
           that this algorithm will not find any primes of small norm.

           See module :mod:`sage.rings.number_field.small_primes_of_degree_one`
           for details.

        INPUT:


        -  ``num_integer_primes (default: 10000)`` - an
           integer. We try to find primes of absolute norm no greater than the
           num_integer_primes-th prime number. For example, if
           num_integer_primes is 2, the largest norm found will be 3, since
           the second prime is 3.

        -  ``max_iterations (default: 100)`` - an integer. We
           test max_iterations integers to find small primes before raising
           StopIteration.


        EXAMPLES::

            sage: K.<z> = CyclotomicField(10)
            sage: it = K.primes_of_degree_one_iter()
            sage: Ps = [ it.next() for i in range(3) ]
            sage: Ps # random
            [Fractional ideal (z^3 + z + 1), Fractional ideal (3*z^3 - z^2 + z - 1), Fractional ideal (2*z^3 - 3*z^2 + z - 2)]
            sage: [ P.norm() for P in Ps ] # random
            [11, 31, 41]
            sage: [ P.residue_class_degree() for P in Ps ]
            [1, 1, 1]
        """
        from sage.rings.number_field.small_primes_of_degree_one import Small_primes_of_degree_one_iter
        return Small_primes_of_degree_one_iter(self, num_integer_primes, max_iterations)

    def primes_of_degree_one_list(self, n, num_integer_primes=10000, max_iterations=100):
        r"""
        Return a list of n prime ideals of absolute degree one and small
        norm.

        .. warning::

           It is possible that there are no primes of `K` of
           absolute degree one of small prime norm, and it possible
           that this algorithm will not find any primes of small norm.

           See module :mod:`sage.rings.number_field.small_primes_of_degree_one`
           for details.

        INPUT:


        -  ``num_integer_primes (default: 10000)`` - an
           integer. We try to find primes of absolute norm no greater than the
           num_integer_primes-th prime number. For example, if
           num_integer_primes is 2, the largest norm found will be 3, since
           the second prime is 3.

        -  ``max_iterations (default: 100)`` - an integer. We
           test max_iterations integers to find small primes before raising
           StopIteration.


        EXAMPLES::

            sage: K.<z> = CyclotomicField(10)
            sage: Ps = K.primes_of_degree_one_list(3)
            sage: Ps  # random output
            [Fractional ideal (-z^3 - z^2 + 1), Fractional ideal (2*z^3 - 2*z^2 + 2*z - 3), Fractional ideal (2*z^3 - 3*z^2 + z - 2)]
            sage: [ P.norm() for P in Ps ]
            [11, 31, 41]
            sage: [ P.residue_class_degree() for P in Ps ]
            [1, 1, 1]
        """
        it = self.primes_of_degree_one_iter()
        return [ it.next() for i in range(n) ]

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        Return whether or not there is a homomorphism defined by the given
        images of generators.

        To do this we just check that the elements of the image of the
        given generator (im_gens always has length 1) satisfies the
        relation of the defining poly of this field.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 - 3)
            sage: k._is_valid_homomorphism_(QQ, [0])
            False
            sage: k._is_valid_homomorphism_(k, [])
            False
            sage: k._is_valid_homomorphism_(k, [a])
            True
            sage: k._is_valid_homomorphism_(k, [-a])
            True
            sage: k._is_valid_homomorphism_(k, [a+1])
            False
        """
        try:
            if len(im_gens) != 1:
                return False
            # We need that elements of the base ring of the polynomial
            # ring map canonically into codomain.
            codomain._coerce_(rational.Rational(1))
            f = self.defining_polynomial()
            return codomain(f(im_gens[0])) == 0
        except (TypeError, ValueError):
            return False

    def pari_polynomial(self, name='x'):
        """
        PARI polynomial with integer coefficients corresponding to the
        polynomial that defines this number field.

        By default, this is a polynomial in the variable "x".  PARI
        prefers integral polynomials, so we clear the denominator.
        Therefore, this is NOT the same as simply converting the defining
        polynomial to PARI.

        EXAMPLES::

            sage: y = polygen(QQ)
            sage: k.<a> = NumberField(y^2 - 3/2*y + 5/3)
            sage: k.pari_polynomial()
            6*x^2 - 9*x + 10
            sage: k.polynomial()._pari_()
            x^2 - 3/2*x + 5/3
            sage: k.pari_polynomial('a')
            6*a^2 - 9*a + 10

        This fails with arguments which are not a valid PARI variable name::

            sage: k = QuadraticField(-1)
            sage: k.pari_polynomial('I')
            Traceback (most recent call last):
            ...
            PariError: I already exists with incompatible valence
            sage: k.pari_polynomial('i')
            i^2 + 1
            sage: k.pari_polynomial('theta')
            Traceback (most recent call last):
            ...
            PariError: theta already exists with incompatible valence
        """
        try:
            return self.__pari_polynomial.change_variable_name(name)
        except AttributeError:
            polypari = self.polynomial()._pari_with_name(name)
            polypari /= polypari.content()   # make polypari integral
            self.__pari_polynomial = polypari
            return self.__pari_polynomial

    def pari_nf(self, important=True):
        """
        PARI number field corresponding to this field.

        This is the number field constructed using nfinit(). This is the same
        as the number field got by doing pari(self) or gp(self).

        INPUT:

        - ``important`` -- (default: True) bool.  If False, raise a
          ``RuntimeError`` if we need to do a difficult discriminant
          factorization.  Useful when the PARI nf structure is useful
          but not strictly required, such as for factoring polynomials
          over this number field.

        EXAMPLES::

            sage: k.<a> = NumberField(x^4 - 3*x + 7); k
            Number Field in a with defining polynomial x^4 - 3*x + 7
            sage: k.pari_nf()[:4]
            [y^4 - 3*y + 7, [0, 2], 85621, 1]
            sage: pari(k)[:4]
            [y^4 - 3*y + 7, [0, 2], 85621, 1]

        ::

            sage: k.<a> = NumberField(x^4 - 3/2*x + 5/3); k
            Number Field in a with defining polynomial x^4 - 3/2*x + 5/3
            sage: k.pari_nf()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce number field defined by non-integral polynomial to PARI.
            sage: pari(k)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce number field defined by non-integral polynomial to PARI.
            sage: gp(k)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce number field defined by non-integral polynomial to PARI.

        With ``important=False``, we simply bail out if we cannot
        easily factor the discriminant::

            sage: p = next_prime(10^40); q = next_prime(10^41)
            sage: K.<a> = NumberField(x^2 - p*q)
            sage: K.pari_nf(important=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Unable to factor discriminant with trial division

        Next, we illustrate the ``maximize_at_primes`` and ``assume_disc_small``
        parameters of the NumberField contructor. The following would take a
        very long time without the ``maximize_at_primes`` option::

            sage: K.<a> = NumberField(x^2 - p*q, maximize_at_primes=[p])
            sage: K.pari_nf()
            [y^2 - 100000000000000000000...]

        Since the discriminant is square-free, this also works::

            sage: K.<a> = NumberField(x^2 - p*q, assume_disc_small=True)
            sage: K.pari_nf()
            [y^2 - 100000000000000000000...]
        """
        if self.absolute_polynomial().denominator() != 1:
            raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
        try:
            return self.__pari_nf
        except AttributeError:
            f = self.pari_polynomial("y")
            self.__pari_nf = pari([f, self._pari_integral_basis(important=important)]).nfinit()
            return self.__pari_nf

    def pari_zk(self):
        """
        Integral basis of the PARI number field corresponding to this field.

        This is the same as pari_nf().getattr('zk'), but much faster.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 17)
            sage: k.pari_zk()
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
            sage: k.pari_nf().getattr('zk')
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
        """
        return self.pari_nf().nf_get_zk()

    def _pari_(self):
        """
        Converts this number field to PARI.

        This only works if the defining polynomial of this number field is
        integral and monic.

        EXAMPLES::

            sage: k = NumberField(x^2 + x + 1, 'a')
            sage: k._pari_()
            [y^2 + y + 1, [0, 1], -3, 1, ... [1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
            sage: pari(k)
            [y^2 + y + 1, [0, 1], -3, 1, ...[1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
        """
        return self.pari_nf()

    def _pari_init_(self):
        """
        Converts this number field to PARI.

        This only works if the defining polynomial of this number field is
        integral and monic.

        EXAMPLES::

            sage: k = NumberField(x^2 + x + 1, 'a')
            sage: k._pari_init_()
            '[y^2 + y + 1, [0, 1], -3, 1, ... [1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]'
            sage: gp(k)
            [y^2 + y + 1, [0, 1], -3, 1, ...[1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
        """
        return str(self.pari_nf())

    def pari_bnf(self, proof=None, units=True):
        """
        PARI big number field corresponding to this field.

        INPUT:

        - ``proof`` -- If False, assume GRH.  If True, run PARI's
          ``bnfcertify()`` to make sure that the results are correct.

        - ``units`` -- (default: True) If True, insist on having
          fundamental units.  If False, the units may or may not be
          computed.

        OUTPUT:

        The PARI ``bnf`` structure of this number field.

        .. warning::

           Even with ``proof=True``, I wouldn't trust this to mean
           that everything computed involving this number field is
           actually correct.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: len(k.pari_bnf())
            10
            sage: k.pari_bnf()[:4]
            [[;], matrix(0,7), [;], ...]
            sage: len(k.pari_nf())
            9
            sage: k.<a> = NumberField(x^7 + 7); k
            Number Field in a with defining polynomial x^7 + 7
            sage: dummy = k.pari_bnf(proof=True)
        """
        proof = get_flag(proof, "number_field")
        # First compute bnf
        try:
            bnf = self.__pari_bnf
        except AttributeError:
            if self.absolute_polynomial().denominator() != 1:
                raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
            f = self.pari_polynomial("y")
            if units:
                self.__pari_bnf = f.bnfinit(1)
            else:
                self.__pari_bnf = f.bnfinit()
            bnf = self.__pari_bnf
        # Certify if needed
        if proof and not getattr(self, "__pari_bnf_certified", False):
            if bnf.bnfcertify() != 1:
                raise ValueError, "The result is not correct according to bnfcertify"
            self.__pari_bnf_certified = True
        return bnf

    def pari_rnfnorm_data(self, L, proof=True):
        """
        Return the PARI rnfisnorminit() data corresponding to the
        extension L/self.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K = NumberField(x^2 - 2, 'alpha')
            sage: L = K.extension(x^2 + 5, 'gamma')
            sage: ls = K.pari_rnfnorm_data(L) ; len(ls)
            8

            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: P.<X> = K[]
            sage: L.<b> = NumberField(X^3 + a)
            sage: ls = K.pari_rnfnorm_data(L); len(ls)
            8
        """
        if L.base_field() != self:
            raise ValueError, "L must be an extension of self"

        relpoly = L.defining_polynomial()
        Kbnf = self.pari_bnf(proof=proof)
        return Kbnf.rnfisnorminit(relpoly._pari_with_name())

    def _gap_init_(self):
        """
        Create a gap object representing self and return its name

        EXAMPLE::

            sage: z = QQ['z'].0
            sage: K.<zeta> = NumberField(z^2 - 2)
            sage: K._gap_init_() # the following variable name $sage1 represents the F.base_ring() in gap and is somehow random
            'CallFuncList(function() local z,E; z:=Indeterminate($sage1,"z"); E:=AlgebraicExtension($sage1,z^2 - 2,"zeta"); return E; end,[])'
            sage: k = gap(K)
            sage: k
            <algebraic extension over the Rationals of degree 2>
            sage: k.GeneratorsOfDivisionRing()
            [ zeta ]

        The following tests that it is possible to use a defining
        polynomial in the variable ``E``, even though by default
        ``E`` is used as a local variable in the above GAP
        ``CallFuncList``::

            sage: P.<E> = QQ[]
            sage: L.<tau> = NumberField(E^3 - 2)
            sage: l = gap(L); l
            <algebraic extension over the Rationals of degree 3>
            sage: l.GeneratorsOfField()
            [ tau ]
            sage: gap(tau)^3
            !2

        """
        if not self.is_absolute():
            raise NotImplementedError, "Currently, only simple algebraic extensions are implemented in gap"
        G = sage.interfaces.gap.gap
        q = self.polynomial()
        if q.variable_name()!='E':
            return 'CallFuncList(function() local %s,E; %s:=Indeterminate(%s,"%s"); E:=AlgebraicExtension(%s,%s,"%s"); return E; end,[])'%(q.variable_name(),q.variable_name(),G(self.base_ring()).name(),q.variable_name(),G(self.base_ring()).name(),self.polynomial().__repr__(),str(self.gen()))
        else:
            return 'CallFuncList(function() local %s,F; %s:=Indeterminate(%s,"%s"); F:=AlgebraicExtension(%s,%s,"%s"); return F; end,[])'%(q.variable_name(),q.variable_name(),G(self.base_ring()).name(),q.variable_name(),G(self.base_ring()).name(),self.polynomial().__repr__(),str(self.gen()))

    def characteristic(self):
        """
        Return the characteristic of this number field, which is of course
        0.

        EXAMPLES::

            sage: k.<a> = NumberField(x^99 + 2); k
            Number Field in a with defining polynomial x^99 + 2
            sage: k.characteristic()
            0
        """
        return ZZ.zero()

    def class_group(self, proof=None, names='c'):
        r"""
        Return the class group of the ring of integers of this number
        field.

        INPUT:


        -  ``proof`` - if True then compute the class group
           provably correctly. Default is True. Call number_field_proof to
           change this default globally.

        -  ``names`` - names of the generators of this class
           group.


        OUTPUT: The class group of this number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: G = K.class_group(); G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
            sage: G.0
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: G.gens()
            (Fractional ideal class (2, 1/2*a - 1/2),)

        ::

            sage: G.number_field()
            Number Field in a with defining polynomial x^2 + 23
            sage: G is K.class_group()
            True
            sage: G is K.class_group(proof=False)
            False
            sage: G.gens()
            (Fractional ideal class (2, 1/2*a - 1/2),)

        There can be multiple generators::

            sage: k.<a> = NumberField(x^2 + 20072)
            sage: G = k.class_group(); G
            Class group of order 76 with structure C38 x C2 of Number Field in a with defining polynomial x^2 + 20072
            sage: G.0 # random
            Fractional ideal class (41, a + 10)
            sage: G.0^38
            Trivial principal fractional ideal class
            sage: G.1 # random
            Fractional ideal class (2, -1/2*a)
            sage: G.1^2
            Trivial principal fractional ideal class

        Class groups of Hecke polynomials tend to be very small::

            sage: f = ModularForms(97, 2).T(2).charpoly()
            sage: f.factor()
            (x - 3) * (x^3 + 4*x^2 + 3*x - 1) * (x^4 - 3*x^3 - x^2 + 6*x - 1)
            sage: [NumberField(g,'a').class_group().order() for g,_ in f.factor()]
            [1, 1, 1]
        """
        proof = proof_flag(proof)
        try:
            return self.__class_group[proof, names]
        except KeyError:
            pass
        except AttributeError:
            self.__class_group = {}
        k = self.pari_bnf(proof)
        cycle_structure = tuple( ZZ(c) for c in k.bnf_get_cyc() )

        # Gens is a list of ideals (the generators)
        gens = tuple( self.ideal(hnf) for hnf in k.bnf_get_gen() )

        G = ClassGroup(cycle_structure, names, self, gens, proof=proof)
        self.__class_group[proof, names] = G
        return G

    def class_number(self, proof=None):
        """
        Return the class number of this number field, as an integer.

        INPUT:


        -  ``proof`` - bool (default: True unless you called
           number_field_proof)


        EXAMPLES::

            sage: NumberField(x^2 + 23, 'a').class_number()
            3
            sage: NumberField(x^2 + 163, 'a').class_number()
            1
            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').class_number(proof=False)
            1539
        """
        proof = proof_flag(proof)
        return self.class_group(proof).order()

    def S_class_group(self, S, proof=None, names='c'):
        """
        Returns the S-class group of this number field over its base field.

        INPUT:

        - ``S`` - a set of primes of the base field

        - ``proof`` - if False, assume the GRH in computing the class group.
          Default is True. Call ``number_field_proof`` to change this
          default globally.

        - ``names`` - names of the generators of this class group.

        OUTPUT:

        The S-class group of this number field.

        EXAMPLE:

        A well known example::

            sage: K.<a> = QuadraticField(-5)
            sage: K.S_class_group([])
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 5

        When we include the prime `(2, a+1)`, the S-class group becomes
        trivial::

            sage: K.S_class_group([K.ideal(2,a+1)])
            S-class group of order 1 of Number Field in a with defining polynomial x^2 + 5

        TESTS::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S);CS
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14
            sage: T = tuple([])
            sage: CT = K.S_class_group(T);CT
            S-class group of order 4 with structure C4 of Number Field in a with defining polynomial x^2 + 14
            sage: K.class_group()
            Class group of order 4 with structure C4 of Number Field in a with defining polynomial x^2 + 14
        """
        proof = proof_flag(proof)
        if all(P.is_principal() for P in S):
            C = self.class_group(proof=proof)
            Slist = zip([g.ideal() for g in C.gens()], C.invariants())
        else:
            Slist = self._S_class_group_and_units(tuple(S), proof=proof)[1]
        return SClassGroup(tuple(s[1] for s in Slist), names, self,
                           tuple(s[0] for s in Slist), tuple(S))

    def S_units(self, S, proof=True):
        """
        Returns a list of generators of the S-units.

        INPUT::

        - ``S`` - a set of primes of the base field

        - ``proof`` - if False, assume the GRH in computing the class group

        OUTPUT:

        A list of generators of the unit group.

       .. note::

            For more functionality see the S_unit_group() function.

        EXAMPLE::

            sage: K.<a> = QuadraticField(-3)
            sage: K.unit_group()
            Unit group with structure C6 of Number Field in a with defining polynomial x^2 + 3
            sage: K.S_units([])
            [-1/2*a + 1/2]
            sage: K.S_units([])[0].multiplicative_order()
            6

        An example in a relative extension (see :trac:`8722`)::

            sage: L.<a,b> = NumberField([x^2 + 1, x^2 - 5])
            sage: p = L.ideal((-1/2*b - 1/2)*a + 1/2*b - 1/2)
            sage: W = L.S_units([p]); [x.norm() for x in W]
            [9, 1, 1]

        Our generators should have the correct parent (:trac:`9367`)::

            sage: _.<x> = QQ[]
            sage: L.<alpha> = NumberField(x^3 + x + 1)
            sage: p = L.S_units([ L.ideal(7) ])
            sage: p[0].parent()
            Number Field in alpha with defining polynomial x^3 + x + 1

        """
        return self._S_class_group_and_units(tuple(S), proof=proof)[0]

    @cached_method
    def _S_class_group_and_units(self, S, proof=True):
        """
        Compute S class group and units.

        INPUT:

        - ``S`` - a tuple of prime ideals of self

        - ``proof`` - if False, assume the GRH in computing the class group

        OUTPUT:

        - ``units, clgp_gens``, where:

        - ``units`` - A list of generators of the unit group.

        - ``clgp_gens`` - A list of generators of the `S`-class group.
          Each generator is represented as a pair ``(gen, order)``,
          where ``gen`` is a fractional ideal of self and ``order`` is
          its order in the `S`-class group.

        EXAMPLE::

            sage: K.<a> = NumberField(x^2+5)
            sage: K._S_class_group_and_units(())
            ([-1], [(Fractional ideal (2, a + 1), 2)])

            sage: K.<a> = NumberField(polygen(QQ))
            sage: K._S_class_group_and_units( (K.ideal(5),) )
            ([5, -1], [])

        TEST::

            sage: K.<a> = NumberField(x^3 - 381 * x + 127)
            sage: K._S_class_group_and_units(tuple(K.primes_above(13)))
            ([-7/13*a^2 - 140/13*a + 36/13,
              14/13*a^2 + 267/13*a - 85/13,
              7/13*a^2 + 127/13*a - 49/13,
              -1,
              1/13*a^2 + 20/13*a - 7/13,
              1/13*a^2 - 19/13*a + 6/13],
             [(Fractional ideal (11, a - 2), 2),
              (Fractional ideal (19, 1/13*a^2 - 45/13*a - 332/13), 2)])
        """
        K_pari = self.pari_bnf(proof=proof)
        S_pari = [p.pari_prime() for p in S]
        result = K_pari.bnfsunit(S_pari)
        units = map(self, result[0]) + self.unit_group().gens_values()
        orders = result[4][1].sage()
        gens = map(self.ideal, result[4][2])
        return units, [(gens[k], orders[k]) for k in range(len(orders)) if orders[k] > 1]

    @cached_method
    def _S_class_group_quotient_matrix(self, S):
        r"""
        Return the matrix of the quotient map from the class group to the
        S-class group. The result is cached.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-21)
            sage: K._S_class_group_quotient_matrix((K.ideal([2, a+1]),))
            [1]
            [0]
            sage: K._S_class_group_quotient_matrix((K.ideal([5, a+2]),))
            [0]
            [1]
            sage: K._S_class_group_quotient_matrix(())
            [1 0]
            [0 1]
            sage: K.<a> = QuadraticField(-105)
            sage: K._S_class_group_quotient_matrix((K.ideal(11, a + 4),))
            [0 0]
            [1 0]
            [0 1]
        """
        from sage.matrix.constructor import matrix
        S_clgp_gens = self._S_class_group_and_units(S)[1]
        a = len(S_clgp_gens)
        c = self.class_group().ngens()
        M = [u[0].ideal_class_log() for u in S_clgp_gens]
        M += [x.ideal_class_log() for x in S]
        M = matrix(ZZ, M)
        A, Q = M.hermite_form(transformation=True)
        assert A[:c] == 1 and A[c:] == 0
        return Q[:c, :a]

    def selmer_group(self, S, m, proof=True):
        r"""
        Compute the Selmer group `K(S,m)`.  This is the subgroup of
        `K^\times/(K^\times)^m` consisting of elements `a` such that
        `K(\sqrt[m]{a})/K` is unramified at all primes of `K` outside
        of `S`.

        INPUT:

        - ``S`` - A set of primes of self.

        - ``m`` - A positive integer.

        - ``proof`` - If False, assume the GRH in computing the class group.

        OUTPUT:

        A list of generators of `K(S,m)`.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-5)
            sage: K.selmer_group((), 2)
            [-1, 2]
            sage: K.selmer_group([K.ideal(2, -a+1)], 2)
            [2, -1]
            sage: K.selmer_group([K.ideal(2, -a+1), K.ideal(3, a+1)], 2)
            [2, a + 1, -1]
            sage: K.selmer_group((K.ideal(2, -a+1),K.ideal(3, a+1)), 4)
            [2, a + 1, -1]
            sage: K.selmer_group([K.ideal(2, -a+1)], 3)
            [2]
            sage: K.selmer_group([K.ideal(2, -a+1), K.ideal(3, a+1)], 3)
            [2, a + 1]
            sage: K.selmer_group([K.ideal(2, -a+1), K.ideal(3, a+1), K.ideal(a)], 3)  # random signs
            [2, a + 1, a]
            sage: K.<a> = NumberField(polygen(QQ))
            sage: K.selmer_group([],5)
            []

        TESTS::

            sage: K.<a> = QuadraticField(-5)
            sage: S = K.selmer_group([K.ideal(2, -a+1), K.ideal(3, a+1), K.ideal(a)], 3)
            sage: S == [2, a + 1, a] or S == [2, a + 1, -a]
            True

        Verify that :trac:`14489` is fixed::

            sage: K.<a> = NumberField(x^3 - 381 * x + 127)
            sage: K.selmer_group(K.primes_above(13), 2)
            [-7/13*a^2 - 140/13*a + 36/13, 14/13*a^2 + 267/13*a - 85/13, 7/13*a^2 + 127/13*a - 49/13, -1, 1/13*a^2 + 20/13*a - 7/13, 1/13*a^2 - 19/13*a + 6/13, -2/13*a^2 - 53/13*a + 92/13, 10/13*a^2 + 44/13*a - 4555/13]
        """
        units, clgp_gens = self._S_class_group_and_units(tuple(S), proof=proof)
        gens = []
        from sage.rings.infinity import Infinity
        for unit in units:
            order = unit.multiplicative_order()
            if order == Infinity or order.gcd(m) != 1:
                gens.append(unit)
        card_S = len(S)
        if card_S != 0:
            from sage.matrix.constructor import Matrix
            H = self.class_group()
            ords = [g.order() for g in H.gens()]
            pari_ords = pari(ords).Col()
            Sords = [H(s).order() for s in S]
            MS = Matrix(ZZ, [H(s).exponents() for s in S]).transpose()
            pari_MS = pari(MS)
        for gen, order in clgp_gens:
            if order.gcd(m) != 1:
                # The ideal I = gen^(order/gcd) is of order gcd in
                # Cl_S[m].  After multiplying by primes in S, the
                # ideal I^gcd = gen^order becomes principal.  We take
                # a generator of this ideal to get the corresponding
                # generator of the m-Selmer group.
                J = gen ** order
                if card_S != 0 and not J.is_principal():
                    B = H(J).exponents()
                    pari_B = (-pari(B)).Col()
                    exps = pari_MS.matsolvemod(pari_ords, pari_B).Vec().sage()
                    Spart = prod([S[i] ** (exps[i] % Sords[i]) for i in range(card_S)])
                    J *= Spart
                gens.append(self(J.gens_reduced()[0]))
        return gens

    def composite_fields(self, other, names=None, both_maps=False, preserve_embedding=True):
        """
        List of all possible composite number fields formed from self and
        other, together with (optionally) embeddings into the compositum;
        see the documentation for both_maps below.

        If preserve_embedding is True and if self and other both have
        embeddings into the same ambient field, or into fields which are
        contained in a common field, only the compositum respecting
        both embeddings is returned.  If one (or both) of self or other
        does not have an embedding or preserve_embedding is False,
        all possible composite number fields are returned.

        INPUT:

        - ``other`` - a number field

        - ``names`` - generator name for composite fields

        - ``both_maps`` - (default: False)  if True, return quadruples
          (F, self_into_F, other_into_F, k) such that self_into_F is an
          embedding of self in F, other_into_F is an embedding of in F,
          and k is an integer such that F.gen() equals
          other_into_F(other.gen()) + k*self_into_F(self.gen())
          or has the value Infinity in which case F.gen() equals
          self_into_F(self.gen()), or is None (which happens when other is a
          relative number field).
          If both self and other have embeddings into an ambient field, then
          F will have an embedding with respect to which both self_into_F
          and other_into_F will be compatible with the ambient embeddings.

        - ``preserve_embedding`` - (default: True) if self and other have
          ambient embeddings, then return only the compatible compositum.

        OUTPUT:

        -  ``list`` - list of the composite fields, possibly with maps.


        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 2)
            sage: K.composite_fields(K)
            [Number Field in a with defining polynomial x^4 - 2,
             Number Field in a0 with defining polynomial x^8 + 28*x^4 + 2500]

        A particular compositum is selected, together with compatible maps
        into the compositum, if the fields are endowed with a real or
        complex embedding::

            sage: K1 = NumberField(x^4 - 2, 'a', embedding=RR(2^(1/4)))
            sage: K2 = NumberField(x^4 - 2, 'a', embedding=RR(-2^(1/4)))
            sage: K1.composite_fields(K2)
            [Number Field in a with defining polynomial x^4 - 2]
            sage: [F, f, g, k], = K1.composite_fields(K2, both_maps=True); F
            Number Field in a with defining polynomial x^4 - 2
            sage: f(K1.0), g(K2.0)
            (a, -a)

        With preserve_embedding set to False, the embeddings are ignored::

            sage: K1.composite_fields(K2, preserve_embedding=False)
            [Number Field in a with defining polynomial x^4 - 2,
             Number Field in a0 with defining polynomial x^8 + 28*x^4 + 2500]

        Changing the embedding selects a different compositum::

            sage: K3 = NumberField(x^4 - 2, 'a', embedding=CC(2^(1/4)*I))
            sage: [F, f, g, k], = K1.composite_fields(K3, both_maps=True); F
            Number Field in a0 with defining polynomial x^8 + 28*x^4 + 2500
            sage: f(K1.0), g(K3.0)
            (1/240*a0^5 - 41/120*a0, 1/120*a0^5 + 19/60*a0)

        If no embeddings are specified, the maps into the composite are chosen arbitrarily::

            sage: Q1.<a> = NumberField(x^4 + 10*x^2 + 1)
            sage: Q2.<b> = NumberField(x^4 + 16*x^2 + 4)
            sage: Q1.composite_fields(Q2, 'c')
            [Number Field in c with defining polynomial x^8 + 64*x^6 + 904*x^4 + 3840*x^2 + 3600]
            sage: F, Q1_into_F, Q2_into_F, k = Q1.composite_fields(Q2, 'c', both_maps=True)[0]
            sage: Q1_into_F
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 + 10*x^2 + 1
              To:   Number Field in c with defining polynomial x^8 + 64*x^6 + 904*x^4 + 3840*x^2 + 3600
              Defn: a |--> 19/14400*c^7 + 137/1800*c^5 + 2599/3600*c^3 + 8/15*c

        This is just one of four embeddings of Q1 into F::
            sage: Hom(Q1, F).order()
            4

        TESTS:

        Let's check that embeddings are being respected::

            sage: x = polygen(ZZ)
            sage: K0.<b> = CyclotomicField(7, 'a').subfields(3)[0][0].change_names()
            sage: K1.<a1> = K0.extension(x^2 - 2*b^2, 'a1').absolute_field()
            sage: K2.<a2> = K0.extension(x^2 - 3*b^2, 'a2').absolute_field()

        We need embeddings, so we redefine::

            sage: L1.<a1> = NumberField(K1.polynomial(), 'a1', embedding=CC.0)
            sage: L2.<a2> = NumberField(K2.polynomial(), 'a2', embedding=CC.0)
            sage: [CDF(a1), CDF(a2)]
            [-0.629384245426, -0.77083512672]

        and we get the same embeddings via the compositum::

            sage: F, L1_into_F, L2_into_F, k = L1.composite_fields(L2, both_maps=True)[0]
            sage: [CDF(L1_into_F(L1.gen())), CDF(L2_into_F(L2.gen()))]
            [-0.629384245426, -0.77083512672]

        Let's check that if only one field has an embedding, the resulting
        fields do not have embeddings::

            sage: L1.composite_fields(K2)[0].coerce_embedding() is None
            True
            sage: L2.composite_fields(K1)[0].coerce_embedding() is None
            True

        We check that other can be a relative number field::

            sage: L.<a, b> = NumberField([x^3 - 5, x^2 + 3])
            sage: CyclotomicField(3, 'w').composite_fields(L, both_maps=True)
            [(Number Field in a with defining polynomial x^3 - 5 over its base field, Ring morphism:
              From: Cyclotomic Field of order 3 and degree 2
              To:   Number Field in a with defining polynomial x^3 - 5 over its base field
              Defn: w |--> -1/2*b - 1/2, Relative number field endomorphism of Number Field in a with defining polynomial x^3 - 5 over its base field
              Defn: a |--> a
                    b |--> b, None)]
        """
        if not isinstance(other, NumberField_generic):
            raise TypeError, "other must be a number field."

        sv = self.variable_name(); ov = other.variable_name()
        if names is None:
            names = sv + (ov if ov != sv else "")
        name = sage.structure.parent_gens.normalize_names(1, names)[0]

        # should we try to preserve embeddings?
        subfields_have_embeddings = preserve_embedding
        if self.coerce_embedding() is None:
            subfields_have_embeddings = False
        if other.coerce_embedding() is None:
            subfields_have_embeddings = False
        if subfields_have_embeddings:
            try:
                from sage.categories.pushout import pushout
                ambient_field = pushout(self.coerce_embedding().codomain(), other.coerce_embedding().codomain())
            except CoercionException:
                ambient_field = None
            if ambient_field is None:
                subfields_have_embeddings = False

        f = self.pari_polynomial()
        g = other.pari_polynomial()
        R = self.absolute_polynomial().parent()

        m = self.degree()
        n = other.absolute_degree()

        if not both_maps and not subfields_have_embeddings:
            # short cut!
            # eliminate duplicates from the fields given by polcompositum
            # and return the resulting number fields.  There is no need to
            # check that the polynomials are irreducible.
            C = []
            for r in f.polcompositum(g):
                if not any(r.nfisisom(s) for s in C):
                    C.append(r)
            C = map(R, C)

            q = sum(1 for r in C if r.degree() != max(m, n))
            if q == 1 and name != sv and name != ov:
                names =[name]
            else:
                names = [name + str(i) for i in range(q)]

            i = 0
            rets = []
            for r in C:
                if r.degree() == m:
                    rets.append(self)
                elif r.degree() == n:
                    rets.append(other)
                else:
                    rets.append(NumberField(r, names[i], check=False))
                    i += 1
            return rets

        # If flag = 1, polcompositum outputs a vector of 4-component vectors
        # [R, a, b, k], where R ranges through the list of all possible compositums
        # as above, and a (resp. b) expresses the root of P (resp. Q) as
        # an element of Q(X )/(R). Finally, k is a small integer such that
        # b + ka = X modulo R.
        # In this case duplicates must only be eliminated if embeddings are going
        # to be preserved.
        C = []
        for v in f.polcompositum(g, 1):
            if subfields_have_embeddings or not any(v[0].nfisisom(u[0]) for u in C):
                C.append(v)

        a = self.gen()
        b = other.gen()

        # If both subfields are provided with embeddings, then we must select
        # the compositum which corresponds to these embeddings.  We do this by
        # evaluating the given polynomials at the corresponding embedded values.
        # For the case we want, the result will be zero, but rounding errors are
        # difficult to predict, so we just take the field which yields the
        # mimumum value.
        if subfields_have_embeddings:
            poly_vals = []
            for r, _, _, k in C:
                r = R(r)
                k = ZZ(k) # essential
                embedding = other.coerce_embedding()(b) + k*self.coerce_embedding()(a)
                poly_vals.append(sage.rings.complex_double.CDF(r(embedding)).abs())
            i = poly_vals.index(min(poly_vals))
            C = [C[i]]

        q = sum(1 for r, _, _, _ in C if R(r).degree() != max(m, n))
        if q == 1 and name != sv and name != ov:
            names =[name, '']
        else:
            names = [name + str(i) for i in range(q + 1)]

        if both_maps and not other.is_absolute():
            other_abs = other.absolute_field('z')
            from_other_abs, to_other_abs = other_abs.structure()

        embedding = None
        i = 0
        rets = []
        for r, a_in_F, b_in_F, k in C:
            r = R(r)
            if r.degree() == m and not both_maps:
                rets.append(self)
            elif r.degree() == n and not both_maps:
                rets.append(other)
            else:
                k = ZZ(k) # essential

                if subfields_have_embeddings:
                    embedding = other.coerce_embedding()(b) + k*self.coerce_embedding()(a)
                F = NumberField(r, names[i], check=False, embedding=embedding)
                i += 1
                if both_maps:
                    if other.is_absolute():
                        if r.degree() == m:
                            self_to_F = self.hom([a])
                            other_to_F = other.hom([(~self.hom([F(a_in_F)]))(F(b_in_F))])
                            F = self
                            k = sage.rings.infinity.Infinity
                            i -= 1
                        elif r.degree() == n:
                            other_to_F = other.hom([b])
                            self_to_F = self.hom([(~other.hom([F(b_in_F)]))(F(a_in_F))])
                            F = other
                            k = ZZ(0)
                            i -= 1
                        else:
                            self_to_F = self.hom([F(a_in_F)])
                            other_to_F = other.hom([F(b_in_F)])
                    else:
                        other_abs_to_F = other_abs.hom([F(b_in_F)])
                        other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(F), other_abs_to_F*to_other_abs)
                        if r.degree() == m:
                            self_to_F = self.hom([a])
                            other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(self), (~self.hom([F(a_in_F)]))*other_abs_to_F*to_other_abs)
                            F = self
                            k = None
                            i -= 1
                        elif r.degree() == n:
                            other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(other), from_other_abs)
                            self_to_F = self.hom([from_other_abs((~other_abs_to_F)(F(a_in_F)))])
                            F = other
                            k = None
                            i -= 1
                        else:
                            self_to_F = self.hom([F(a_in_F)])
                            other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(F), other_abs_to_F*to_other_abs)
                    rets.append( (F, self_to_F, other_to_F, k) )
                else:
                    rets.append(F)
        return rets

    def absolute_degree(self):
        """
        Return the degree of self over `\QQ`.

        EXAMPLES::

            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').absolute_degree()
            3
            sage: NumberField(x + 1, 'a').absolute_degree()
            1
            sage: NumberField(x^997 + 17*x + 3, 'a', check=False).absolute_degree()
            997
        """
        return self.polynomial().degree()

    def degree(self):
        """
        Return the degree of this number field.

        EXAMPLES::

            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').degree()
            3
            sage: NumberField(x + 1, 'a').degree()
            1
            sage: NumberField(x^997 + 17*x + 3, 'a', check=False).degree()
            997
        """
        return self.polynomial().degree()

    def different(self):
        r"""
        Compute the different fractional ideal of this number field.

        The codifferent is the fractional ideal of all `x` in `K`
        such that the the trace of `xy` is an integer for
        all `y \in O_K`.

        The different is the integral ideal which is the inverse of
        the codifferent.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 23)
            sage: d = k.different()
            sage: d
            Fractional ideal (-a)
            sage: d.norm()
            23
            sage: k.disc()
            -23

        The different is cached::

            sage: d is k.different()
            True

        Another example::

            sage: k.<b> = NumberField(x^2 - 123)
            sage: d = k.different(); d
            Fractional ideal (2*b)
            sage: d.norm()
            492
            sage: k.disc()
            492
        """
        try:
            return self.__different
        except AttributeError:
            self.__different = self.ideal(self.pari_nf().nf_get_diff())
            return self.__different

    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the number
        field, or if v is specified, the determinant of the trace pairing
        on the elements of the list v.

        INPUT:


        -  ``v (optional)`` - list of element of this number
           field


        OUTPUT: Integer if v is omitted, and Rational otherwise.

        EXAMPLES::

            sage: K.<t> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K.disc()
            -503
            sage: K.disc([1, t, t^2])
            -2012
            sage: K.disc([1/7, (1/5)*t, (1/3)*t^2])
            -2012/11025
            sage: (5*7*3)^2
            11025
            sage: NumberField(x^2 - 1/2, 'a').discriminant()
            8
        """
        if v == None:
            try:
                return self.__disc
            except AttributeError:
                self.__disc = ZZ(self.pari_polynomial().nfdisc())
                return self.__disc
        else:
            return QQ(self.trace_pairing(v).det())

    def disc(self, v=None):
        """
        Shortcut for self.discriminant.

        EXAMPLES::

            sage: k.<b> = NumberField(x^2 - 123)
            sage: k.disc()
            492
        """
        return self.discriminant(v=v)

    def trace_dual_basis(self, b):
        r"""
        Compute the dual basis of a basis of ``self`` with respect to the trace pairing.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: b = [1, 2*a, 3*a^2]
            sage: T = K.trace_dual_basis(b); T
            [4/31*a^2 - 6/31*a + 13/31, -9/62*a^2 - 1/31*a - 3/31, 2/31*a^2 - 3/31*a + 4/93]
            sage: [(b[i]*T[j]).trace() for i in xrange(3) for j in xrange(3)]
            [1, 0, 0, 0, 1, 0, 0, 0, 1]
        """
        if not len(b) == self.degree():
            raise ValueError, 'Not a basis of the number field.'
        M = self.trace_pairing(b)
        if not M.is_invertible():
            raise ValueError, 'Not a basis of the number field.'
        return [sum([v[i]*b[i] for i in xrange(len(b))]) for v in M.inverse()]

    def elements_of_norm(self, n, proof=None):
        r"""
        Return a list of solutions modulo units of positive norm to
        `Norm(a) = n`, where a can be any integer in this number
        field.

        INPUT:


        -  ``proof`` - default: True, unless you called
           number_field_proof and set it otherwise.


        EXAMPLES::

            sage: K.<a> = NumberField(x^2+1)
            sage: K.elements_of_norm(3)
            []
            sage: K.elements_of_norm(50)
            [-7*a + 1, -5*a - 5, a - 7]
        """
        proof = proof_flag(proof)
        B = self.pari_bnf(proof).bnfisintnorm(n)
        return map(self, B)

    def extension(self, poly, name=None, names=None, check=True, embedding=None):
        """
        Return the relative extension of this field by a given polynomial.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: R.<t> = K[]
            sage: L.<b> = K.extension(t^2 + a); L
            Number Field in b with defining polynomial t^2 + a over its base field

        We create another extension::

            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: y = var('y')
            sage: m.<b> = k.extension(y^2 + 2); m
            Number Field in b with defining polynomial y^2 + 2 over its base field

        Note that b is a root of `y^2 + 2`::

            sage: b.minpoly()
            x^2 + 2
            sage: b.minpoly('z')
            z^2 + 2

        A relative extension of a relative extension::

            sage: k.<a> = NumberField([x^2 + 1, x^3 + x + 1])
            sage: R.<z> = k[]
            sage: L.<b> = NumberField(z^3 + 3 + a); L
            Number Field in b with defining polynomial z^3 + a0 + 3 over its base field
        """
        if not isinstance(poly, polynomial_element.Polynomial):
            try:
                poly = poly.polynomial(self)
            except (AttributeError, TypeError):
                raise TypeError, "polynomial (=%s) must be a polynomial."%repr(poly)
        if not names is None:
            name = names
        if isinstance(name, tuple):
            name = name[0]
        if name is None:
            raise TypeError, "the variable name must be specified."
        from sage.rings.number_field.number_field_rel import NumberField_relative
        return NumberField_relative(self, poly, str(name), check=check, embedding=embedding)

    def factor(self, n):
        r"""
        Ideal factorization of the principal ideal generated by `n`.

        EXAMPLES:

        Here we show how to factor Gaussian integers (up to units).
        First we form a number field defined by `x^2 + 1`::

            sage: K.<I> = NumberField(x^2 + 1); K
            Number Field in I with defining polynomial x^2 + 1

        Here are the factors::

            sage: fi, fj = K.factor(17); fi,fj
            ((Fractional ideal (4*I + 1), 1), (Fractional ideal (-I - 4), 1))

        Now we extract the reduced form of the generators::

            sage: zi = fi[0].gens_reduced()[0]; zi
            4*I + 1
            sage: zj = fj[0].gens_reduced()[0]; zj
            -I - 4

        We recover the integer that was factored in `\ZZ[i]` (up to a unit)::

            sage: zi*zj
            -17*I

        One can also factor elements or ideals of the number field::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.factor(1/3)
            (Fractional ideal (3))^-1
            sage: K.factor(1+a)
            Fractional ideal (a + 1)
            sage: K.factor(1+a/5)
            (Fractional ideal (-3*a - 2)) * (Fractional ideal (a + 1)) * (Fractional ideal (-a - 2))^-1 * (Fractional ideal (a - 2))^-1

        An example over a relative number field::

            sage: pari('setrand(2)')
            sage: L.<b> = K.extension(x^2 - 7)
            sage: f = L.factor(a + 1); f
            (Fractional ideal (-1/2*a*b - a - 1/2)) * (Fractional ideal (1/2*b + 1/2*a - 1)) # 32-bit
            (Fractional ideal (1/2*a*b + a + 1/2)) * (Fractional ideal (-1/2*b - 1/2*a + 1)) # 64-bit
            sage: f.value() == a+1
            True

        It doesn't make sense to factor the ideal (0), so this raises an error::

            sage: L.factor(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'factor'

        AUTHORS:

        - Alex Clemesha (2006-05-20), Francis Clarke (2009-04-21): examples
        """
        return self.ideal(n).factor()

    def prime_factors(self, x):
        """
        Return a list of the prime ideals of self which divide
        the ideal generated by `x`.

        OUTPUT: list of prime ideals (a new list is returned each time this
        function is called)

        EXAMPLES::

            sage: K.<w> = NumberField(x^2 + 23)
            sage: K.prime_factors(w + 1)
            [Fractional ideal (2, 1/2*w - 1/2), Fractional ideal (2, 1/2*w + 1/2), Fractional ideal (3, 1/2*w + 1/2)]
        """
        return self.ideal(x).prime_factors()

    def gen(self, n=0):
        """
        Return the generator for this number field.

        INPUT:


        -  ``n`` - must be 0 (the default), or an exception is
           raised.


        EXAMPLES::

            sage: k.<theta> = NumberField(x^14 + 2); k
            Number Field in theta with defining polynomial x^14 + 2
            sage: k.gen()
            theta
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0:
            raise IndexError, "Only one generator."
        try:
            return self.__gen
        except AttributeError:
            if self.__polynomial != None:
                X = self.__polynomial.parent().gen()
            else:
                X = PolynomialRing(rational_field.RationalField()).gen()
            self.__gen = self._element_class(self, X)
            return self.__gen

    def is_field(self, proof=True):
        """
        Return True since a number field is a field.

        EXAMPLES::

            sage: NumberField(x^5 + x + 3, 'c').is_field()
            True
        """
        return True

    def is_galois(self):
        r"""
        Return True if this number field is a Galois extension of
        `\QQ`.

        EXAMPLES::

            sage: NumberField(x^2 + 1, 'i').is_galois()
            True
            sage: NumberField(x^3 + 2, 'a').is_galois()
            False
        """
        return self.galois_group(type="pari").order() == self.degree()

    @cached_method
    def galois_group(self, type=None, algorithm='pari', names=None):
        r"""
        Return the Galois group of the Galois closure of this number field.

        INPUT:

        -  ``type`` - ``none``, ``gap``, or ``pari``. If None (the default),
           return an explicit group of automorphisms of self as a
           ``GaloisGroup_v2`` object.  Otherwise, return a ``GaloisGroup_v1``
           wrapper object based on a PARI or Gap transitive group object, which
           is quicker to compute, but rather less useful (in particular, it
           can't be made to act on self).  If type = 'gap', the database_gap
           package should be installed.

        -  ``algorithm`` - 'pari', 'kash', 'magma'. (default: 'pari', except
           when the degree is >= 12 when 'kash' is tried.)

        -  ``name`` - a string giving a name for the generator of the Galois
           closure of self, when self is not Galois. This is ignored if type is
           not None.

        Note that computing Galois groups as abstract groups is often much
        faster than computing them as explicit automorphism groups (but of
        course you get less information out!) For more (important!)
        documentation, so the documentation for Galois groups of polynomials
        over `\QQ`, e.g., by typing ``K.polynomial().galois_group?``,
        where `K` is a number field.

        To obtain actual field homomorphisms from the number field to its
        splitting field, use type=None.

        EXAMPLES:

        With type ``None``::

            sage: k.<b> = NumberField(x^2 - 14) # a Galois extension
            sage: G = k.galois_group(); G
            Galois group of Number Field in b with defining polynomial x^2 - 14
            sage: G.gen(0)
            (1,2)
            sage: G.gen(0)(b)
            -b
            sage: G.artin_symbol(k.primes_above(3)[0])
            (1,2)

            sage: k.<b> = NumberField(x^3 - x + 1) # not Galois
            sage: G = k.galois_group(names='c'); G
            Galois group of Galois closure in c of Number Field in b with defining polynomial x^3 - x + 1
            sage: G.gen(0)
            (1,2)(3,5)(4,6)

        With type ``'pari'``::

            sage: NumberField(x^3-2, 'a').galois_group(type="pari")
            Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in a with defining polynomial x^3 - 2

        ::

            sage: NumberField(x-1, 'a').galois_group(type="gap")    # optional - database_gap
            Galois group Transitive group number 1 of degree 1 of the Number Field in a with defining polynomial x - 1
            sage: NumberField(x^2+2, 'a').galois_group(type="gap")  # optional - database_gap
            Galois group Transitive group number 1 of degree 2 of the Number Field in a with defining polynomial x^2 + 2
            sage: NumberField(x^3-2, 'a').galois_group(type="gap")  # optional - database_gap
            Galois group Transitive group number 2 of degree 3 of the Number Field in a with defining polynomial x^3 - 2

        ::

            sage: x = polygen(QQ)
            sage: NumberField(x^3 + 2*x + 1, 'a').galois_group(pari_group=False)    # optional - database_gap
            Galois group Transitive group number 2 of degree 3 of the Number Field in a with defining polynomial x^3 + 2*x + 1
            sage: NumberField(x^3 + 2*x + 1, 'a').galois_group(algorithm='magma')   # optional - magma, , database_gap
            Galois group Transitive group number 2 of degree 3 of the Number Field in a with defining polynomial x^3 + 2*x + 1

        EXPLICIT GALOIS GROUP: We compute the Galois group as an explicit
        group of automorphisms of the Galois closure of a field.

        ::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<b1> = K.galois_closure(); L
            Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: G = End(L); G
            Automorphism group of Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: G.list()
            [
            Ring endomorphism of Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
              Defn: b1 |--> b1,
            ...
            Ring endomorphism of Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
              Defn: b1 |--> -2/63*b1^4 - 31/63*b1
            ]
            sage: G[1](b1)
            1/36*b1^4 + 1/18*b1
        """
        from galois_group import GaloisGroup_v1, GaloisGroup_v2

        if type is None:
            return GaloisGroup_v2(self, names)

        elif type=="pari":
            return GaloisGroup_v1(self.absolute_polynomial().galois_group(pari_group=True, algorithm=algorithm), self)
        elif type=="gap":
            return GaloisGroup_v1(self.absolute_polynomial().galois_group(pari_group=False, algorithm=algorithm), self)
        else:
            raise ValueError, "Galois group type must be None, 'pari', or 'gap'."

    def _normalize_prime_list(self, v):
        """
        Internal function to convert into a tuple of primes either None or
        a single prime or a list.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K._normalize_prime_list(None)
            ()
            sage: K._normalize_prime_list(3)
            (3,)
            sage: K._normalize_prime_list([3,5])
            (3, 5)
        """
        if v is None:
            v = []
        elif not isinstance(v, (list, tuple)):
            v = [v]
        return tuple(map(ZZ, v))

    def power_basis(self):
        r"""
        Return a power basis for this number field over its base field.

        If this number field is represented as `k[t]/f(t)`, then
        the basis returned is `1, t, t^2, \ldots, t^{d-1}` where
        `d` is the degree of this number field over its base
        field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K.power_basis()
            [1, a, a^2, a^3, a^4]

        ::

            sage: L.<b> = K.extension(x^2 - 2)
            sage: L.power_basis()
            [1, b]
            sage: L.absolute_field('c').power_basis()
            [1, c, c^2, c^3, c^4, c^5, c^6, c^7, c^8, c^9]

        ::

            sage: M = CyclotomicField(15)
            sage: M.power_basis()
            [1, zeta15, zeta15^2, zeta15^3, zeta15^4, zeta15^5, zeta15^6, zeta15^7]
        """
        g = self.gen()
        return [ g**i for i in range(self.relative_degree()) ]

    def integral_basis(self, v=None):
        """
        Returns a list containing a ZZ-basis for the full ring of integers
        of this number field.

        INPUT:


        -  ``v`` - None, a prime, or a list of primes. See the
           documentation for self.maximal_order.


        EXAMPLES::

            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K.integral_basis()
            [1, a, a^2, a^3, a^4]

        Next we compute the ring of integers of a cubic field in which 2 is
        an "essential discriminant divisor", so the ring of integers is not
        generated by a single element.

        ::

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K.integral_basis()
            [1, 1/2*a^2 + 1/2*a, a^2]

        ALGORITHM: Uses the pari library (via _pari_integral_basis).
        """
        return self.maximal_order(v=v).basis()

    def _pari_integral_basis(self, v=None, important=True):
        """
        Internal function returning an integral basis of this number field in
        PARI format.

        INPUT:

        -  ``v`` -- None, a prime, or a list of primes. See the
           documentation for self.maximal_order.

        - ``important`` -- (default:True) bool.  If False, raise a
          ``RuntimeError`` if we need to do a difficult discriminant
          factorization.  Useful when the integral basis is useful
          but not strictly required.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K._pari_integral_basis()
            [1, y, y^2, y^3, y^4]

        Next we compute the ring of integers of a cubic field in which 2 is
        an "essential discriminant divisor", so the ring of integers is not
        generated by a single element.

        ::

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K._pari_integral_basis()
            [1, y, 1/2*y^2 - 1/2*y]
            sage: K.integral_basis()
            [1, 1/2*a^2 + 1/2*a, a^2]
        """
        if (v is None or len(v) == 0) and self._maximize_at_primes:
            v = self._maximize_at_primes

        v = self._normalize_prime_list(v)
        try:
            return self._integral_basis_dict[v]
        except (AttributeError, KeyError):
            f = self.pari_polynomial("y")
            if len(v) > 0:
                m = self._pari_disc_factorization_matrix(v)
                B = f.nfbasis(fa = m)
            elif self._assume_disc_small:
                B = f.nfbasis(1)
            elif not important:
                # Trial divide the discriminant
                m = self.pari_polynomial().poldisc().abs().factor(limit=0)
                # Since we only need a *squarefree* factorization for
                # primes with exponent 1, we need trial division up to D^(1/3)
                # instead of D^(1/2).
                trialdivlimit2 = pari(pari._primelimit()**2)
                trialdivlimit3 = pari(pari._primelimit()**3)
                if all([ p < trialdivlimit2 or (e == 1 and p < trialdivlimit3) or p.isprime() for p,e in zip(m[0],m[1]) ]):
                    B = f.nfbasis(fa = m)
                else:
                    raise RuntimeError, "Unable to factor discriminant with trial division"
            else:
                B = f.nfbasis()

            self._integral_basis_dict[v] = B
            return B

    def _pari_disc_factorization_matrix(self, v):
        """
        Returns a PARI matrix representation for the partial
        factorization of the discriminant of the defining polynomial
        of self, defined by the list of primes in the Python list v.
        This function is used internally by the number fields code.

        EXAMPLES::

            sage: x = polygen(QQ,'x')
            sage: f = x^3 + 17*x + 393; f.discriminant().factor()
            -1 * 5^2 * 29 * 5779
            sage: K.<a> = NumberField(f)
            sage: fa = K._pari_disc_factorization_matrix([5,29]); fa
            [5, 2; 29, 1]
            sage: fa.type()
            't_MAT'
        """
        f = self.pari_polynomial()
        m = pari.matrix(len(v), 2)
        d = f.poldisc()
        for i in range(len(v)):
            p = pari(ZZ(v[i]))
            m[i,0] = p
            m[i,1] = d.valuation(p)
        return m


    def reduced_basis(self, prec=None):
        r"""
        This function returns an LLL-reduced basis for the
        Minkowski-embedding of the maximal order of a number field.

        INPUT:

        -  ``self`` - number field, the base field

        -  ``prec (default: None)`` - the precision with which to
           compute the Minkowski embedding.


        OUTPUT:

        An LLL-reduced basis for the Minkowski-embedding of the
        maximal order of a number field, given by a sequence of (integral)
        elements from the field.

        .. note::

           In the non-totally-real case, the LLL routine we call is
           currently PARI's qflll(), which works with floating point
           approximations, and so the result is only as good as the
           precision promised by PARI. The matrix returned will always
           be integral; however, it may only be only "almost" LLL-reduced
           when the precision is not sufficiently high.

        EXAMPLES::

            sage: F.<t> = NumberField(x^6-7*x^4-x^3+11*x^2+x-1)
            sage: F.maximal_order().basis()
            [1/2*t^5 + 1/2*t^4 + 1/2*t^2 + 1/2, t, t^2, t^3, t^4, t^5]
            sage: F.reduced_basis()
            [-1, -1/2*t^5 + 1/2*t^4 + 3*t^3 - 3/2*t^2 - 4*t - 1/2, t, 1/2*t^5 + 1/2*t^4 - 4*t^3 - 5/2*t^2 + 7*t + 1/2, 1/2*t^5 - 1/2*t^4 - 2*t^3 + 3/2*t^2 - 1/2, 1/2*t^5 - 1/2*t^4 - 3*t^3 + 5/2*t^2 + 4*t - 5/2]
            sage: CyclotomicField(12).reduced_basis()
            [1, zeta12^2, zeta12, zeta12^3]
        """
        if self.is_totally_real():
            try:
                return self.__reduced_basis
            except AttributeError:
                pass
        else:
            try:
                if self.__reduced_basis_precision >= prec:
                    return self.__reduced_basis
            except AttributeError:
                pass

        from sage.matrix.constructor import matrix

        d = self.absolute_degree()
        Z_basis = self.integral_basis()

        ## If self is totally real, then we can use (x*y).trace() as
        ## the inner product on the Minkowski embedding, which is
        ## faster than computing all the conjugates, etc ...
        if self.is_totally_real():
            T = pari(matrix(ZZ, d, d, [[(x*y).trace() for x in Z_basis]
                                       for y in Z_basis])).qflllgram()
            self.__reduced_basis = [ sum([ ZZ(T[i][j]) * Z_basis[j]
                                           for j in range(d)])
                                     for i in range(d)]
        else:
            M = self.Minkowski_embedding(self.integral_basis(), prec=prec)
            T = pari(M).qflll().python()
            self.__reduced_basis = [ self(v.list()) for v in T.columns() ]
            if prec is None:
                ## this is the default choice for Minkowski_embedding
                self.__reduced_basis_prec = 53
            else:
                self.__reduced_basis_prec = prec

        return self.__reduced_basis


    def reduced_gram_matrix(self, prec=None):
        r"""
        This function returns the Gram matrix of an LLL-reduced basis for
        the Minkowski embedding of the maximal order of a number field.

        INPUT:


        -  ``self`` - number field, the base field

        -  ``prec (default: None)`` - the precision with which
           to calculate the Minkowski embedding. (See NOTE below.)


        OUTPUT: The Gram matrix `[\langle x_i,x_j \rangle]` of an LLL reduced
        basis for the maximal order of self, where the integral basis for
        self is given by `\{x_0, \dots, x_{n-1}\}`. Here `\langle , \rangle` is
        the usual inner product on `\RR^n`, and self is embedded in `\RR^n` by
        the Minkowski embedding. See the docstring for
        :meth:`NumberField_absolute.Minkowski_embedding` for more information.

        .. note::

           In the non-totally-real case, the LLL routine we call is
           currently PARI's qflll(), which works with floating point
           approximations, and so the result is only as good as the
           precision promised by PARI. In particular, in this case,
           the returned matrix will *not* be integral, and may not
           have enough precision to recover the correct gram matrix
           (which is known to be integral for theoretical
           reasons). Thus the need for the prec flag above.

        If the following run-time error occurs: "PariError: not a definite
        matrix in lllgram (42)" try increasing the prec parameter,

        EXAMPLES::

            sage: F.<t> = NumberField(x^6-7*x^4-x^3+11*x^2+x-1)
            sage: F.reduced_gram_matrix()
            [ 6  3  0  2  0  1]
            [ 3  9  0  1  0 -2]
            [ 0  0 14  6 -2  3]
            [ 2  1  6 16 -3  3]
            [ 0  0 -2 -3 16  6]
            [ 1 -2  3  3  6 19]
            sage: Matrix(6, [(x*y).trace() for x in F.integral_basis() for y in F.integral_basis()])
            [2550  133  259  664 1368 3421]
            [ 133   14    3   54   30  233]
            [ 259    3   54   30  233  217]
            [ 664   54   30  233  217 1078]
            [1368   30  233  217 1078 1371]
            [3421  233  217 1078 1371 5224]

        ::

            sage: var('x')
            x
            sage: F.<alpha> = NumberField(x^4+x^2+712312*x+131001238)
            sage: F.reduced_gram_matrix(prec=128)
            [   4.0000000000000000000000000000000000000   0.00000000000000000000000000000000000000 -2.1369320000000000000000000000000000000e6 -3.3122478000000000000000000000000000000e7]
            [  0.00000000000000000000000000000000000000    46721.539331563218381658483353092335550 -2.2467769057394530109094755223395819322e7 -3.4807276041138450473611629088647496430e8]
            [-2.1369320000000000000000000000000000000e6 -2.2467769057394530109094755223395819322e7 7.0704243186034907491782135494859351061e12 1.1256636615786237006027526953641297995e14]
            [-3.3122478000000000000000000000000000000e7 -3.4807276041138450473611629088647496430e8 1.1256636615786237006027526953641297995e14 1.7923838231014970520503146603069479547e15]
        """
        if self.is_totally_real():
            try:
                return self.__reduced_gram_matrix
            except AttributeError:
                pass
        else:
            try:
                if self.__reduced_gram_matrix_prec >= prec:
                    return self.__reduced_gram_matrix
            except AttributeError:
                pass

        from sage.matrix.constructor import matrix
        from sage.misc.flatten import flatten
        d = self.absolute_degree()

        if self.is_totally_real():
            B = self.reduced_basis()
            self.__reduced_gram_matrix = matrix(ZZ, d, d,
                                                [[(x*y).trace() for x in B]
                                                 for y in B])
        else:
            M = self.Minkowski_embedding(prec=prec)
            T = matrix(d, flatten([ a.vector().list()
                                    for a in self.reduced_basis(prec=prec) ]))
            A = M*(T.transpose())
            self.__reduced_gram_matrix = A.transpose()*A
            if prec is None:
                ## this is the default choice for Minkowski_embedding
                self.__reduced_gram_matrix_prec = 53
            else:
                self.__reduced_gram_matrix_prec = prec

        return self.__reduced_gram_matrix


    #******************************************************
    # Supplementary algorithm to enumerate lattice points
    #******************************************************

    def _positive_integral_elements_with_trace(self, C):
        r"""
        Find all totally positive integral elements in self whose
        trace is between C[0] and C[1], inclusive.

        .. note::

           This is currently only implemented in the case that self is
           totally real, since it requires exact computation of
           :meth:`.reduced_gram_matrix`.

        EXAMPLES::

            sage: K.<alpha> = NumberField(ZZ['x'].0^2-2)
            sage: K._positive_integral_elements_with_trace([0,5])
            [alpha + 2, -alpha + 2, 2, 1]
            sage: L.<beta> = NumberField(ZZ['x'].0^2+1)
            sage: L._positive_integral_elements_with_trace([5,11])
            Traceback (most recent call last):
            ...
            NotImplementedError: exact computation of LLL reduction only implemented in the totally real case
            sage: L._positive_integral_elements_with_trace([-5,1])
            Traceback (most recent call last):
            ...
            ValueError: bounds must be positive
        """
        if C[0] < 0:
            raise ValueError, "bounds must be positive"

        if not self.is_totally_real():
            raise NotImplementedError, "exact computation of LLL reduction only implemented in the totally real case"

        Z_F = self.maximal_order()
        B = self.reduced_basis()
        T = self.reduced_gram_matrix()
        P = pari(T).qfminim((C[1]**2)*(1./2), 10**6)[2]

        S = []
        for p in P:
            theta = sum([ p.list()[i]*B[i] for i in range(self.degree())])
            if theta.trace() < 0:
                theta *= -1
            if theta.trace() >= C[0] and theta.trace() <= C[1]:
                if self(theta).is_totally_positive():
                    S.append(self(theta))
        return S


    def zeta_function(self, prec=53,
                      max_imaginary_part=0,
                      max_asymp_coeffs=40):
        r"""
        Return the Zeta function of this number field.

        This actually returns an interface to Tim Dokchitser's program for
        computing with the Dedekind zeta function zeta_F(s) of the number
        field F.

        INPUT:


        -  ``prec`` - integer (bits precision)

        -  ``max_imaginary_part`` - real number

        -  ``max_asymp_coeffs`` - integer


        OUTPUT: The zeta function of this number field.

        EXAMPLES::

            sage: K.<a> = NumberField(ZZ['x'].0^2+ZZ['x'].0-1)
            sage: Z = K.zeta_function()
            sage: Z
            Zeta function associated to Number Field in a with defining polynomial x^2 + x - 1
            sage: Z(-1)
            0.0333333333333333
            sage: L.<a, b, c> = NumberField([x^2 - 5, x^2 + 3, x^2 + 1])
            sage: Z = L.zeta_function()
            sage: Z(5)
            1.00199015670185
        """
        from sage.lfunctions.all import Dokchitser
        key = (prec, max_imaginary_part, max_asymp_coeffs)
        r1 = self.signature()[0]
        r2 = self.signature()[1]
        zero = [0]
        one = [1]
        Z = Dokchitser(conductor = abs(self.absolute_discriminant()),
                       gammaV = (r1+r2)*zero + r2*one,
                       weight = 1,
                       eps = 1,
                       poles = [1],
                       prec = prec)
        s = 'nf = nfinit(%s);'%self.absolute_polynomial()
        s += 'dzk = dirzetak(nf,cflength());'
        Z.init_coeffs('dzk[k]', pari_precode=s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        Z.check_functional_equation()
        Z.rename('Zeta function associated to %s'%self)
        return Z

    def narrow_class_group(self, proof=None):
        r"""
        Return the narrow class group of this field.

        INPUT:


        -  ``proof`` - default: None (use the global proof
           setting, which defaults to True).


        EXAMPLES::

            sage: NumberField(x^3+x+9, 'a').narrow_class_group()
            Multiplicative Abelian group isomorphic to C2
        """
        proof = proof_flag(proof)
        try:
            return self.__narrow_class_group
        except AttributeError:
            k = self.pari_bnf(proof)
            s = str(k.bnfnarrow())
            s = s.replace(";",",")
            s = eval(s)
            self.__narrow_class_group = sage.groups.abelian_gps.abelian_group.AbelianGroup(s[1])
        return self.__narrow_class_group

    def ngens(self):
        """
        Return the number of generators of this number field (always 1).

        OUTPUT: the python integer 1.

        EXAMPLES::

            sage: NumberField(x^2 + 17,'a').ngens()
            1
            sage: NumberField(x + 3,'a').ngens()
            1
            sage: k.<a> = NumberField(x + 3)
            sage: k.ngens()
            1
            sage: k.0
            -3
        """
        return 1

    def order(self):
        """
        Return the order of this number field (always +infinity).

        OUTPUT: always positive infinity

        EXAMPLES::

            sage: NumberField(x^2 + 19,'a').order()
            +Infinity
        """
        return infinity.infinity

    def absolute_polynomial_ntl(self):
        r"""
        Alias for :meth:`~polynomial_ntl`. Mostly for internal use.

        EXAMPLES::

            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').absolute_polynomial_ntl()
            ([-27 34 51], 51)
        """
        return self.polynomial_ntl()

    def polynomial_ntl(self):
        """
        Return defining polynomial of this number field as a pair, an ntl
        polynomial and a denominator.

        This is used mainly to implement some internal arithmetic.

        EXAMPLES::

            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').polynomial_ntl()
            ([-27 34 51], 51)
        """
        try:
            return (self.__polynomial_ntl, self.__denominator_ntl)
        except AttributeError:
            self.__denominator_ntl = ntl.ZZ()
            den = self.polynomial().denominator()
            self.__denominator_ntl.set_from_sage_int(ZZ(den))
            self.__polynomial_ntl = ntl.ZZX((self.polynomial()*den).list())
        return (self.__polynomial_ntl, self.__denominator_ntl)

    def polynomial(self):
        """
        Return the defining polynomial of this number field.

        This is exactly the same as
        ``self.defining_polynomal()``.

        EXAMPLES::

            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').polynomial()
            x^2 + 2/3*x - 9/17
        """
        return self.__polynomial

    def defining_polynomial(self):   # do not overload this -- overload polynomial instead
        """
        Return the defining polynomial of this number field.

        This is exactly the same as ``self.polynomal()``.

        EXAMPLES::

            sage: k5.<z> = CyclotomicField(5)
            sage: k5.defining_polynomial()
            x^4 + x^3 + x^2 + x + 1
            sage: y = polygen(QQ,'y')
            sage: k.<a> = NumberField(y^9 - 3*y + 5); k
            Number Field in a with defining polynomial y^9 - 3*y + 5
            sage: k.defining_polynomial()
            y^9 - 3*y + 5
        """
        return self.polynomial()

    def polynomial_ring(self):
        """
        Return the polynomial ring that we view this number field as being
        a quotient of (by a principal ideal).

        EXAMPLES: An example with an absolute field::

            sage: k.<a> = NumberField(x^2 + 3)
            sage: y = polygen(QQ, 'y')
            sage: k.<a> = NumberField(y^2 + 3)
            sage: k.polynomial_ring()
            Univariate Polynomial Ring in y over Rational Field

        An example with a relative field::

            sage: y = polygen(QQ, 'y')
            sage: M.<a> = NumberField([y^3 + 97, y^2 + 1]); M
            Number Field in a0 with defining polynomial y^3 + 97 over its base field
            sage: M.polynomial_ring()
            Univariate Polynomial Ring in y over Number Field in a1 with defining polynomial y^2 + 1
        """
        return self.relative_polynomial().parent()

    def polynomial_quotient_ring(self):
        """
        Return the polynomial quotient ring isomorphic to this number
        field.

        EXAMPLES::

            sage: K = NumberField(x^3 + 2*x - 5, 'alpha')
            sage: K.polynomial_quotient_ring()
            Univariate Quotient Polynomial Ring in alpha over Rational Field with modulus x^3 + 2*x - 5
        """
        return self.polynomial_ring().quotient(self.relative_polynomial(), self.variable_name())

    def regulator(self, proof=None):
        """
        Return the regulator of this number field.

        Note that PARI computes the regulator to higher precision than the
        Sage default.

        INPUT:


        -  ``proof`` - default: True, unless you set it
           otherwise.


        EXAMPLES::

            sage: NumberField(x^2-2, 'a').regulator()
            0.881373587019543
            sage: NumberField(x^4+x^3+x^2+x+1, 'a').regulator()
            0.962423650119207
        """
        proof = proof_flag(proof)
        try:
            return self.__regulator
        except AttributeError:
            from sage.rings.all import RealField
            k = self.pari_bnf(proof)
            self.__regulator = RealField(53)(k.bnf_get_reg())
            return self.__regulator

    def residue_field(self, prime, names=None, check=True):
        """
        Return the residue field of this number field at a given prime, ie
        `O_K / p O_K`.

        INPUT:


        -  ``prime`` - a prime ideal of the maximal order in
           this number field, or an element of the field which generates a
           principal prime ideal.

        -  ``names`` - the name of the variable in the residue
           field

        -  ``check`` - whether or not to check the primality of
           prime.


        OUTPUT: The residue field at this prime.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: K.residue_field(P)
            Residue field in abar of Fractional ideal (61, a^2 + 30)

        ::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.residue_field(1+i)
            Residue field of Fractional ideal (i + 1)

        TESTS::

            sage: L.<b> = NumberField(x^2 + 5)
            sage: L.residue_field(P)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (61, a^2 + 30) is not an ideal of Number Field in b with defining polynomial x^2 + 5
            sage: L.residue_field(2)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (2) is not a prime ideal

        ::

            sage: L.residue_field(L.prime_above(5)^2)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (5) is not a prime ideal
        """
        from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal
        if is_NumberFieldIdeal(prime) and prime.number_field() is not self:
            raise ValueError, "%s is not an ideal of %s"%(prime,self)
        # This allows principal ideals to be specified using a generator:
        try:
            prime = self.ideal(prime)
        except TypeError:
            pass

        if not is_NumberFieldIdeal(prime) or prime.number_field() is not self:
            raise ValueError, "%s is not an ideal of %s"%(prime,self)
        if check and not prime.is_prime():
            raise ValueError, "%s is not a prime ideal"%prime
        from sage.rings.residue_field import ResidueField
        return ResidueField(prime, names=names, check=False)

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES::

            sage: NumberField(x^2+1, 'a').signature()
            (0, 1)
            sage: NumberField(x^3-2, 'a').signature()
            (1, 1)
        """
        r1, r2 = self.pari_nf().nf_get_sign()
        return (ZZ(r1), ZZ(r2))

    def trace_pairing(self, v):
        """
        Return the matrix of the trace pairing on the elements of the list
        `v`.

        EXAMPLES::

            sage: K.<zeta3> = NumberField(x^2 + 3)
            sage: K.trace_pairing([1,zeta3])
            [ 2  0]
            [ 0 -6]
        """
        import sage.matrix.matrix_space
        A = sage.matrix.matrix_space.MatrixSpace(self.base_ring(), len(v))(0)
        for i in range(len(v)):
            for j in range(i,len(v)):
                t = (self(v[i]*v[j])).trace()
                A[i,j] = t
                A[j,i] = t
        return A

    def uniformizer(self, P, others="positive"):
        """
        Returns an element of self with valuation 1 at the prime ideal P.

        INPUT:


        -  ``self`` - a number field

        -  ``P`` - a prime ideal of self

        -  ``others`` - either "positive" (default), in which
           case the element will have non-negative valuation at all other
           primes of self, or "negative", in which case the element will have
           non-positive valuation at all other primes of self.


        .. note::

           When P is principal (e.g. always when self has class number
           one) the result may or may not be a generator of P!

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 5); K
            Number Field in a with defining polynomial x^2 + 5
            sage: P,Q = K.ideal(3).prime_factors()
            sage: P
            Fractional ideal (3, a + 1)
            sage: pi = K.uniformizer(P); pi
            a + 1
            sage: K.ideal(pi).factor()
            (Fractional ideal (2, a + 1)) * (Fractional ideal (3, a + 1))
            sage: pi = K.uniformizer(P,'negative'); pi
            1/2*a + 1/2
            sage: K.ideal(pi).factor()
            (Fractional ideal (2, a + 1))^-1 * (Fractional ideal (3, a + 1))

        ::

            sage: K = CyclotomicField(9)
            sage: Plist=K.ideal(17).prime_factors()
            sage: pilist = [K.uniformizer(P) for P in Plist]
            sage: [pi.is_integral() for pi in pilist]
            [True, True, True]
            sage: [pi.valuation(P) for pi,P in zip(pilist,Plist)]
            [1, 1, 1]
            sage: [ pilist[i] in Plist[i] for i in range(len(Plist)) ]
            [True, True, True]

        ::

            sage: K.<t> = NumberField(x^4 - x^3 - 3*x^2 - x + 1)
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(2))]
            [2]
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(3))]
            [t - 1]
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(5))]
            [t^2 - t + 1, t + 2, t - 2]
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(7))]
            [t^2 + 3*t + 1]
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(67))]
            [t + 23, t + 26, t - 32, t - 18]

        ALGORITHM:

            Use PARI. More precisely, use the second component of
            ``idealprimedec`` in the "positive" case. Use `idealappr`
            with exponent of -1 and invert the result in the "negative"
            case.
        """
        if not is_NumberFieldIdeal(P):
            P = self.ideal(P)
        P = P.pari_prime()
        if others == "positive":
            return self(P[1])
        elif others == "negative":
            nf = self.pari_nf()
            F = pari.matrix(1, 2, [P, -1])
            return ~self(nf.idealappr(F, 1))
        else:
            raise ValueError("others must be 'positive' or 'negative'")

    def units(self, proof=None):
        """
        Return generators for the unit group modulo torsion.

        INPUT:

        - ``proof`` (bool, default True) flag passed to ``pari``.

        .. note::

            For more functionality see the unit_group() function.

        ALGORITHM: Uses PARI's bnfunit command.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: K.units()
            [6/275*a^3 - 9/55*a^2 + 14/11*a - 2]

        For big number fields, provably computing the unit group can
        take a very long time.  In this case, one can ask for the
        conjectural unit group (correct if the Generalized Riemann
        Hypothesis is true)::

            sage: K = NumberField(x^17 + 3, 'a')
            sage: K.units(proof=True)  # takes forever, not tested
            ...
            sage: K.units(proof=False)  # result not independently verified
            [a^9 + a - 1, a^16 - a^15 + a^14 - a^12 + a^11 - a^10 - a^8 + a^7 - 2*a^6 + a^4 - 3*a^3 + 2*a^2 - 2*a + 1, 2*a^16 - a^14 - a^13 + 3*a^12 - 2*a^10 + a^9 + 3*a^8 - 3*a^6 + 3*a^5 + 3*a^4 - 2*a^3 - 2*a^2 + 3*a + 4, a^15 + a^14 + 2*a^11 + a^10 - a^9 + a^8 + 2*a^7 - a^5 + 2*a^3 - a^2 - 3*a + 1, a^16 + a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 - 2, 2*a^16 - 3*a^15 + 3*a^14 - 3*a^13 + 3*a^12 - a^11 + a^9 - 3*a^8 + 4*a^7 - 5*a^6 + 6*a^5 - 4*a^4 + 3*a^3 - 2*a^2 - 2*a + 4, a^15 - a^12 + a^10 - a^9 - 2*a^8 + 3*a^7 + a^6 - 3*a^5 + a^4 + 4*a^3 - 3*a^2 - 2*a + 2, 2*a^16 + a^15 - a^11 - 3*a^10 - 4*a^9 - 4*a^8 - 4*a^7 - 5*a^6 - 7*a^5 - 8*a^4 - 6*a^3 - 5*a^2 - 6*a - 7]
        """
        proof = proof_flag(proof)

        # if we have cached provable results, return them immediately
        try:
            return self.__units
        except AttributeError:
            pass

        # if proof==False and we have cached results, return them immediately
        if not proof:
            try:
                return self.__units_no_proof
            except AttributeError:
                pass

        # get PARI to compute the units
        B = self.pari_bnf(proof).bnfunit()
        if proof:
            # cache the provable results and return them
            self.__units = map(self, B)
            return self.__units
        else:
            # cache the conjectural results and return them
            self.__units_no_proof = map(self, B)
            return self.__units_no_proof

    def unit_group(self, proof=None):
        """
        Return the unit group (including torsion) of this number field.

        ALGORITHM: Uses PARI's bnfunit command.

        INPUT:

        - ``proof`` (bool, default True) flag passed to ``pari``.

        .. note::

           The group is cached.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: U = K.unit_group(); U
            Unit group with structure C10 x Z of Number Field in a with defining polynomial x^4 - 10*x^3 + 100*x^2 - 375*x + 1375
            sage: U.gens()
            (u0, u1)
            sage: U.gens_values()
            [-7/275*a^3 + 1/11*a^2 - 9/11*a - 1, 6/275*a^3 - 9/55*a^2 + 14/11*a - 2]
            sage: U.invariants()
            (10, 0)
            sage: [u.multiplicative_order() for u in U.gens()]
            [10, +Infinity]

        For big number fields, provably computing the unit group can
        take a very long time.  In this case, one can ask for the
        conjectural unit group (correct if the Generalized Riemann
        Hypothesis is true)::

            sage: K = NumberField(x^17 + 3, 'a')
            sage: K.unit_group(proof=True)  # takes forever, not tested
            ...
            sage: U = K.unit_group(proof=False)
            sage: U
            Unit group with structure C2 x Z x Z x Z x Z x Z x Z x Z x Z of Number Field in a with defining polynomial x^17 + 3
            sage: U.gens()
            (u0, u1, u2, u3, u4, u5, u6, u7, u8)
            sage: U.gens_values()  # result not independently verified
            [-1, a^9 + a - 1, a^16 - a^15 + a^14 - a^12 + a^11 - a^10 - a^8 + a^7 - 2*a^6 + a^4 - 3*a^3 + 2*a^2 - 2*a + 1, 2*a^16 - a^14 - a^13 + 3*a^12 - 2*a^10 + a^9 + 3*a^8 - 3*a^6 + 3*a^5 + 3*a^4 - 2*a^3 - 2*a^2 + 3*a + 4, a^15 + a^14 + 2*a^11 + a^10 - a^9 + a^8 + 2*a^7 - a^5 + 2*a^3 - a^2 - 3*a + 1, a^16 + a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 - 2, 2*a^16 - 3*a^15 + 3*a^14 - 3*a^13 + 3*a^12 - a^11 + a^9 - 3*a^8 + 4*a^7 - 5*a^6 + 6*a^5 - 4*a^4 + 3*a^3 - 2*a^2 - 2*a + 4, a^15 - a^12 + a^10 - a^9 - 2*a^8 + 3*a^7 + a^6 - 3*a^5 + a^4 + 4*a^3 - 3*a^2 - 2*a + 2, 2*a^16 + a^15 - a^11 - 3*a^10 - 4*a^9 - 4*a^8 - 4*a^7 - 5*a^6 - 7*a^5 - 8*a^4 - 6*a^3 - 5*a^2 - 6*a - 7]
        """
        proof = proof_flag(proof)

        try:
            return self._unit_group
        except AttributeError:
            pass

        if proof == False:
            try:
                return self._unit_group_no_proof
            except AttributeError:
                pass

        U = UnitGroup(self,proof)
        if proof:
            self._unit_group = U
        else:
            self._unit_group_no_proof = U
        return U

    def S_unit_group(self, proof=None, S=None):
        """
        Return the S-unit group (including torsion) of this number field.

        ALGORITHM: Uses PARI's bnfsunit command.

        INPUT:

        - ``proof`` (bool, default True) flag passed to ``pari``.
        - ``S`` - list or tuple of prime ideals, or an ideal, or a single
          ideal or element from which an ideal can be constructed, in
          which case the support is used.  If None, the global unit
          group is constructed; otherwise, the S-unit group is
          constructed.

        .. note::

           The group is cached.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3)
            sage: U = K.S_unit_group(S=a); U
            S-unit group with structure C10 x Z x Z x Z of Number Field in a with defining polynomial x^4 - 10*x^3 + 100*x^2 - 375*x + 1375 with S = (Fractional ideal (11, 1/275*a^3 + 4/55*a^2 - 5/11*a + 9), Fractional ideal (5, 1/275*a^3 + 4/55*a^2 - 5/11*a + 5))
            sage: U.gens()
            (u0, u1, u2, u3)
            sage: U.gens_values()
            [-7/275*a^3 + 1/11*a^2 - 9/11*a - 1, 6/275*a^3 - 9/55*a^2 + 14/11*a - 2, 14/275*a^3 - 21/55*a^2 + 29/11*a - 6, 1/275*a^3 + 4/55*a^2 - 5/11*a + 5]  # 32-bit
            [-7/275*a^3 + 1/11*a^2 - 9/11*a - 1, 6/275*a^3 - 9/55*a^2 + 14/11*a - 2, -14/275*a^3 + 21/55*a^2 - 29/11*a + 6, 1/275*a^3 + 4/55*a^2 - 5/11*a + 5]  # 64-bit
            sage: U.invariants()
            (10, 0, 0, 0)
            sage: [u.multiplicative_order() for u in U.gens()]
            [10, +Infinity, +Infinity, +Infinity]
            sage: U.primes()
            (Fractional ideal (11, 1/275*a^3 + 4/55*a^2 - 5/11*a + 9), Fractional ideal (5, 1/275*a^3 + 4/55*a^2 - 5/11*a + 5))

        With the default value of `S`, the S-unit group is the same as
        the global unit group::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.unit_group(proof=False)
            sage: U == K.S_unit_group(proof=False)
            True

        The value of `S` may be specified as a list of prime ideals,
        or an ideal, or an element of the field::

            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.S_unit_group(proof=False, S=K.ideal(6).prime_factors()); U
            S-unit group with structure C2 x Z x Z x Z x Z of Number Field in a with defining polynomial x^3 + 3 with S = (Fractional ideal (-a^2 + a - 1), Fractional ideal (a + 1), Fractional ideal (a))
            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.S_unit_group(proof=False, S=K.ideal(6)); U
            S-unit group with structure C2 x Z x Z x Z x Z of Number Field in a with defining polynomial x^3 + 3 with S = (Fractional ideal (-a^2 + a - 1), Fractional ideal (a + 1), Fractional ideal (a))
            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.S_unit_group(proof=False, S=6); U
            S-unit group with structure C2 x Z x Z x Z x Z of Number Field in a with defining polynomial x^3 + 3 with S = (Fractional ideal (-a^2 + a - 1), Fractional ideal (a + 1), Fractional ideal (a))

            sage: U
            S-unit group with structure C2 x Z x Z x Z x Z of Number Field in a with defining polynomial x^3 + 3 with S = (Fractional ideal (-a^2 + a - 1), Fractional ideal (a + 1), Fractional ideal (a))
            sage: U.primes()
            (Fractional ideal (-a^2 + a - 1),
            Fractional ideal (a + 1),
            Fractional ideal (a))
            sage: U.gens()
            (u0, u1, u2, u3, u4)
            sage: U.gens_values()
            [-1, a^2 - 2, -a^2 + a - 1, a + 1, a]

        The exp and log methods can be used to create `S`-units from
        sequences of exponents, and recover the exponents::

            sage: U.gens_orders()
            (2, 0, 0, 0, 0)
            sage: u = U.exp((3,1,4,1,5)); u
            -6*a^2 + 18*a - 54
            sage: u.norm().factor()
            -1 * 2^9 * 3^5
            sage: U.log(u)
            (1, 1, 4, 1, 5)

        """
        proof = proof_flag(proof)

        # process the parameter S:
        if not S:
            S = ()
        else:
            if type(S)==list:
                S = tuple(S)
            if not type(S)==tuple:
                try:
                    S = tuple(self.ideal(S).prime_factors())
                except (NameError, TypeError, ValueError):
                    raise ValueError("Cannot make a set of primes from %s"%(S,))
            else:
                try:
                    S = tuple(self.ideal(P) for P in S)
                except (NameError, TypeError, ValueError):
                    raise ValueError("Cannot make a set of primes from %s"%(S,))
                if not all([P.is_prime() for P in S]):
                    raise ValueError("Not all elements of %s are prime ideals"%(S,))

        try:
            return self._S_unit_group_cache[S]
        except AttributeError:
            self._S_unit_group_cache = {}
        except KeyError:
            pass

        if proof == False:
            try:
                return self._S_unit_group_no_proof_cache[S]
            except AttributeError:
                self._S_unit_group_no_proof_cache = {}
            except KeyError:
                pass

        U = UnitGroup(self,proof,S=S)
        if proof:
            self._S_unit_group_cache[S] = U
        else:
            self._S_unit_group_no_proof_cache[S] = U
        return U

    def zeta(self, n=2, all=False):
        """
        Return one, or a list of all, primitive n-th root of unity in this field.

        INPUT:

        -  ``n`` - positive integer

        - ``all`` - bool.  If False (default), return a primitive
          `n`-th root of unity in this field, or raise a ValueError
          exception if there are none.  If True, return a list of
          all primitive `n`-th roots of unity in this field
          (possibly empty).

        .. note::

           To obtain the maximal order of a root of unity in this field,
           use self.number_of_roots_of_unity().

        .. note::

           We do not create the full unit group since that can be
           expensive, but we do use it if it is already known.

        EXAMPLES::

            sage: K.<z> = NumberField(x^2 + 3)
            sage: K.zeta(1)
            1
            sage: K.zeta(2)
            -1
            sage: K.zeta(2, all=True)
            [-1]
            sage: K.zeta(3)
            -1/2*z - 1/2
            sage: K.zeta(3, all=True)
            [-1/2*z - 1/2, 1/2*z - 1/2]
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: There are no 4th roots of unity in self.

        ::

            sage: r.<x> = QQ[]
            sage: K.<b> = NumberField(x^2+1)
            sage: K.zeta(4)
            -b
            sage: K.zeta(4,all=True)
            [-b, b]
            sage: K.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: There are no 3rd roots of unity in self.
            sage: K.zeta(3,all=True)
            []
        """
        try:
            return self._unit_group.zeta(n,all)
        except AttributeError:
            pass
        try:
            return self._unit_group_no_proof.zeta(n,all)
        except AttributeError:
            pass

        K = self
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n (=%s) must be positive"%n
        if n == 1:
            if all:
                return [K(1)]
            else:
                return K(1)
        elif n == 2:
            if all:
                return [K(-1)]
            else:
                return K(-1)
        N = K.zeta_order()
        if n.divides(N):
            z = K.primitive_root_of_unity() ** (N//n)
            if all:
                return [z**i for i in n.coprime_integers(n)]
            else:
                return z
        else:
            if all:
                return []
            else:
                raise ValueError, "There are no %s roots of unity in self."%n.ordinal_str()

    def zeta_order(self):
        r"""
        Return the number of roots of unity in this field.

        .. note::

           We do not create the full unit group since that can be
           expensive, but we do use it if it is already known.

        EXAMPLES::

            sage: F.<alpha> = NumberField(x**22+3)
            sage: F.zeta_order()
            6
            sage: F.<alpha> = NumberField(x**2-7)
            sage: F.zeta_order()
            2
        """
        try:
            return self._unit_group.zeta_order()
        except AttributeError:
            pass
        try:
            return self._unit_group_no_proof.zeta_order()
        except AttributeError:
            pass

        return ZZ(self.pari_nf().nfrootsof1()[0])

    number_of_roots_of_unity = zeta_order

    def primitive_root_of_unity(self):
        """
        Return a generator of the roots of unity in this field.

        OUTPUT: a primitive root of unity. No guarantee is made about
        which primitive root of unity this returns, not even for
        cyclotomic fields.

        .. note::

           We do not create the full unit group since that can be
           expensive, but we do use it if it is already known.

        ALGORITHM:

        We use the PARI function ``nfrootsof1`` in all cases. This is
        required (even for cyclotomic fields) in order to be consistent
        with the full unit group, which is also computed by PARI.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: z = K.primitive_root_of_unity(); z
            -i
            sage: z.multiplicative_order()
            4

            sage: K.<a> = NumberField(x^2+x+1)
            sage: z = K.primitive_root_of_unity(); z
            a + 1
            sage: z.multiplicative_order()
            6

            sage: x = polygen(QQ)
            sage: F.<a,b> = NumberField([x^2 - 2, x^2 - 3])
            sage: y = polygen(F)
            sage: K.<c> = F.extension(y^2 - (1 + a)*(a + b)*a*b)
            sage: K.primitive_root_of_unity()
            -1

        We do not special-case cyclotomic fields, so we do not always
        get the most obvious primitive root of unity::

            sage: K.<a> = CyclotomicField(3)
            sage: z = K.primitive_root_of_unity(); z
            a + 1
            sage: z.multiplicative_order()
            6

            sage: K = CyclotomicField(7)
            sage: z = K.primitive_root_of_unity(); z
            zeta7^5 + zeta7^4 + zeta7^3 + zeta7^2 + zeta7 + 1
            sage: z.multiplicative_order()
            14

        TESTS:

        Check for :trac:`15027`. We use a new variable name::

            sage: K.<f> = QuadraticField(-3)
            sage: K.primitive_root_of_unity()
            -1/2*f + 1/2
            sage: UK = K.unit_group()
            sage: K.primitive_root_of_unity()
            -1/2*f + 1/2
        """
        try:
            return self._unit_group.torsion_generator().value()
        except AttributeError:
            pass
        try:
            return self._unit_group_no_proof.torsion_generator().value()
        except AttributeError:
            pass

        pK = self.pari_nf()
        n, z = pK.nfrootsof1()
        return self(z)

    def roots_of_unity(self):
        """
        Return all the roots of unity in this field, primitive or not.

        EXAMPLES::

            sage: K.<b> = NumberField(x^2+1)
            sage: zs = K.roots_of_unity(); zs
            [-b, -1, b, 1]
            sage: [ z**K.number_of_roots_of_unity() for z in zs ]
            [1, 1, 1, 1]
        """
        z = self.primitive_root_of_unity()
        n = self.zeta_order()
        return [ z**k for k in range(1, n+1) ]

    def zeta_coefficients(self, n):
        """
        Compute the first n coefficients of the Dedekind zeta function of
        this field as a Dirichlet series.

        EXAMPLE::

            sage: x = QQ['x'].0
            sage: NumberField(x^2+1, 'a').zeta_coefficients(10)
            [1, 1, 0, 1, 2, 0, 0, 1, 1, 2]
        """
        return self.pari_nf().dirzetak(n)


class NumberField_absolute(NumberField_generic):

    def __init__(self, polynomial, name, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None):
        """
        Function to initialize an absolute number field.

        EXAMPLE::

            sage: K = NumberField(x^17 + 3, 'a'); K
            Number Field in a with defining polynomial x^17 + 3
            sage: type(K)
            <class 'sage.rings.number_field.number_field.NumberField_absolute_with_category'>
            sage: TestSuite(K).run()
        """
        NumberField_generic.__init__(self, polynomial, name, latex_name, check, embedding,
                                     assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
        self._element_class = number_field_element.NumberFieldElement_absolute
        self._zero_element = self(0)
        self._one_element =  self(1)

    def _coerce_from_other_number_field(self, x):
        """
        Coerce a number field element x into this number field.

        Unless `x` is in ``QQ``, this requires ``x.parent()`` and
        ``self`` to have compatible embeddings: either they both embed
        in a common field, or there is an embedding of ``x.parent()``
        into ``self`` or the other way around.  If no compatible
        embeddings are found and `x` is not in ``QQ``, then raise
        ``TypeError``.  This guarantees that these conversions respect
        the field operations and conversions between several fields
        commute.

        REMARK:

        The name of this method was chosen for historical reasons.
        In fact, what it does is not a coercion but a conversion.

        INPUT:

        ``x`` -- an element of some number field

        OUTPUT:

        An element of ``self`` corresponding to ``x``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: L.<b> = NumberField(x^2 + 1)
            sage: K._coerce_from_other_number_field(L(2/3))
            2/3
            sage: L._coerce_from_other_number_field(K(0))
            0
            sage: K._coerce_from_other_number_field(b)
            Traceback (most recent call last):
            ...
            TypeError: No compatible natural embeddings found for Number Field in a with defining polynomial x^3 + 2 and Number Field in b with defining polynomial x^2 + 1

        Two number fields both containing `i`::

            sage: K.<a> = NumberField(x^4 + 6*x^2 + 1, embedding = CC(-2.4*I))
            sage: L.<b> = NumberField(x^4 + 8*x^2 + 4, embedding = CC(2.7*I))
            sage: Ki = 1/2*a^3 + 5/2*a; Ki.minpoly()
            x^2 + 1
            sage: L(Ki)
            -1/4*b^3 - 3/2*b
            sage: K(L(Ki)) == Ki
            True
            sage: Q.<i> = QuadraticField(-1)
            sage: Q(Ki)
            i
            sage: Q(L(Ki))
            i
            sage: L( (Ki+2)^1000 )
            737533628...075020804*b^3 + 442520177...450124824*b + 793311113...453515313

        This fails if we don't specify the embeddings::

            sage: K.<a> = NumberField(x^4 + 6*x^2 + 1)
            sage: L.<b> = NumberField(x^4 + 8*x^2 + 4)
            sage: L(1/2*a^3 + 5/2*a)
            Traceback (most recent call last):
            ...
            TypeError: No compatible natural embeddings found for Number Field in b with defining polynomial x^4 + 8*x^2 + 4 and Number Field in a with defining polynomial x^4 + 6*x^2 + 1

        Embeddings can also be `p`-adic::

            sage: F = Qp(73)
            sage: K.<a> = NumberField(x^4 + 6*x^2 + 1, embedding = F(1290990671961076190983179596556712119))
            sage: L.<b> = NumberField(x^4 + 8*x^2 + 4, embedding = F(1773398470280167815153042237103591466))
            sage: L(2*a^3 + 10*a + 3)
            b^3 + 6*b + 3

        If we take the same non-Galois number field with two different
        embeddings, conversion fails::

            sage: K.<a> = NumberField(x^3 - 4*x + 1, embedding = 0.254)
            sage: L.<b> = NumberField(x^3 - 4*x + 1, embedding = 1.86)
            sage: L(a)
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert a to Number Field in b with defining polynomial x^3 - 4*x + 1 (using the specified embeddings)

        Subfields automatically come with an embedding::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: L.<b>, phi = K.subfield(-a)
            sage: phi(b)
            -a
            sage: K(b)
            -a
            sage: L(a)
            -b

        Below we create two subfields of `K` which both contain `i`.
        Since `L2` and `L3` both embed in `K`, conversion works::

            sage: K.<z> = NumberField(x^8 - x^4 + 1)
            sage: i = (x^2+1).roots(ring=K)[0][0]
            sage: r2 = (x^2-2).roots(ring=K)[0][0]
            sage: r3 = (x^2-3).roots(ring=K)[0][0]
            sage: L2.<a2>, phi2 = K.subfield(r2+i)
            sage: L3.<a3>, phi3 = K.subfield(r3+i)
            sage: i_in_L2 = L2(i); i_in_L2
            1/6*a2^3 + 1/6*a2
            sage: i_in_L3 = L3(i); i_in_L3
            1/8*a3^3
            sage: L2(i_in_L3) == i_in_L2
            True
            sage: L3(i_in_L2) == i_in_L3
            True

        TESTS:

        The following was fixed in :trac:`8800`::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-5,embedding=0)
            sage: L.<b> = K.extension(x^2+a)
            sage: F,R = L.construction()
            sage: F(R) == L   #indirect doctest
            True

        AUTHORS:

        - Jeroen Demeyer (2011-09-30): :trac:`11869`

        """
        # Special case for x in QQ.  This is common, so should be fast.
        xpol = x.polynomial()
        if xpol.degree() <= 0:
            return self._element_class(self, xpol[0])
        # Convert from L to K
        K = self
        L = x.parent()
        # Find embeddings for K and L.  If no embedding is given, simply
        # take the identity map as "embedding".  This handles the case
        # where one field is created as subfield of the other.
        Kgen = K.gen_embedding()
        if Kgen is None:
            Kgen = K.gen()
        KF = Kgen.parent()
        Lgen = L.gen_embedding()
        if Lgen is None:
            Lgen = L.gen()
        LF = Lgen.parent()

        # Do not use CDF or RDF because of constraints on the
        # exponent of floating-point numbers
        from sage.rings.all import RealField, ComplexField
        CC = ComplexField(53)
        RR = RealField(53)

        # Find a common field F into which KF and LF both embed.
        if CC.has_coerce_map_from(KF) and CC.has_coerce_map_from(LF):
            # We postpone converting Kgen and Lgen to F until we know the
            # floating-point precision required.
            F = CC
        elif KF is LF:
            F = KF
        elif KF.has_coerce_map_from(LF):
            F = KF
            Lgen = F(Lgen)
        elif LF.has_coerce_map_from(KF):
            F = LF
            Kgen = F(Kgen)
        else:
            raise TypeError, "No compatible natural embeddings found for %s and %s"%(KF,LF)

        # List of candidates for K(x)
        f = x.minpoly()
        ys = f.roots(ring=K, multiplicities=False)
        if not ys:
            raise ValueError, "Cannot convert %s to %s (regardless of embeddings)"%(x,K)

        # Define a function are_roots_equal to determine whether two
        # roots of f are equal.  A simple a == b does not suffice for
        # inexact fields because of floating-point errors.
        if F.is_exact():
            are_roots_equal = lambda a,b: a == b
        else:
            ### Compute a lower bound on the distance between the roots of f.
            ### This essentially gives the precision to work with.

            # A function
            # log2abs: F --> RR
            #          x |-> log2(abs(x))
            # This should work for all fields F with an absolute value.
            # The p-adic absolute value goes into QQ, so we need the RR().
            log2abs = lambda x: RR(F(x).abs()).log2()

            # Compute half Fujiwara's bound on the roots of f
            n = f.degree()
            log_half_root_bound = log2abs(f[0]/2)/n
            for i in range(1,n):
                bd = log2abs(f[i])/(n-i)
                if bd > log_half_root_bound:
                    log_half_root_bound = bd
            # Twice the bound on the roots of f, in other words an upper
            # bound for the distance between two roots.
            log_double_root_bound = log_half_root_bound + 2.0  # 2.0 = log2(4)
            # Now we compute the minimum distance between two roots of f
            # using the fact that the discriminant of f is the product of
            # all root distances.
            # We use pari to compute the discriminant to work around #11872.
            log_root_diff = log2abs(pari(f).poldisc())*0.5 - (n*(n-1)*0.5 - 1.0)*log_double_root_bound
            # Let eps be 1/128 times the minimal root distance.
            # This implies: If two roots of f are at distance <= eps, then
            # they are equal.  The factor 128 is arbitrary, it is an extra
            # safety margin.
            eps = (log_root_diff - 7.0).exp2()
            are_roots_equal = lambda a,b: (a-b).abs() <= eps
            if F is CC:
                # Adjust the precision of F, sufficient to represent all
                # the temporaries in the computation with a precision
                # of eps, plus some extra bits.
                H = [log_double_root_bound - 1.0]
                for e in [x] + ys:
                    H += [log2abs(c) for c in e.polynomial().coefficients()]
                prec = (max(H) + RR(n+1).log2() - log_root_diff).ceil() + 12 + n
                F = ComplexField(prec=prec)
                Kgen = F(Kgen)
                Lgen = F(Lgen)

        # Embed x and the y's in F
        emb_x = x.polynomial()(Lgen)
        for y in ys:
            emb_y = y.polynomial()(Kgen)
            if are_roots_equal(emb_x, emb_y):
                return y
        raise ValueError, "Cannot convert %s to %s (using the specified embeddings)"%(x,K)

    def _coerce_non_number_field_element_in(self, x):
        """
        Coerce a non-number field element x into this number field.

        INPUT:
            x -- a non number field element x, e.g., a list, integer,
            rational, or polynomial.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2/3)
            sage: K._coerce_non_number_field_element_in(-7/8)
            -7/8
            sage: K._coerce_non_number_field_element_in([1,2,3])
            3*a^2 + 2*a + 1

        The list is just turned into a polynomial in the generator.
            sage: K._coerce_non_number_field_element_in([0,0,0,1,1])
            -2/3*a - 2/3

        Any polynomial whose coefficients can be coerced to rationals will
        coerce, e.g., this one in characteristic 7.
            sage: f = GF(7)['y']([1,2,3]); f
            3*y^2 + 2*y + 1
            sage: K._coerce_non_number_field_element_in(f)
            3*a^2 + 2*a + 1

        But not this one over a field of order 27.
            sage: F27.<g> = GF(27)
            sage: f = F27['z']([g^2, 2*g, 1]); f
            z^2 + 2*g*z + g^2
            sage: K._coerce_non_number_field_element_in(f)
            Traceback (most recent call last):
            ...
            TypeError: <type 'sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX'>

        One can also coerce an element of the polynomial quotient ring
        that's isomorphic to the number field:
            sage: K.<a> = NumberField(x^3 + 17)
            sage: b = K.polynomial_quotient_ring().random_element()
            sage: K(b)
            -1/2*a^2 - 4
        """
        if isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              list)):
            return self._element_class(self, x)

        if isinstance(x, sage.rings.polynomial.polynomial_quotient_ring_element.PolynomialQuotientRingElement)\
               and (x in self.polynomial_quotient_ring()):
            y = self.polynomial_ring().gen()
            return x.lift().subs({y:self.gen()})

        try:
            if isinstance(x, polynomial_element.Polynomial):
                return self._element_class(self, x)

            return self._element_class(self, x._rational_())
        except (TypeError, AttributeError), msg:
            pass
        raise TypeError, type(x)

    def _coerce_from_str(self, x):
        r"""
        Coerce a string representation of an element of this
        number field into this number field.

        INPUT:
            x -- string

        EXAMPLES:
            sage: k.<theta25> = NumberField(x^3+(2/3)*x+1)
            sage: k._coerce_from_str('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1

        This function is called by the coerce method when it gets a string
        as input:
            sage: k('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1
        """
        w = sage.misc.all.sage_eval(x,locals=self.gens_dict())
        if not (is_Element(w) and w.parent() is self):
            return self(w)
        else:
            return w

    def _coerce_map_from_(self, R):
        """
        Canonical coercion of a ring R into self.

        Currently any ring coercing into the base ring canonically coerces
        into this field, as well as orders in any number field coercing into
        this field, and of course the field itself as well.

        Two embedded number fields may mutually coerce into each other, if
        the pushout of the two ambient fields exists and if it is possible
        to construct an :class:`~sage.rings.number_field.number_field_morphisms.EmbeddedNumberFieldMorphism`.

        EXAMPLES::

            sage: S.<y> = NumberField(x^3 + x + 1)
            sage: S.coerce(int(4)) # indirect doctest
            4
            sage: S.coerce(long(7))
            7
            sage: S.coerce(-Integer(2))
            -2
            sage: z = S.coerce(-7/8); z, type(z)
            (-7/8, <type 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>)
            sage: S.coerce(y) is y
            True

        Fields with embeddings into an ambient field coerce naturally by the given embedding::

            sage: CyclotomicField(15).coerce(CyclotomicField(5).0 - 17/3)
            zeta15^3 - 17/3
            sage: K.<a> = CyclotomicField(16)
            sage: K(CyclotomicField(4).0)
            a^4
            sage: QuadraticField(-3, 'a').coerce_map_from(CyclotomicField(3))
            Generic morphism:
              From: Cyclotomic Field of order 3 and degree 2
              To:   Number Field in a with defining polynomial x^2 + 3
              Defn: zeta3 -> 1/2*a - 1/2

        Two embedded number fields with mutual coercions (testing against a
        bug that was fixed in :trac:`8800`)::

            sage: K.<r4> = NumberField(x^4-2)
            sage: L1.<r2_1> = NumberField(x^2-2, embedding = r4**2)
            sage: L2.<r2_2> = NumberField(x^2-2, embedding = -r4**2)
            sage: r2_1+r2_2    # indirect doctest
            0
            sage: (r2_1+r2_2).parent() is L1
            True
            sage: (r2_2+r2_1).parent() is L2
            True

        Coercion of an order (testing against a bug that was fixed in
        :trac:`8800`)::

            sage: K.has_coerce_map_from(L1)
            True
            sage: L1.has_coerce_map_from(K)
            False
            sage: K.has_coerce_map_from(L1.maximal_order())
            True
            sage: L1.has_coerce_map_from(K.maximal_order())
            False

        There are situations for which one might imagine conversion
        could make sense (at least after fixing choices), but of course
        there will be no coercion from the Symbolic Ring to a Number Field::

            sage: K.<a> = QuadraticField(2)
            sage: K.coerce(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Symbolic Ring to Number Field in a with defining polynomial x^2 - 2

        TESTS::

            sage: K.<a> = NumberField(polygen(QQ)^3-2)
            sage: type(K.coerce_map_from(QQ))
            <type 'sage.structure.coerce_maps.DefaultConvertMap_unique'>

        Make sure we still get our optimized morphisms for special fields::

            sage: K.<a> = NumberField(polygen(QQ)^2-2)
            sage: type(K.coerce_map_from(QQ))
            <type 'sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element'>

        """
        if R in [int, long, ZZ, QQ, self.base()]:
            return self._generic_convert_map(R)
        from sage.rings.number_field.order import is_NumberFieldOrder
        if is_NumberFieldOrder(R) and self.has_coerce_map_from(R.number_field()):
            return self._generic_convert_map(R)
        if is_NumberField(R) and R != QQ:
            if R.coerce_embedding() is not None:
                if self.coerce_embedding() is not None:
                    try:
                        from sage.categories.pushout import pushout
                        ambient_field = pushout(R.coerce_embedding().codomain(), self.coerce_embedding().codomain())
                        if ambient_field is not None:
                            try:
                                # the original ambient field
                                return number_field_morphisms.EmbeddedNumberFieldMorphism(R, self, ambient_field)
                            except ValueError: # no embedding found
                                # there might be one in the alg. completion
                                return number_field_morphisms.EmbeddedNumberFieldMorphism(R, self, ambient_field.algebraic_closure() if hasattr(ambient_field,'algebraic_closure') else ambient_field)
                    except (ValueError, TypeError, NotImplementedError, sage.structure.coerce_exceptions.CoercionException),msg:
                        # no success with the pushout
                        try:
                            return number_field_morphisms.EmbeddedNumberFieldMorphism(R, self)
                        except (TypeError, ValueError):
                            pass
                else:
                    # R is embedded, self isn't. So, we could only have
                    # the forgetful coercion. But this yields to non-commuting
                    # coercions, as was pointed out at ticket #8800
                    return None

    def base_field(self):
        """
        Returns the base field of self, which is always QQ

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.base_field()
            Rational Field
        """
        return QQ

    def is_absolute(self):
        """
        Returns True since self is an absolute field.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.is_absolute()
            True
        """
        return True

    def absolute_polynomial(self):
        r"""
        Return absolute polynomial that defines this absolute field. This
        is the same as ``self.polynomial()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.absolute_polynomial ()
            x^2 + 1
        """
        return self.polynomial()

    def absolute_generator(self):
        r"""
        An alias for
        :meth:`sage.rings.number_field.number_field.NumberField_generic.gen`.
        This is provided for consistency with relative fields, where the
        element returned by
        :meth:`sage.rings.number_field.number_field_rel.NumberField_relative.gen`
        only generates the field over its base field (not necessarily over
        `\QQ`).

        EXAMPLE::

            sage: K.<a> = NumberField(x^2 - 17)
            sage: K.absolute_generator()
            a
        """
        return self.gen()


    def __reduce__(self):
        """
        TESTS::

            sage: Z = var('Z')
            sage: K.<w> = NumberField(Z^3 + Z + 1)
            sage: L = loads(dumps(K))
            sage: L
            Number Field in w with defining polynomial Z^3 + Z + 1
            sage: print L == K
            True
        """
        return NumberField_absolute_v1, (self.polynomial(), self.variable_name(), self.latex_variable_name(), self.gen_embedding())

    def optimized_representation(self, names=None, both_maps=True):
        """
        Return a field isomorphic to self with a better defining polynomial
        if possible, along with field isomorphisms from the new field to
        self and from self to the new field.

        EXAMPLES: We construct a compositum of 3 quadratic fields, then
        find an optimized representation and transform elements back and
        forth.

        ::

            sage: K = NumberField([x^2 + p for p in [5, 3, 2]],'a').absolute_field('b'); K
            Number Field in b with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
            sage: L, from_L, to_L = K.optimized_representation()
            sage: L    # your answer may different, since algorithm is random
            Number Field in a14 with defining polynomial x^8 + 4*x^6 + 7*x^4 + 36*x^2 + 81
            sage: to_L(K.0)   # random
            4/189*a14^7 - 1/63*a14^6 + 1/27*a14^5 + 2/9*a14^4 - 5/27*a14^3 + 8/9*a14^2 + 3/7*a14 + 3/7
            sage: from_L(L.0)   # random
            1/1152*a1^7 + 1/192*a1^6 + 23/576*a1^5 + 17/96*a1^4 + 37/72*a1^3 + 5/6*a1^2 + 55/24*a1 + 3/4

        The transformation maps are mutually inverse isomorphisms.

        ::

            sage: from_L(to_L(K.0))
            b
            sage: to_L(from_L(L.0))     # random
            a14
        """
        return self.optimized_subfields(degree=self.degree(), name=names, both_maps=both_maps)[0]

    def optimized_subfields(self, degree=0, name=None, both_maps=True):
        """
        Return optimized representations of many (but *not* necessarily
        all!) subfields of self of the given degree, or of all possible degrees if
        degree is 0.

        EXAMPLES::

            sage: K = NumberField([x^2 + p for p in [5, 3, 2]],'a').absolute_field('b'); K
            Number Field in b with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
            sage: L = K.optimized_subfields(name='b')
            sage: L[0][0]
            Number Field in b0 with defining polynomial x - 1
            sage: L[1][0]
            Number Field in b1 with defining polynomial x^2 - 3*x + 3
            sage: [z[0] for z in L]          # random -- since algorithm is random
            [Number Field in b0 with defining polynomial x - 1,
             Number Field in b1 with defining polynomial x^2 - x + 1,
             Number Field in b2 with defining polynomial x^4 - 5*x^2 + 25,
             Number Field in b3 with defining polynomial x^4 - 2*x^2 + 4,
             Number Field in b4 with defining polynomial x^8 + 4*x^6 + 7*x^4 + 36*x^2 + 81]

        We examine one of the optimized subfields in more detail::

            sage: M, from_M, to_M = L[2]
            sage: M                             # random
            Number Field in b2 with defining polynomial x^4 - 5*x^2 + 25
            sage: from_M     # may be slightly random
            Ring morphism:
              From: Number Field in b2 with defining polynomial x^4 - 5*x^2 + 25
              To:   Number Field in a1 with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
              Defn: b2 |--> -5/1152*a1^7 + 1/96*a1^6 - 97/576*a1^5 + 17/48*a1^4 - 95/72*a1^3 + 17/12*a1^2 - 53/24*a1 - 1

        The to_M map is None, since there is no map from K to M::

            sage: to_M

        We apply the from_M map to the generator of M, which gives a
        rather large element of `K`::

            sage: from_M(M.0)          # random
            -5/1152*a1^7 + 1/96*a1^6 - 97/576*a1^5 + 17/48*a1^4 - 95/72*a1^3 + 17/12*a1^2 - 53/24*a1 - 1

        Nevertheless, that large-ish element lies in a degree 4 subfield::

            sage: from_M(M.0).minpoly()   # random
            x^4 - 5*x^2 + 25
        """
        return self._subfields_helper(degree=degree,name=name,
                                      both_maps=both_maps,optimize=True)

    def change_names(self, names):
        r"""
        Return number field isomorphic to self but with the given generator
        name.

        INPUT:


        -  ``names`` - should be exactly one variable name.


        Also, ``K.structure()`` returns from_K and to_K,
        where from_K is an isomorphism from K to self and to_K is an
        isomorphism from self to K.

        EXAMPLES::

            sage: K.<z> = NumberField(x^2 + 3); K
            Number Field in z with defining polynomial x^2 + 3
            sage: L.<ww> = K.change_names()
            sage: L
            Number Field in ww with defining polynomial x^2 + 3
            sage: L.structure()[0]
            Isomorphism given by variable name change map:
              From: Number Field in ww with defining polynomial x^2 + 3
              To:   Number Field in z with defining polynomial x^2 + 3
            sage: L.structure()[0](ww + 5/3)
            z + 5/3
        """
        return self.absolute_field(names)

    def subfields(self, degree=0, name=None):
        """
        Return all subfields of self of the given degree,
        or of all possible degrees if degree is 0.  The subfields are returned as
        absolute fields together with an embedding into self.  For the case of the
        field itself, the reverse isomorphism is also provided.

        EXAMPLES::

            sage: K.<a> = NumberField( [x^3 - 2, x^2 + x + 1] )
            sage: K = K.absolute_field('b')
            sage: S = K.subfields()
            sage: len(S)
            6
            sage: [k[0].polynomial() for k in S]
            [x - 3,
             x^2 - 3*x + 9,
             x^3 - 3*x^2 + 3*x + 1,
             x^3 - 3*x^2 + 3*x + 1,
             x^3 - 3*x^2 + 3*x - 17,
             x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1]
            sage: R.<t> = QQ[]
            sage: L = NumberField(t^3 - 3*t + 1, 'c')
            sage: [k[1] for k in L.subfields()]
            [Ring morphism:
              From: Number Field in c0 with defining polynomial t
              To:   Number Field in c with defining polynomial t^3 - 3*t + 1
              Defn: 0 |--> 0,
             Ring morphism:
              From: Number Field in c1 with defining polynomial t^3 - 3*t + 1
              To:   Number Field in c with defining polynomial t^3 - 3*t + 1
              Defn: c1 |--> c]
            sage: len(L.subfields(2))
            0
            sage: len(L.subfields(1))
            1
        """
        return self._subfields_helper(degree=degree, name=name,
                                      both_maps=True, optimize=False)

    def _subfields_helper(self, degree=0, name=None, both_maps=True, optimize=False):
        """
        Internal function: common code for optimized_subfields() and subfields().

        TESTS:

        Let's make sure embeddings are being respected::

            sage: K.<a> = NumberField(x^4 - 23, embedding=50)
            sage: K, CDF(a)
            (Number Field in a with defining polynomial x^4 - 23, 2.18993870309)
            sage: Ss = K.subfields(); len(Ss) # indirect doctest
            3
            sage: diffs = [ S.coerce_embedding()(S.gen()) - CDF(S_into_K(S.gen())) for S, S_into_K, _ in Ss ]
            sage: all(abs(diff) < 1e-5 for diff in diffs)
            True

            sage: L1, _, _ = K.subfields(2)[0]; L1, CDF(L1.gen()) # indirect doctest
            (Number Field in a0 with defining polynomial x^2 - 23, -4.79583152331)

            If we take a different embedding of the large field, we get a
            different embedding of the degree 2 subfield:

            sage: K.<a> = NumberField(x^4 - 23, embedding=-50)
            sage: L2, _, _ = K.subfields(2)[0]; L2, CDF(L2.gen()) # indirect doctest
            (Number Field in a0 with defining polynomial x^2 - 23, -4.79583152331)
        """
        if name is None:
            name = self.variable_names()
        name = sage.structure.parent_gens.normalize_names(1, name)[0]
        try:
            return self.__subfields[name, degree, both_maps, optimize]
        except AttributeError:
            self.__subfields = {}
        except KeyError:
            pass
        f = pari(self.polynomial())
        if optimize:
            v = f.polred(2)
            elts = v[0]
            polys = v[1]
        else:
            v = f.nfsubfields(degree)
            elts = [x[1] for x in v]
            polys = [x[0] for x in v]

        R = self.polynomial_ring()

        embedding = None
        ans = []
        for i in range(len(elts)):
            f = R(polys[i])
            if not (degree == 0 or f.degree() == degree):
                continue
            a = self(elts[i])
            if self.coerce_embedding() is not None:
                embedding = self.coerce_embedding()(a)
            K = NumberField(f, names=name + str(i), embedding=embedding)

            from_K = K.hom([a])    # check=False here ??   would be safe unless there are bugs.

            if both_maps and K.degree() == self.degree():
                g = K['x'](self.polynomial())
                v = g.roots()
                a = from_K(K.gen())
                for i in range(len(v)):
                    r = g.roots()[i][0]
                    to_K = self.hom([r])    # check=False here ??
                    if to_K(a) == K.gen():
                        break
            else:
                to_K = None
            ans.append((K, from_K, to_K))
        ans = Sequence(ans, immutable=True, cr=ans!=[])
        self.__subfields[name, degree, both_maps, optimize] = ans
        return ans


    def maximal_order(self, v=None):
        """
        Return the maximal order, i.e., the ring of integers, associated to
        this number field.

        INPUT:


        -  ``v`` - (default: None) None, a prime, or a list of
           primes.

           - if v is None, return the maximal order.

           - if v is a prime, return an order that is p-maximal.

           - if v is a list, return an order that is maximal at each
             prime in the list v.


        EXAMPLES: In this example, the maximal order cannot be generated by
        a single element.

        ::

            sage: k.<a> = NumberField(x^3 + x^2 - 2*x+8)
            sage: o = k.maximal_order()
            sage: o
            Maximal Order in Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8

        We compute `p`-maximal orders for several `p`. Note
        that computing a `p`-maximal order is much faster in
        general than computing the maximal order::

            sage: p = next_prime(10^22); q = next_prime(10^23)
            sage: K.<a> = NumberField(x^3 - p*q)
            sage: K.maximal_order([3]).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]
            sage: K.maximal_order([2]).basis()
            [1, a, a^2]
            sage: K.maximal_order([p]).basis()
            [1, a, a^2]
            sage: K.maximal_order([q]).basis()
            [1, a, a^2]
            sage: K.maximal_order([p,3]).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

        An example with bigger discriminant::

            sage: p = next_prime(10^97); q = next_prime(10^99)
            sage: K.<a> = NumberField(x^3 - p*q)
            sage: K.maximal_order(prime_range(10000)).basis()
            [1, a, a^2]
        """
        v = self._normalize_prime_list(v)

        try:
            return self.__maximal_order[v]
        except AttributeError:
            self.__maximal_order = {}
        except KeyError:
            pass

        B = map(self, self._pari_integral_basis(v=v))

        if len(v) == 0 or v is None:
            is_maximal = True
        else:
            is_maximal = False

        import sage.rings.number_field.order as order
        O = order.absolute_order_from_module_generators(B,
                 check_integral=False, check_rank=False,
                 check_is_ring=False, is_maximal=is_maximal)

        self.__maximal_order[v] = O
        return O

    def order(self, *args, **kwds):
        r"""
        Return the order with given ring generators in the maximal order of
        this number field.

        INPUT:

        -  ``gens`` - list of elements of self; if no
           generators are given, just returns the cardinality of this number
           field (oo) for consistency.

        -  ``check_is_integral`` - bool (default: True),
           whether to check that each generator is integral.

        -  ``check_rank`` - bool (default: True), whether to
           check that the ring generated by gens is of full rank.

        -  ``allow_subfield`` - bool (default: False), if True
           and the generators do not generate an order, i.e., they generate a
           subring of smaller rank, instead of raising an error, return an
           order in a smaller number field.


        EXAMPLES::

            sage: k.<i> = NumberField(x^2 + 1)
            sage: k.order(2*i)
            Order in Number Field in i with defining polynomial x^2 + 1
            sage: k.order(10*i)
            Order in Number Field in i with defining polynomial x^2 + 1
            sage: k.order(3)
            Traceback (most recent call last):
            ...
            ValueError: the rank of the span of gens is wrong
            sage: k.order(i/2)
            Traceback (most recent call last):
            ...
            ValueError: each generator must be integral

        Alternatively, an order can be constructed by adjoining elements to
        `\ZZ`:

            sage: K.<a> = NumberField(x^3 - 2)
            sage: ZZ[a]
            Order in Number Field in a0 with defining polynomial x^3 - 2

        TESTS:

        We verify that :trac:`2480` is fixed::

            sage: K.<a> = NumberField(x^4 + 4*x^2 + 2)
            sage: B = K.integral_basis()
            sage: K.order(*B)
            Order in Number Field in a with defining polynomial x^4 + 4*x^2 + 2
            sage: K.order(B)
            Order in Number Field in a with defining polynomial x^4 + 4*x^2 + 2
            sage: K.order(gens=B)
            Order in Number Field in a with defining polynomial x^4 + 4*x^2 + 2
        """
        # set gens appropriately from the arguments
        gens = kwds.pop('gens', args)

        if len(gens) == 0:
            return NumberField_generic.order(self)
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens = map(self, gens)
        import sage.rings.number_field.order as order
        return order.absolute_order_from_ring_generators(gens, **kwds)

    def vector_space(self):
        """
        Return a vector space V and isomorphisms self --> V and V --> self.

        OUTPUT:


        -  ``V`` - a vector space over the rational numbers

        -  ``from_V`` - an isomorphism from V to self

        -  ``to_V`` - an isomorphism from self to V


        EXAMPLES::

            sage: k.<a> = NumberField(x^3 + 2)
            sage: V, from_V, to_V  = k.vector_space()
            sage: from_V(V([1,2,3]))
            3*a^2 + 2*a + 1
            sage: to_V(1 + 2*a + 3*a^2)
            (1, 2, 3)
            sage: V
            Vector space of dimension 3 over Rational Field
            sage: to_V
            Isomorphism map:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Vector space of dimension 3 over Rational Field
            sage: from_V(to_V(2/3*a - 5/8))
            2/3*a - 5/8
            sage: to_V(from_V(V([0,-1/7,0])))
            (0, -1/7, 0)
        """
        try:
            return self.__vector_space
        except AttributeError:
            V = QQ**self.degree()
            from_V = maps.MapVectorSpaceToNumberField(V, self)
            to_V   = maps.MapNumberFieldToVectorSpace(self, V)
            self.__vector_space = (V, from_V, to_V)
            return self.__vector_space

    def absolute_vector_space(self):
        r"""
        Return vector space over `\QQ` corresponding to this
        number field, along with maps from that space to this number field
        and in the other direction.

        For an absolute extension this is identical to
        ``self.vector_space()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 5)
            sage: K.absolute_vector_space()
            (Vector space of dimension 3 over Rational Field,
             Isomorphism map:
              From: Vector space of dimension 3 over Rational Field
              To:   Number Field in a with defining polynomial x^3 - 5,
             Isomorphism map:
              From: Number Field in a with defining polynomial x^3 - 5
              To:   Vector space of dimension 3 over Rational Field)
        """
        return self.vector_space()


    def _galois_closure_and_embedding(self, names=None):
        r"""
        Return number field `K` that is the Galois closure of self and an
        embedding of self into `K`.

        INPUT:

        - ``names`` - variable name for Galois closure

        .. warning::

           This is an internal function; see :meth:`galois_closure`.

        EXAMPLES:

        For medium-sized galois groups of fields with small discriminants,
        this computation is feasible::

            sage: K.<a> = NumberField(x^6 + 4*x^2 + 2)
            sage: K.galois_group(type='pari').order()
            48
            sage: L, phi = K._galois_closure_and_embedding('c')  # long time (25s on sage.math, 2011)
            sage: phi.domain() is K, phi.codomain() is L  # long time
            (True, True)
            sage: L  # long time
            Number Field in c with defining polynomial x^48 + 672*x^44 - 29904*x^42 + 5573568*x^40 - 71988672*x^38 + 5657686832*x^36 + 24204531456*x^34 + 12891821550720*x^32 - 696282917339072*x^30 + 12802184716989696*x^28 + 987385486040115456*x^26 + 39174898963577334624*x^24 + 8992357198620665856*x^22 + 8323003980399007710720*x^20 - 250984831575663605326592*x^18 + 23158113523989042998289408*x^16 + 449057170502643207256722432*x^14 + 28571480781210190985371945728*x^12 - 14430449019771278412479533056*x^10 + 11623371608161089275002854217728*x^8 - 101534009256831666166090628596736*x^6 + 2206697807875188993828211869560832*x^4 - 4684747774733973147488568884219904*x^2 + 2708520329592370027664581735403776
            sage: K.defining_polynomial()( phi(K.gen()) )  # long time
            0
        """
        try:
            # compose with variable renaming
            L = self.__galois_closure.change_names(names)
            L_to_orig, orig_to_L = L.structure()
            # "flatten" the composition by hand
            self_into_L = self.hom([ (orig_to_L * self.__galois_closure_embedding)(self.gen()) ])
            return (L, self_into_L)
        except AttributeError:
            pass

        G = self.galois_group(type='pari')
        K = self
        self_into_K = maps.IdentityMap(self)

        while K.degree() < G.order():
            # take the one of largest degree
            newK, K_into_newK, _, _ = K.composite_fields(self, names=names, both_maps=True, preserve_embedding=False)[-1]
            if newK.degree() <= K.degree():
                raise ValueError, "Compositum field degree did not increase"
            self_into_K = K_into_newK * self_into_K
            K = newK

        L = K.change_names(names)
        L_to_orig, orig_to_L = L.structure()
        self_into_L = orig_to_L * self_into_K
        self_into_L = self.hom([ self_into_L(self.gen()) ]) # "flatten" the composition by hand

        self.__galois_closure = L
        self.__galois_closure_embedding = self_into_L
        return (self.__galois_closure, self.__galois_closure_embedding)

    def galois_closure(self, names=None, map=False):
        """
        Return number field `K` that is the Galois closure of self,
        i.e., is generated by all roots of the defining polynomial of
        self, and possibly an embedding of self into `K`.

        INPUT:

        - ``names`` - variable name for Galois closure

        - ``map`` - (default: False) also return an embedding of self into `K`

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 2)
            sage: M = K.galois_closure('b'); M
            Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
            sage: L.<a2> = K.galois_closure(); L
            Number Field in a2 with defining polynomial x^8 + 28*x^4 + 2500
            sage: K.galois_group(names=("a3")).order()
            8

        ::

            sage: phi = K.embeddings(L)[0]
            sage: phi(K.0)
            1/120*a2^5 + 19/60*a2
            sage: phi(K.0).minpoly()
            x^4 - 2

            sage: L, phi = K.galois_closure('b', map=True)
            sage: L
            Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
            sage: phi
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 - 2
              To:   Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
              Defn: a |--> 1/240*b^5 - 41/120*b

        TESTS:

        Let's make sure we're renaming correctly::

            sage: L, phi = K.galois_closure('cc', map=True)
            sage: L
            Number Field in cc with defining polynomial x^8 + 28*x^4 + 2500
            sage: phi
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 - 2
              To:   Number Field in cc with defining polynomial x^8 + 28*x^4 + 2500
              Defn: a |--> 1/240*cc^5 - 41/120*cc
        """
        L, self_into_L = self._galois_closure_and_embedding(names)
        if map:
            return (L, self_into_L)
        else:
            return L

    def automorphisms(self):
        r"""
        Compute all Galois automorphisms of self.

        This uses PARI's ``nfgaloisconj`` and is much faster than root finding
        for many fields.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 10000)
            sage: K.automorphisms()
            [
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 10000
              Defn: a |--> a,
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 10000
              Defn: a |--> -a
            ]

        Here's a larger example, that would take some time if we found
        roots instead of using PARI's specialized machinery::

            sage: K = NumberField(x^6 - x^4 - 2*x^2 + 1, 'a')
            sage: len(K.automorphisms())
            2

        `L` is the Galois closure of `K`::

            sage: L = NumberField(x^24 - 84*x^22 + 2814*x^20 - 15880*x^18 - 409563*x^16 - 8543892*x^14 + 25518202*x^12 + 32831026956*x^10 - 672691027218*x^8 - 4985379093428*x^6 + 320854419319140*x^4 + 817662865724712*x^2 + 513191437605441, 'a')
            sage: len(L.automorphisms())
            24
        """
        try:
            # this should be concordant with embeddings
            return self.__embeddings[self]
        except AttributeError:
            self.__embeddings = {}
        except KeyError:
            pass
        embs = map(self, self.pari_nf().nfgaloisconj())
        embs.sort()
        v = [ self.hom([ e ], check=False) for e in embs ]
        put_natural_embedding_first(v)
        self.__embeddings[self] = Sequence(v, cr=(v != []), immutable=True,
                                        check=False, universe=self.Hom(self))
        return self.__embeddings[self]

    def embeddings(self, K):
        """
        Compute all field embeddings of self into the field K (which need
        not even be a number field, e.g., it could be the complex numbers).
        This will return an identical result when given K as input again.

        If possible, the most natural embedding of self into K is put first
        in the list.

        INPUT:


        -  ``K`` - a number field


        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<a1> = K.galois_closure(); L
            Number Field in a1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: K.embeddings(L)[0]
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in a1 with defining polynomial x^6 + 40*x^3 + 1372
              Defn: a |--> 1/84*a1^4 + 13/42*a1
            sage: K.embeddings(L)  is K.embeddings(L)
            True

        We embed a quadratic field into a cyclotomic field::

            sage: L.<a> = QuadraticField(-7)
            sage: K = CyclotomicField(7)
            sage: L.embeddings(K)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^2 + 7
              To:   Cyclotomic Field of order 7 and degree 6
              Defn: a |--> 2*zeta7^4 + 2*zeta7^2 + 2*zeta7 + 1,
            Ring morphism:
              From: Number Field in a with defining polynomial x^2 + 7
              To:   Cyclotomic Field of order 7 and degree 6
              Defn: a |--> -2*zeta7^4 - 2*zeta7^2 - 2*zeta7 - 1
            ]

        We embed a cubic field in the complex numbers::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: K.embeddings(CC)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Complex Field with 53 bits of precision
              Defn: a |--> -0.62996052494743... - 1.09112363597172*I,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Complex Field with 53 bits of precision
              Defn: a |--> -0.62996052494743... + 1.09112363597172*I,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Complex Field with 53 bits of precision
              Defn: a |--> 1.25992104989487
            ]

        Test that :trac:`15053` is fixed::

            sage: K = NumberField(x^3 - 2, 'a')
            sage: K.embeddings(GF(3))
            []
        """
        try:
            # this should be concordant with automorphisms
            return self.__embeddings[K]
        except AttributeError:
            self.__embeddings = {}
        except KeyError:
            pass
        if K is self:
            return self.automorphisms()
        if K.characteristic() != 0:
            return Sequence([], immutable=True, check=False, universe=self.Hom(K))

        f = K['x'](self.defining_polynomial())
        r = f.roots(); r.sort()
        v = [self.hom([e[0]], check=False) for e in r]
        # If there is an embedding that preserves variable names
        # then it is most natural, so we put it first.
        put_natural_embedding_first(v)

        self.__embeddings[K] = Sequence(v, cr=v!=[], immutable=True,
                                        check=False, universe=self.Hom(K))
        return self.__embeddings[K]

    def Minkowski_embedding(self, B=None, prec=None):
        r"""
        Return an nxn matrix over RDF whose columns are the images of the
        basis `\{1, \alpha, \dots, \alpha^{n-1}\}` of self over
        `\QQ` (as vector spaces), where here
        `\alpha` is the generator of self over
        `\QQ`, i.e. self.gen(0). If B is not None, return
        the images of the vectors in B as the columns instead. If prec is
        not None, use RealField(prec) instead of RDF.

        This embedding is the so-called "Minkowski embedding" of a number
        field in `\RR^n`: given the `n` embeddings
        `\sigma_1, \dots, \sigma_n` of self in
        `\CC`, write `\sigma_1, \dots, \sigma_r`
        for the real embeddings, and
        `\sigma_{r+1}, \dots, \sigma_{r+s}` for choices of one of
        each pair of complex conjugate embeddings (in our case, we simply
        choose the one where the image of `\alpha` has positive
        real part). Here `(r,s)` is the signature of self. Then the
        Minkowski embedding is given by

        .. math::

            x \mapsto ( \sigma_1(x), \dots,
                     \sigma_r(x), \sqrt{2}\Re(\sigma_{r+1}(x)),
                     \sqrt{2}\Im(\sigma_{r+1}(x)), \dots,
                     \sqrt{2}\Re(\sigma_{r+s}(x)),
                     \sqrt{2}\Im(\sigma_{r+s}(x)))

        Equivalently, this is an embedding of self in `\RR^n` so
        that the usual norm on `\RR^n` coincides with
        `|x| = \sum_i |\sigma_i(x)|^2` on self.

        TODO: This could be much improved by implementing homomorphisms
        over VectorSpaces.

        EXAMPLES::

            sage: F.<alpha> = NumberField(x^3+2)
            sage: F.Minkowski_embedding()
            [ 1.00000000000000 -1.25992104989487  1.58740105196820]
            [ 1.41421356237... 0.8908987181... -1.12246204830...]
            [0.000000000000000  1.54308184421...  1.94416129723...]
            sage: F.Minkowski_embedding([1, alpha+2, alpha^2-alpha])
            [ 1.00000000000000 0.740078950105127  2.84732210186307]
            [ 1.41421356237...  3.7193258428... -2.01336076644...]
            [0.000000000000000  1.54308184421... 0.40107945302...]
            sage: F.Minkowski_embedding() * (alpha + 2).vector().column()
            [0.740078950105127]
            [ 3.7193258428...]
            [ 1.54308184421...]
        """
        n = self.degree()
        if prec is None:
            R = sage.rings.real_double.RDF
        else:
            R = sage.rings.real_mpfr.RealField(prec)
        r,s = self.signature()
        places = self.places(prec=prec)

        if B is None:
            B = [ (self.gen(0))**i for i in range(n) ]

        A = ZZ['x']
        f = A.gen(0)**2-2
        sqrt2 = f.roots(R)[1][0]

        d = {}

        for col in range(n):

            for row in range(r):
                d[(row,col)] = places[row](B[col])

            for i in range(s):
                z = places[r+i](B[col])
                d[(r+2*i,col)] = z.real()*sqrt2
                d[(r+2*i+1,col)] = z.imag()*sqrt2


        M = sage.matrix.all.matrix(d)

        return M


    def places(self, all_complex=False, prec=None):
        """
        Return the collection of all infinite places of self.

        By default, this returns the set of real places as
        homomorphisms into RIF first, followed by a choice of one of
        each pair of complex conjugate homomorphisms into CIF.

        On the other hand, if prec is not None, we simply return places
        into RealField(prec) and ComplexField(prec) (or RDF, CDF if
        prec=53). One can also use ``prec=infinity``, which returns embeddings
        into the field `\overline{\QQ}` of algebraic numbers (or its subfield
        `\mathbb{A}` of algebraic reals); this permits exact computatation, but
        can be extremely slow.

        There is an optional flag all_complex, which defaults to False. If
        all_complex is True, then the real embeddings are returned as
        embeddings into CIF instead of RIF.

        EXAMPLES::

            sage: F.<alpha> = NumberField(x^3-100*x+1) ; F.places()
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 - 100*x + 1
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> -10.00499625499181184573367219280,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 - 100*x + 1
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> 0.01000001000003000012000055000273,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 - 100*x + 1
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> 9.994996244991781845613530439509]

        ::

            sage: F.<alpha> = NumberField(x^3+7) ; F.places()
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> -1.912931182772389101199116839549,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 53 bits of precision
            Defn: alpha |--> 0.956465591386195 + 1.65664699997230*I]

        ::

            sage: F.<alpha> = NumberField(x^3+7) ; F.places(all_complex=True)
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 53 bits of precision
            Defn: alpha |--> -1.91293118277239,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 53 bits of precision
            Defn: alpha |--> 0.956465591386195 + 1.65664699997230*I]
            sage: F.places(prec=10)
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Real Field with 10 bits of precision
            Defn: alpha |--> -1.9,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 10 bits of precision
            Defn: alpha |--> 0.96 + 1.7*I]
        """
        if prec is None:
            R = RIF
            C = CIF

        elif prec == 53:
            R = sage.rings.real_double.RDF
            C = sage.rings.complex_double.CDF

        elif prec == sage.rings.infinity.Infinity:
            R = sage.rings.all.AA
            C = sage.rings.all.QQbar

        else:
            R = sage.rings.real_mpfr.RealField(prec)
            C = sage.rings.complex_field.ComplexField(prec)

        ## first, find the intervals with roots, and see how much
        ## precision we need to approximate the roots
        ##
        all_intervals = [ x[0] for x in self.defining_polynomial().roots(C) ]

        ## first, set up the real places
        if all_complex:
            real_intervals = [ x for x in all_intervals if x.imag().is_zero() ]
        else:
            real_intervals = [ x[0] for x in self.defining_polynomial().roots(R) ]

        if prec is None:
            real_places = [ self.hom([i.center()], check=False) for i in real_intervals ]

            complex_places = [ self.hom([i.center()], check=False) for i in
                               all_intervals if i.imag() > 0 ]
        else:
            real_places = [ self.hom([i], check=False) for i in real_intervals ]

            complex_places = [ self.hom([i], check=False) for i in
                               all_intervals if i.imag() > 0 ]

        return real_places + complex_places

    def real_places(self, prec=None):
        """
        Return all real places of self as homomorphisms into RIF.

        EXAMPLES::

            sage: F.<alpha> = NumberField(x^4-7) ; F.real_places()
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^4 - 7
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> -1.626576561697785743211232345494,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^4 - 7
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> 1.626576561697785743211232345494]
        """
        return self.places(prec=prec)[0:self.signature()[0]]

    def relativize(self, alpha, names):
        r"""
        Given an element in self or an embedding of a subfield into self,
        return a relative number field `K` isomorphic to self that is relative
        over the absolute field `\QQ(\alpha)` or the domain of `alpha`, along
        with isomorphisms from `K` to self and from self to `K`.

        INPUT:

        - ``alpha`` - an element of self  or an embedding of a subfield into
          self
        - ``names`` - 2-tuple of names of generator for output field K and the
          subfield QQ(alpha) names[0] generators K and names[1] QQ(alpha).

        OUTPUT:

        K -- relative number field

        Also, ``K.structure()`` returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an isomorphism
        from self to K.

        EXAMPLES::

            sage: K.<a> = NumberField(x^10 - 2)
            sage: L.<c,d> = K.relativize(a^4 + a^2 + 2); L
            Number Field in c with defining polynomial x^2 - 1/5*d^4 + 8/5*d^3 - 23/5*d^2 + 7*d - 18/5 over its base field
            sage: c.absolute_minpoly()
            x^10 - 2
            sage: d.absolute_minpoly()
            x^5 - 10*x^4 + 40*x^3 - 90*x^2 + 110*x - 58
            sage: (a^4 + a^2 + 2).minpoly()
            x^5 - 10*x^4 + 40*x^3 - 90*x^2 + 110*x - 58
            sage: from_L, to_L = L.structure()
            sage: to_L(a)
            c
            sage: to_L(a^4 + a^2 + 2)
            d
            sage: from_L(to_L(a^4 + a^2 + 2))
            a^4 + a^2 + 2

        The following demonstrates distinct embeddings of a subfield into a
        larger field::

            sage: K.<a> = NumberField(x^4 + 2*x^2 + 2)
            sage: K0 = K.subfields(2)[0][0]; K0
            Number Field in a0 with defining polynomial x^2 - 2*x + 2
            sage: rho, tau = K0.embeddings(K)
            sage: L0 = K.relativize(rho(K0.gen()), 'b'); L0
            Number Field in b0 with defining polynomial x^2 - b1 + 2 over its base field
            sage: L1 = K.relativize(rho, 'b'); L1
            Number Field in b0 with defining polynomial x^2 - a0 + 2 over its base field
            sage: L2 = K.relativize(tau, 'b'); L2
            Number Field in b0 with defining polynomial x^2 + a0 over its base field
            sage: L0.base_field() is K0
            False
            sage: L1.base_field() is K0
            True
            sage: L2.base_field() is K0
            True

        Here we see that with the different embeddings, the relative norms are
        different::

            sage: a0 = K0.gen()
            sage: L1_into_K, K_into_L1 = L1.structure()
            sage: L2_into_K, K_into_L2 = L2.structure()
            sage: len(K.factor(41))
            4
            sage: w1 = -a^2 + a + 1; P = K.ideal([w1])
            sage: Pp = L1.ideal(K_into_L1(w1)).ideal_below(); Pp == K0.ideal([4*a0 + 1])
            True
            sage: Pp == w1.norm(rho)
            True

            sage: w2 = a^2 + a - 1; Q = K.ideal([w2])
            sage: Qq = L2.ideal(K_into_L2(w2)).ideal_below(); Qq == K0.ideal([-4*a0 + 9])
            True
            sage: Qq == w2.norm(tau)
            True

            sage: Pp == Qq
            False

        TESTS:

        We can relativize over the whole field::

            sage: K.<a> = NumberField(x^4 + 2*x^2 + 2)
            sage: K.relativize(K.gen(), 'a')
            Number Field in a0 with defining polynomial x - a1 over its base field
            sage: K.relativize(2*K.gen(), 'a')
            Number Field in a0 with defining polynomial x - 1/2*a1 over its base field

        We can relativize over the prime field::

            sage: L = K.relativize(K(1), 'a'); L
            Number Field in a0 with defining polynomial x^4 + 2*x^2 + 2 over its base field
            sage: L.base_field()
            Number Field in a1 with defining polynomial x - 1
            sage: L.base_field().base_field()
            Rational Field

            sage: L = K.relativize(K(2), 'a'); L
            Number Field in a0 with defining polynomial x^4 + 2*x^2 + 2 over its base field
            sage: L.base_field()
            Number Field in a1 with defining polynomial x - 2
            sage: L.base_field().base_field()
            Rational Field

            sage: L = K.relativize(K(0), 'a'); L
            Number Field in a0 with defining polynomial x^4 + 2*x^2 + 2 over its base field
            sage: L.base_field()
            Number Field in a1 with defining polynomial x
            sage: L.base_field().base_field()
            Rational Field

        We can relativize over morphisms returned by self.subfields()::

            sage: L = NumberField(x^4 + 1, 'a')
            sage: [L.relativize(h, 'c') for (f,h,i) in L.subfields()]
            [Number Field in c0 with defining polynomial x^4 + 1 over its base field, Number Field in c0 with defining polynomial x^2 - 1/2*a1 over its base field, Number Field in c0 with defining polynomial x^2 - a2*x - 1 over its base field, Number Field in c0 with defining polynomial x^2 - a3*x + 1 over its base field, Number Field in c0 with defining polynomial x - a4 over its base field]

        We can relativize over a relative field::

            sage: K.<z> = CyclotomicField(16)
            sage: L, L_into_K, _ = K.subfields(4)[0]; L
            Number Field in z0 with defining polynomial x^4 + 16
            sage: F, F_into_L, _ = L.subfields(2)[0]; F
            Number Field in z00 with defining polynomial x^2 + 64

            sage: L_over_F = L.relativize(F_into_L, 'c'); L_over_F
            Number Field in c0 with defining polynomial x^2 - 1/2*z00 over its base field
            sage: L_over_F_into_L, _ = L_over_F.structure()

            sage: K_over_rel = K.relativize(L_into_K * L_over_F_into_L, 'a'); K_over_rel
            Number Field in a0 with defining polynomial x^2 - 1/2*c0 over its base field
            sage: K_over_rel.base_field() is L_over_F
            True
            sage: K_over_rel.structure()
            (Relative number field morphism:
              From: Number Field in a0 with defining polynomial x^2 - 1/2*c0 over its base field
              To:   Cyclotomic Field of order 16 and degree 8
              Defn: a0 |--> z
                    c0 |--> 2*z^2
                    z00 |--> 8*z^4, Ring morphism:
              From: Cyclotomic Field of order 16 and degree 8
              To:   Number Field in a0 with defining polynomial x^2 - 1/2*c0 over its base field
              Defn: z |--> a0)

        We can relativize over a really large field::

            sage: K.<a> = CyclotomicField(3^3*2^3)
            sage: R = K.relativize(a^(3^2), 't'); R
            Number Field in t0 with defining polynomial x^9 - t1 over its base field
            sage: R.structure()
            (Relative number field morphism:
              From: Number Field in t0 with defining polynomial x^9 - t1 over its base field
              To:   Cyclotomic Field of order 216 and degree 72
              Defn: t0 |--> a
                    t1 |--> a^9,
             Ring morphism:
              From: Cyclotomic Field of order 216 and degree 72
              To:   Number Field in t0 with defining polynomial x^9 - t1 over its base field
              Defn: a |--> t0)
        """
        # step 1: construct the abstract field generated by alpha.w
        # step 2: make a relative extension of it.
        # step 3: construct isomorphisms
        from sage.all import vector, matrix

        names = sage.structure.parent_gens.normalize_names(2, names)

        from sage.categories.map import is_Map
        if is_Map(alpha):
            # alpha better be a morphism with codomain self
            if alpha.codomain() != self:
                raise ValueError, "Co-domain of morphism must be self"
            L = alpha.domain()
            alpha = alpha(L.gen()) # relativize over phi's domain
            f = L.defining_polynomial() # = alpha.minpoly()
        else:
            # alpha must be an element coercible to self
            alpha = self(alpha)
            f = alpha.minpoly()
            L = NumberField(f, names[1])

        # now we do some linear algebra to find the minpoly of self.gen() over L
        L_into_self = L.hom([alpha])

        extdeg = self.absolute_degree() // L.absolute_degree() # [ L : self ]
        a = self.gen()

        # we will find a linear relation between small powers of a over L
        basis = [ a**i * b for i in range(extdeg) for b in map(L_into_self, L.power_basis()) ]
        basis.append(a**extdeg) # this one makes the basis no longer a basis
        mat = matrix([ b.vector() for b in basis ])
        soln_space = mat.left_kernel(mat.row_space()(0))
        # the solution space is one dimensional and the last entry is non-zero
        # because a satisfies no smaller linear relation
        assert soln_space.dimension() == 1
        (reln, ) = soln_space.basis()
        assert reln[-1] != 0
        reln = reln * ~reln[-1]

        # now we need to get those coeffs in L
        coeff_mat = matrix(extdeg, f.degree(), list(reln)[:-1]) # easy way to divide into the correct lengths
        coeffs_in_L = [ r*vector(L.power_basis()) for r in coeff_mat.rows() ]
        # f is the minimal polynomial of a over L
        f = L['x'](coeffs_in_L + [1])
        # sanity check...
        mp_in_self = self['x'](map(L_into_self, f.coeffs()))
        assert mp_in_self(a) == 0

        M = L.extension(f, names[0])
        beta = M(L.gen())
        to_M = self.hom([M.gen(0)], M, check=True)  # be paranoid
        assert to_M(alpha) == beta
        # Bingo.
        # We have now constructed a relative
        # number field M, and an isomorphism
        # self --> M that sends alpha to
        # the generator of the intermediate field.

        # this strange step avoids computing an embedding of the base_field L
        # into self.
        base_hom = L.hom([alpha], self)
        from_M = M.Hom(self)([self.gen()], base_hom=base_hom, check=True)

        # at some point in the future, we will worry about caching since relative
        # number fields should be cached but aren't -- yet
        M._set_structure(from_M, to_M)
        return M

    # Synonyms so that terminology appropriate to relative number fields
    # can be applied to an absolute number field:

    def absolute_degree(self):
        """
        A synonym for degree.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_degree()
            2
        """
        return self.degree()

    def relative_degree(self):
        """
        A synonym for degree.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_degree()
            2
        """
        return self.degree()

    def absolute_polynomial(self):
        """
        A synonym for polynomial.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_polynomial()
            x^2 + 1
        """
        return self.polynomial()

    def relative_polynomial(self):
        """
        A synonym for polynomial.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_polynomial()
            x^2 + 1
        """
        return self.polynomial()

    def absolute_vector_space(self):
        """
        A synonym for vector_space.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_vector_space()
            (Vector space of dimension 2 over Rational Field,
             Isomorphism map:
              From: Vector space of dimension 2 over Rational Field
              To:   Number Field in i with defining polynomial x^2 + 1,
             Isomorphism map:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Vector space of dimension 2 over Rational Field)
        """
        return self.vector_space()

    def relative_vector_space(self):
        """
        A synonym for vector_space.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_vector_space()
            (Vector space of dimension 2 over Rational Field,
             Isomorphism map:
              From: Vector space of dimension 2 over Rational Field
              To:   Number Field in i with defining polynomial x^2 + 1,
             Isomorphism map:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Vector space of dimension 2 over Rational Field)
        """
        return self.vector_space()

    def absolute_discriminant(self):
        """
        A synonym for discriminant.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_discriminant()
            -4
        """
        return self.discriminant()

    def relative_discriminant(self):
        """
        A synonym for discriminant.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_discriminant()
            -4
        """
        return self.discriminant()

    def absolute_different(self):
        """
        A synonym for different.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_different()
            Fractional ideal (2)
        """
        return self.different()

    def relative_different(self):
        """
        A synonym for different.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_different()
            Fractional ideal (2)
        """
        return self.different()

    def hilbert_symbol(self, a, b, P = None):
        r"""
        Returns the Hilbert symbol `(a,b)_P` for a prime P of self
        and non-zero elements a and b of self.
        If P is omitted, return the global Hilbert symbol `(a,b)` instead.

        INPUT:

        - ``a``, ``b`` -- elements of self

        - ``P`` -- (default: None) If `P` is ``None``, compute the global
          symbol.  Otherwise, `P` should be either a prime ideal of self
          (which may also be given as a generator or set of generators)
          or a real or complex embedding.

        OUTPUT:

        If a or b is zero, returns 0.

        If a and b are non-zero and P is specified, returns
        the Hilbert symbol `(a,b)_P`, which is 1 if the equation
        `a x^2 + b y^2 = 1` has a solution in the completion of
        self at P, and is -1 otherwise.

        If a and b are non-zero and P is unspecified, returns 1
        if the equation has a solution in self and -1 otherwise.

        EXAMPLES:

        Some global examples::

            sage: K.<a> = NumberField(x^2 - 23)
            sage: K.hilbert_symbol(0, a+5)
            0
            sage: K.hilbert_symbol(a, 0)
            0
            sage: K.hilbert_symbol(-a, a+1)
            1
            sage: K.hilbert_symbol(-a, a+2)
            -1
            sage: K.hilbert_symbol(a, a+5)
            -1

        That the latter two are unsolvable should be visible in local
        obstructions.  For the first, this is a prime ideal above 19.
        For the second, the ramified prime above 23::

            sage: K.hilbert_symbol(-a, a+2, a+2)
            -1
            sage: K.hilbert_symbol(a, a+5, K.ideal(23).factor()[0][0])
            -1

        More local examples::

            sage: K.hilbert_symbol(a, 0, K.ideal(5))
            0
            sage: K.hilbert_symbol(a, a+5, K.ideal(5))
            1
            sage: K.hilbert_symbol(a+1, 13, (a+6)*K.maximal_order())
            -1
            sage: [emb1, emb2] = K.embeddings(AA)
            sage: K.hilbert_symbol(a, -1, emb1)
            -1
            sage: K.hilbert_symbol(a, -1, emb2)
            1

        Ideals P can be given by generators::

            sage: K.<a> = NumberField(x^5 - 23)
            sage: pi = 2*a^4 + 3*a^3 + 4*a^2 + 15*a + 11
            sage: K.hilbert_symbol(a, a+5, pi)
            1
            sage: rho = 2*a^4 + 3*a^3 + 4*a^2 + 15*a + 11
            sage: K.hilbert_symbol(a, a+5, rho)
            1

        This also works for non-principal ideals::

            sage: K.<a> = QuadraticField(-5)
            sage: P = K.ideal(3).factor()[0][0]
            sage: P.gens_reduced()  # random, could be the other factor
            (3, a + 1)
            sage: K.hilbert_symbol(a, a+3, P)
            1
            sage: K.hilbert_symbol(a, a+3, [3, a+1])
            1

        Primes above 2::

            sage: K.<a> = NumberField(x^5 - 23)
            sage: O = K.maximal_order()
            sage: p = [p[0] for p in (2*O).factor() if p[0].norm() == 16][0]
            sage: K.hilbert_symbol(a, a+5, p)
            1
            sage: K.hilbert_symbol(a, 2, p)
            1
            sage: K.hilbert_symbol(-1, a-2, p)
            -1

        Various real fields are allowed::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: K.hilbert_symbol(a/3, 1/2, K.embeddings(RDF)[0])
            1
            sage: K.hilbert_symbol(a/5, -1, K.embeddings(RR)[0])
            -1
            sage: [K.hilbert_symbol(a, -1, e) for e in K.embeddings(AA)]
            [-1]

        Real embeddings are not allowed to be disguised as complex embeddings::

            sage: K.<a> = QuadraticField(5)
            sage: K.hilbert_symbol(-1, -1, K.embeddings(CC)[0])
            Traceback (most recent call last):
            ...
            ValueError: Possibly real place (=Ring morphism:
              From: Number Field in a with defining polynomial x^2 - 5
              To:   Complex Field with 53 bits of precision
              Defn: a |--> -2.23606797749979) given as complex embedding in hilbert_symbol. Is it real or complex?
            sage: K.hilbert_symbol(-1, -1, K.embeddings(QQbar)[0])
            Traceback (most recent call last):
            ...
            ValueError: Possibly real place (=Ring morphism:
              From: Number Field in a with defining polynomial x^2 - 5
              To:   Algebraic Field
              Defn: a |--> -2.236067977499790?) given as complex embedding in hilbert_symbol. Is it real or complex?
            sage: K.<b> = QuadraticField(-5)
            sage: K.hilbert_symbol(-1, -1, K.embeddings(CDF)[0])
            1
            sage: K.hilbert_symbol(-1, -1, K.embeddings(QQbar)[0])
            1

        a and b do not have to be integral or coprime::

            sage: K.<i> = QuadraticField(-1)
            sage: O = K.maximal_order()
            sage: K.hilbert_symbol(1/2, 1/6, 3*O)
            1
            sage: p = 1+i
            sage: K.hilbert_symbol(p, p, p)
            1
            sage: K.hilbert_symbol(p, 3*p, p)
            -1
            sage: K.hilbert_symbol(3, p, p)
            -1
            sage: K.hilbert_symbol(1/3, 1/5, 1+i)
            1
            sage: L = QuadraticField(5, 'a')
            sage: L.hilbert_symbol(-3, -1/2, 2)
            1

        Various other examples::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: K.hilbert_symbol(-6912, 24, -a^2-a-2)
            1
            sage: K.<a> = NumberField(x^5-23)
            sage: P = K.ideal(-1105*a^4 + 1541*a^3 - 795*a^2 - 2993*a + 11853)
            sage: Q = K.ideal(-7*a^4 + 13*a^3 - 13*a^2 - 2*a + 50)
            sage: b = -a+5
            sage: K.hilbert_symbol(a,b,P)
            1
            sage: K.hilbert_symbol(a,b,Q)
            1
            sage: K.<a> = NumberField(x^5-23)
            sage: P = K.ideal(-1105*a^4 + 1541*a^3 - 795*a^2 - 2993*a + 11853)
            sage: K.hilbert_symbol(a, a+5, P)
            1
            sage: K.hilbert_symbol(a, 2, P)
            1
            sage: K.hilbert_symbol(a+5, 2, P)
            -1

        AUTHOR:

        - Aly Deines (2010-08-19): part of the doctests

        - Marco Streng (2010-12-06)
        """
        if a.is_zero() or b.is_zero():
            return 0
        a = self(a)
        b = self(b)
        if P is None:
            return pari(self).nfhilbert(a, b)

        from sage.rings.morphism import is_RingHomomorphism
        if is_RingHomomorphism(P):
            if not P.domain() == self:
                raise ValueError, "Domain of P (=%s) should be self (=%s) in self.hilbert_symbol" % (P, self)
            codom = P.codomain()
            from sage.rings.complex_field import is_ComplexField
            from sage.rings.complex_interval_field import is_ComplexIntervalField
            from sage.rings.real_mpfr import is_RealField
            from sage.rings.real_mpfi import is_RealIntervalField
            from sage.rings.all import (AA, CDF, QQbar, RDF)
            if is_ComplexField(codom) or is_ComplexIntervalField(codom) or \
                                         codom is CDF or codom is QQbar:
                if P(self.gen()).imag() == 0:
                    raise ValueError, "Possibly real place (=%s) given as complex embedding in hilbert_symbol. Is it real or complex?" % P
                return 1
            if is_RealField(codom) or codom is RDF or codom is AA:
                if P(a) > 0 or P(b) > 0:
                    return 1
                return -1
        if not is_NumberFieldIdeal(P):
            P = self.ideal(P)
        if not P.number_field() == self:
            raise ValueError, "P (=%s) should be an ideal of self (=%s) in hilbert_symbol, not of %s" % (P, self, P.number_field())
        if not P.is_prime():
            raise ValueError, "Non-prime ideal P (=%s) in hilbert_symbol" % P
        return pari(self).nfhilbert(a, b, P.pari_prime())

    def hilbert_conductor(self,a,b):
        """
        This is the product of all (finite) primes where the Hilbert symbol is -1.
        What is the same, this is the (reduced) discriminant of the quaternion
        algebra `(a,b)` over a number field.

        INPUT:

        -``a``,``b`` -- elements of the number field `self`

        OUTPUT:

        - squarefree ideal of the ring of integers of `self`

        EXAMPLES::

            sage: F.<a> = NumberField(x^2-x-1)
            sage: F.hilbert_conductor(2*a,F(-1))
            Fractional ideal (2)
            sage: K.<b> = NumberField(x^3-4*x+2)
            sage: K.hilbert_conductor(K(2),K(-2))
            Fractional ideal (1)
            sage: K.hilbert_conductor(K(2*b),K(-2))
            Fractional ideal (-b^2 - b + 2)


        AUTHOR:

        - Aly Deines

        """
        a, b = self(a), self(b)
        d = self.ideal(1)
        for p in union(union( self.ideal(2).prime_factors(), self.ideal(a).prime_factors()), self.ideal(b).prime_factors()):
            if self.hilbert_symbol(a,b,p) == -1:
                d *= p
        return d


class NumberField_cyclotomic(NumberField_absolute):
    """
    Create a cyclotomic extension of the rational field.

    The command CyclotomicField(n) creates the n-th cyclotomic field,
    obtained by adjoining an n-th root of unity to the rational field.

    EXAMPLES::

        sage: CyclotomicField(3)
        Cyclotomic Field of order 3 and degree 2
        sage: CyclotomicField(18)
        Cyclotomic Field of order 18 and degree 6
        sage: z = CyclotomicField(6).gen(); z
        zeta6
        sage: z^3
        -1
        sage: (1+z)^3
        6*zeta6 - 3

    ::

        sage: K = CyclotomicField(197)
        sage: loads(K.dumps()) == K
        True
        sage: loads((z^2).dumps()) == z^2
        True

    ::

        sage: cf12 = CyclotomicField( 12 )
        sage: z12 = cf12.0
        sage: cf6 = CyclotomicField( 6 )
        sage: z6 = cf6.0
        sage: FF = Frac( cf12['x'] )
        sage: x = FF.0
        sage: z6*x^3/(z6 + x)
        zeta12^2*x^3/(x + zeta12^2)

    ::

        sage: cf6 = CyclotomicField(6) ; z6 = cf6.gen(0)
        sage: cf3 = CyclotomicField(3) ; z3 = cf3.gen(0)
        sage: cf3(z6)
        zeta3 + 1
        sage: cf6(z3)
        zeta6 - 1
        sage: type(cf6(z3))
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
        sage: cf1 = CyclotomicField(1) ; z1 = cf1.0
        sage: cf3(z1)
        1
        sage: type(cf3(z1))
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
    """
    def __init__(self, n, names, embedding=None, assume_disc_small=False, maximize_at_primes=None):
        """
        A cyclotomic field, i.e., a field obtained by adjoining an n-th
        root of unity to the rational numbers.

        EXAMPLES::

            sage: k = CyclotomicField(3)
            sage: type(k)
            <class 'sage.rings.number_field.number_field.NumberField_cyclotomic_with_category'>

        TESTS::

            sage: TestSuite(k).run()
            sage: type(CyclotomicField(4).zero_element())
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: type(CyclotomicField(6).one_element())
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: type(CyclotomicField(15).zero_element())
            <type 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>
        """
        f = QQ['x'].cyclotomic_polynomial(n)
        if names[0].startswith('zeta'):
            latex_name = "\\zeta_{%s}"%n
        else:
            latex_name = None
        self.__n = n = integer.Integer(n)
        NumberField_absolute.__init__(self, f,
                                      name= names,
                                      latex_name=latex_name,
                                      check=False,
                                      embedding = embedding,
                                      assume_disc_small=assume_disc_small,
                                      maximize_at_primes=maximize_at_primes)
        if n%2:
            self.__zeta_order = 2*n
        else:
            self.__zeta_order = n
        ## quadratic number fields require this:
        if f.degree() == 2:
            # define a boolean flag as for NumberField_quadratic to know, which
            # square root we choose (True means no embedding or positive
            # imaginary value).
            # Note that the test is done with NumberFieldElement and not with
            # NumberFieldElement_quadratic which requires somehow this flag.
            self._standard_embedding = not CDF.has_coerce_map_from(self) or CDF(self.gen()).imag() > 0

            self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic
            if n == 4:
                self._D = ZZ(-1)
                self._NumberField_generic__gen = self._element_class(self, (QQ(0), QQ(1)))
            else:
                ## n is 3 or 6
                self._D = ZZ(-3)
                one_half = ZZ(1)/ZZ(2)
                if n == 3:
                    self._NumberField_generic__gen = self._element_class(self, (one_half-1, one_half))
                else:
                    self._NumberField_generic__gen = self._element_class(self, (one_half, one_half))

            # NumberField_absolute.__init__(...) set _zero_element and
            # _one_element to NumberFieldElement_absolute values, which is
            # wrong (and dangerous; such elements can actually be used to
            # crash Sage: see #5316).  Overwrite them with correct values.
            self._zero_element = self._element_class(self, (QQ(0),QQ(0)))
            self._one_element =  self._element_class(self, (QQ(1),QQ(0)))

        zeta = self.gen()
        zeta._set_multiplicative_order(n)

    def __reduce__(self):
        """
        TESTS::

            sage: K.<zeta7> = CyclotomicField(7)
            sage: L = loads(dumps(K))
            sage: L
            Cyclotomic Field of order 7 and degree 6
            sage: L == K
            True
        """
        return NumberField_cyclotomic_v1, (self.__n, self.variable_name(), self.gen_embedding())

    def construction(self):
        F,R = NumberField_generic.construction(self)
        F.cyclotomic = self.__n
        return F,R

    def _magma_init_(self, magma):
        """
        Function returning a string to create this cyclotomic field in
        Magma.

        .. note::

           The Magma generator name is also initialized to be the same
           as for the Sage field.

        EXAMPLES::

            sage: K=CyclotomicField(7,'z')
            sage: K._magma_init_(magma)                                # optional - magma
            'SageCreateWithNames(CyclotomicField(7),["z"])'
            sage: K=CyclotomicField(7,'zeta')
            sage: K._magma_init_(magma)                                # optional - magma
            'SageCreateWithNames(CyclotomicField(7),["zeta"])'
        """
        s = 'CyclotomicField(%s)'%self.__n
        return magma._with_names(s, self.variable_names())

    def _gap_init_(self):
        """
        Return a string that provides a representation of ``self`` in GAP.

        TESTS:

            sage: K = CyclotomicField(8)
            sage: gap(K)   # indirect doctest
            CF(8)
            sage: gap(K.0)
            E(8)
            sage: K(gap(K.0^5)); K(gap(K.0^5))==K.0^5
            -zeta8
            True

        The following was the motivating example to introduce
        a genuine representation of cyclotomic fields in the
        GAP interface -- see ticket #5618.
        ::

            sage: H = AlternatingGroup(4)
            sage: g = H.list()[1]
            sage: K = H.subgroup([g])
            sage: z = CyclotomicField(3).an_element(); z
            zeta3
            sage: c = K.character([1,z,z**2]); c
            Character of Subgroup of (Alternating group of order 4!/2 as a permutation group) generated by [(2,3,4)]
            sage: c(g^2); z^2
            -zeta3 - 1
            -zeta3 - 1

        """
        return 'CyclotomicField(%s)'%self.__n

    def _libgap_(self):
        """
        Return a LibGAP representation of ``self``.

        TESTS:

            sage: K = CyclotomicField(8)
            sage: K._libgap_()
            CF(8)
            sage: libgap(K)   # indirect doctest
            CF(8)
        """
        from sage.libs.gap.libgap import libgap
        return libgap.CyclotomicField(self.__n)

    def _repr_(self):
        r"""
        Return string representation of this cyclotomic field.

        The "order" of the cyclotomic field `\QQ(\zeta_n)`
        in the string output refers to the order of the `\zeta_n`,
        i.e., it is the integer `n`. The degree is the degree of
        the field as an extension of `\QQ`.

        EXAMPLES::

            sage: CyclotomicField(4)._repr_()
            'Cyclotomic Field of order 4 and degree 2'
            sage: CyclotomicField(400)._repr_()
            'Cyclotomic Field of order 400 and degree 160'
        """
        return "Cyclotomic Field of order %s and degree %s"%(
                self.__n, self.degree())

    def _n(self):
        """
        Return the n used to create this cyclotomic field.

        EXAMPLES::

            sage: CyclotomicField(3).zeta_order()
            6
            sage: CyclotomicField(3)._n()
            3
        """
        return self.__n

    def _latex_(self):
        """
        Return the latex representation of this cyclotomic field.

        EXAMPLES::

            sage: Z = CyclotomicField(4)
            sage: Z.gen()
            zeta4
            sage: latex(Z) # indirect doctest
            \Bold{Q}(\zeta_{4})

        Latex printing respects the generator name::

            sage: k.<a> = CyclotomicField(4)
            sage: latex(k)
            \Bold{Q}[a]/(a^{2} + 1)
            sage: k
            Cyclotomic Field of order 4 and degree 2
            sage: k.gen()
            a

        TESTS:

        We check that the bug reported on :trac:`8938` is fixed::

            sage: C5.<z> = CyclotomicField(5)
            sage: P.<s, t> = C5[]
            sage: f = (z^2 + z)*s
            sage: f
            (z^2 + z)*s
            sage: latex(f)
            \left(z^{2} + z\right) s
        """
        v = self.latex_variable_name()
        if v.startswith('\\zeta_'):
            return "%s(%s)"%(latex(QQ), v)
        else:
            return NumberField_generic._latex_(self)

    def _coerce_map_from_(self, K):
        r"""
        The cyclotomic field `\Q(\zeta_n)` coerces into the cyclotomic field
        `\Q(\zeta_m)` iff `n'|m`, where `n'` is the odd part of `n` if `4 \not
        | n` and `n'=n` otherwise.

        The morphism is consistant with the chosen embedding into `\CC`.

        If `K` is not a cyclotomic field, the normal coercion rules for number
        fields are used.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(12)
            sage: L.<b> = CyclotomicField(132)
            sage: L.coerce_map_from(K) # indirect doctest
            Generic morphism:
              From: Cyclotomic Field of order 12 and degree 4
              To:   Cyclotomic Field of order 132 and degree 40
              Defn: a -> b^11
            sage: a + b
            b^11 + b
            sage: L.coerce_map_from(CyclotomicField(4, 'z'))
            Generic morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Cyclotomic Field of order 132 and degree 40
              Defn: z -> b^33
            sage: L.coerce_map_from(CyclotomicField(5, 'z')) is None
            True

            sage: K.<a> = CyclotomicField(3)
            sage: L.<b> = CyclotomicField(6)
            sage: L.coerce_map_from(K)
            Generic morphism:
              From: Cyclotomic Field of order 3 and degree 2
              To:   Cyclotomic Field of order 6 and degree 2
              Defn: a -> b - 1
            sage: K.coerce_map_from(L)
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 3 and degree 2
              Defn: b -> a + 1

            sage: CyclotomicField(33).coerce_map_from(CyclotomicField(66))
            Generic morphism:
              From: Cyclotomic Field of order 66 and degree 20
              To:   Cyclotomic Field of order 33 and degree 20
              Defn: zeta66 -> -zeta33^17
            sage: CyclotomicField(15).coerce_map_from(CyclotomicField(6))
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 15 and degree 8
              Defn: zeta6 -> zeta15^5 + 1

        Check that #12632 is fixed::

            sage: K1 = CyclotomicField(1); K2 = CyclotomicField(2)
            sage: K1.coerce_map_from(K2)
            Generic morphism:
              From: Cyclotomic Field of order 2 and degree 1
              To:   Cyclotomic Field of order 1 and degree 1
              Defn: zeta2 -> -1

        Check that custom embeddings are respected (:trac:`13765`)::

            sage: z105 = CDF(exp(2*pi*I/105))
            sage: Ka.<a> = CyclotomicField(105, embedding=z105^11)
            sage: Kb.<b> = CyclotomicField(35, embedding=z105^6)
            sage: Ka.coerce_map_from(Kb)
            Generic morphism:
              From: Cyclotomic Field of order 35 and degree 24
              To:   Cyclotomic Field of order 105 and degree 48
              Defn: b -> -a^44 - a^42 + a^39 + a^37 + a^35 - a^29 - a^27 - a^25 + a^24 - a^23 + a^22 - a^21 + a^20 + a^18 + a^16 - a^12 - a^10 - a^8 - a^6 + a^5 + a^3 + a
            sage: CC(b)
            0.936234870639737 + 0.351374824081343*I
            sage: CC(-a^44 - a^42 + a^39 + a^37 + a^35 - a^29 - a^27 - a^25 + a^24 - a^23 + a^22 - a^21 + a^20 + a^18 + a^16 - a^12 - a^10 - a^8 - a^6 + a^5 + a^3 + a)
            0.936234870639731 + 0.351374824081341*I

            sage: z15 = CDF(exp(2*pi*I/15))
            sage: CyclotomicField(15).coerce_map_from(CyclotomicField(6, embedding=-z15^5))
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 15 and degree 8
              Defn: zeta6 -> -zeta15^5

            sage: CyclotomicField(15, embedding=z15^4).coerce_map_from(CyclotomicField(6, embedding=-z15^5))
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 15 and degree 8
              Defn: zeta6 -> -zeta15^5
        """
        if isinstance(K, NumberField_cyclotomic):
            if (self.coerce_embedding() is None or K.coerce_embedding() is None):
                return None
            ambient_field = self.coerce_embedding().codomain()
            if not ambient_field.has_coerce_map_from(K.coerce_embedding().codomain()):
                return None
            Kn = K.__n
            n = self.__n
            if Kn.divides(n):
                return number_field_morphisms.CyclotomicFieldEmbedding(K, self)
            if Kn == 2 and n == 1:
                # see #12632
                return number_field_morphisms.NumberFieldEmbedding(K, self, -self.gen())
            if Kn % 4 == 2 and (Kn//2).divides(n):
                e = self._log_gen(ambient_field(-K.gen()))
                return number_field_morphisms.NumberFieldEmbedding(K, self, -self.gen() ** e)
            else:
                return None

        elif self.degree() == 2:
            if K is ZZ:
                return number_field_element_quadratic.Z_to_quadratic_field_element(self)
            if K is QQ:
                return number_field_element_quadratic.Q_to_quadratic_field_element(self)

        return NumberField_absolute._coerce_map_from_(self, K)

    def _log_gen(self, x):
        """
        Returns an integer `e` such that `self.gen()^e == x`, or `None`
        if no such integer exists. This is primarily used to construct
        embedding-respecting coercions.

        If `x` is complex, the result is either an integer `e` such
        that the absolute value of `self.gen()^e-x` is small or
        `None` if no such `e` is found.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(5)
            sage: K._log_gen(CDF(a))
            1
            sage: K._log_gen(CDF(a^4))
            4

            sage: zeta105 = CC(exp(2*pi*i/105))
            sage: K.<a> = CyclotomicField(105, embedding=zeta105^13)
            sage: zeta105^13, CC(a)
            (0.712376096951345 + 0.701797902883992*I, 0.712376096951345 + 0.701797902883991*I)
            sage: K._log_gen(zeta105^26)
            2
            sage: K._log_gen(zeta105)
            97
            sage: zeta105, CC(a^97)
            (0.998210129767735 + 0.0598041539450342*I, 0.998210129767736 + 0.0598041539450313*I)
            sage: K._log_gen(zeta105^3)
            81
            sage: zeta105^3, CC(a)^81
            (0.983929588598630 + 0.178556894798637*I, 0.983929588598631 + 0.178556894798635*I)

            sage: K.<a> = CyclotomicField(5, embedding=None)
            sage: K._log_gen(CDF(.5, -.8)) is None
            True

            sage: zeta5 = cyclotomic_polynomial(5).change_ring(Qp(11)).roots()[0][0]
            sage: zeta5 ^ 5
            1 + O(11^20)
            sage: K.<a> = CyclotomicField(5, embedding=zeta5^2)
            sage: K._log_gen(zeta5)
            3
        """
        if not x.parent().has_coerce_map_from(self):
            return None
        if CDF.has_coerce_map_from(x.parent()):
            x = CDF(x)
        gen = x.parent().coerce(self.gen())
        n = self._n()
        two_pi = 2*RDF.pi()
        if x.parent() is CDF:
            # Let zeta = e^(2*pi*i/n)
            a = (n * x.arg() / two_pi).round()         # x = zeta^a
            b = (n * gen.arg() / two_pi).round()      # gen = zeta^b
            e = mod(a/b, n).lift()          # e is the expected result
            if abs(gen**e-x) < 1/n:        # a sanity check
                return e
        else:
            gen_pow_e = 1
            for e in range(n):
                if gen_pow_e == x:
                    return e
                gen_pow_e *= gen

    def _element_constructor_(self, x):
        """
        Create an element of this cyclotomic field from `x`.

        EXAMPLES:

        The following example illustrates coercion from the
        cyclotomic field Q(zeta_42) to the cyclotomic field Q(zeta_6), in
        a case where such coercion is defined::

            sage: k42 = CyclotomicField(42)
            sage: k6 = CyclotomicField(6)
            sage: a = k42.gen(0)
            sage: b = a^7
            sage: b
            zeta42^7
            sage: k6(b) # indirect doctest
            zeta6
            sage: b^2
            zeta42^7 - 1
            sage: k6(b^2)
            zeta6 - 1

        Conversion of elements of the :class:`~sage.rings.universal_cyclotomic_field.universal_cyclotomic_field.UniversalCyclotomicField`::

            sage: CF = CyclotomicField(5)
            sage: UCF.<E> = UniversalCyclotomicField()
            sage: CF(E(5))
            zeta5

            sage: CF = CyclotomicField(10)
            sage: CF(E(5))
            zeta10^2

       Coercion of GAP cyclotomic elements is also supported::

            sage: CyclotomicField(18)(gap('E(3)')) # indirect doctest
            zeta18^3 - 1

        Converting from rings of integers::

            sage: K.<z> = CyclotomicField(7)
            sage: O = K.maximal_order()
            sage: K(O.1)
            z
            sage: K(O.1^2 + O.1 - 2)
            z^2 + z - 2
        """
        from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import UniversalCyclotomicField

        if isinstance(x, number_field_element.NumberFieldElement):
            if isinstance(x.parent(), NumberField_cyclotomic):
                return self._coerce_from_other_cyclotomic_field(x)
            else:
                return NumberField_absolute._element_constructor_(self, x)
        elif sage.interfaces.gap.is_GapElement(x):
            return self._coerce_from_gap(x)
        elif isinstance(x,UniversalCyclotomicField.Element):
            return self._coerce_from_universal_cyclotomic_field(x)
        elif isinstance(x,str):
            return self._coerce_from_str(x)
        else:
            return self._coerce_non_number_field_element_in(x)

    # TODO:
    # The following is very nice and much more flexible / powerful.
    # However, it is simply not *consistent*, since it totally
    # breaks the doctests in eisenstein_submodule.py.
    # FIX THIS.

##     def _will_be_better_coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
##         """
##         Coerce an element x of a cyclotomic field into self, if at all possible.

##         INPUT:
##             x -- number field element

##             only_canonical -- bool (default: False); Attempt to work,
##                    even in some cases when x is not in a subfield of
##                    the cyclotomics (as long as x is a root of unity).

##         EXAMPLES:
##             sage: k5 = CyclotomicField(5)
##             sage: k3 = CyclotomicField(3)
##             sage: k15 = CyclotomicField(15)
##             sage: k15._coerce_from_other_cyclotomic_field(k3.gen())
##             zeta15^5
##             sage: k15._coerce_from_other_cyclotomic_field(k3.gen()^2 + 17/3)
##             -zeta15^5 + 14/3
##             sage: k3._coerce_from_other_cyclotomic_field(k15.gen()^5)
##             zeta3
##             sage: k3._coerce_from_other_cyclotomic_field(-2/3 * k15.gen()^5 + 2/3)
##             -2/3*zeta3 + 2/3
##         """

##         K = x.parent()

##         if K is self:
##             return x
##         elif K == self:
##             return self._element_class(self, x.polynomial())
##         n = K.zeta_order()
##         m = self.zeta_order()
##         print n, m, x


##         self_gen = self.gen()

##         if m % n == 0:   # easy case
##             # pass this off to a method in the element class
##             # it can be done very quickly and easily by the cython<->NTL
##             # interface there
##             return x._lift_cyclotomic_element(self)

##         # Whatever happens below, it has to be consistent with
##         #  zeta_r |---> (zeta_s)^m

##         if m % 2 and not n%2:
##             m *= 2
##             self_gen = -self_gen

##         if only_canonical and m % n:
##             raise TypeError, "no canonical coercion"

##         if not is_CyclotomicField(K):
##             raise TypeError, "x must be in a cyclotomic field"

##         v = x.list()

##         # Find the smallest power r >= 1 of the generator g of K that is in self,
##         # i.e., find the smallest r such that g^r has order dividing m.

##         d = sage.rings.arith.gcd(m,n)
##         r = n // d

##         # Since we use the power basis for cyclomotic fields, if every
##         # v[i] with i not divisible by r is 0, then we're good.

##         # If h generates self and has order m, then the element g^r
##         # maps to the power of self of order gcd(m,n)., i.e., h^(m/gcd(m,n))
##         #
##         z = self_gen**(m // d)
##         w = self(1)

##         a = self(0)
##         for i in range(len(v)):
##             if i%r:
##                 if v[i]:
##                     raise TypeError, "element does not belong to cyclotomic field"
##             else:
##                 a += w*v[i]
##                 w *= z
##         return a

    def _coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
        """
        Coerce an element x of a cyclotomic field into self, if at all
        possible.

        INPUT:


        -  ``x`` - number field element

        -  ``only_canonical`` - bool (default: False); Attempt
           to work, even in some cases when x is not in a subfield of the
           cyclotomics (as long as x is a root of unity).


        EXAMPLES::

            sage: K = CyclotomicField(24) ; L = CyclotomicField(48)
            sage: L._coerce_from_other_cyclotomic_field(K.0+1)
            zeta48^2 + 1
            sage: K(L.0**2)
            zeta24
        """
        K = x.parent()
        if K is self:
            return x
        elif K == self:
            return self._element_class(self, x.polynomial())
        n = K._n()
        m = self._n()
        if m % n == 0:   # easy case
            # pass this off to a method in the element class
            # it can be done very quickly and easily by the
            # Cython<->NTL interface there
            return x._lift_cyclotomic_element(self)
        else:
            if only_canonical:
                raise TypeError
            n = x.multiplicative_order()
            m = self.zeta_order()
            if m % n == 0:
                # Harder case.  E.g., x = (zeta_42)^7 and
                # self.__zeta = zeta_6, so it is possible to
                # coerce x in, but not zeta_42 in.
                # Algorithm:
                #    1. Compute self.__zeta as an element
                #       of K = parent of x.  Call this y.
                #    2. Write x as a power r of y.
                #       TODO: we do step two STUPIDLY.
                #    3. Return self.__zeta to the power r.
                y = K(self.zeta(m))
                z = y
                for r in xrange(y.multiplicative_order()):
                    if z == x:
                        return self.zeta(m)**(r+1)
                    z *= y
            raise TypeError, "Cannot coerce %s into %s"%(x,self)
        return self._element_class(self, g)


    def _coerce_from_gap(self, x):
        """
        Attempt to coerce a GAP number field element into this cyclotomic
        field.

        EXAMPLES::

            sage: k5.<z> = CyclotomicField(5)
            sage: gap('E(5)^7 + 3')
            -3*E(5)-2*E(5)^2-3*E(5)^3-3*E(5)^4
            sage: w = gap('E(5)^7 + 3')
            sage: z^7 + 3
            z^2 + 3
            sage: k5(w) # indirect doctest
            z^2 + 3

        It may be that GAP uses a name for the generator of the cyclotomic field.
        We can deal with this case, if this name coincides with the name in Sage::

            sage: F=CyclotomicField(8)
            sage: z=F.gen()
            sage: a=gap(z+1/z); a
            E(8)-E(8)^3
            sage: F(a)
            -zeta8^3 + zeta8

        Matrices over cyclotomic fields are correctly dealt with it as well::

            sage: b=gap(Matrix(F,[[z^2,1],[0,a+1]])); b
            [ [ E(4), 1 ], [ 0, 1+E(8)-E(8)^3 ] ]
            sage: b[1,2]
            1
            sage: F(b[1,2])
            1
            sage: Matrix(b,F)
            [             zeta8^2                    1]
            [                   0 -zeta8^3 + zeta8 + 1]
        """
        s = str(x)
        i = s.find('E(')
        if i == -1:
            try:
                # it may be that a number field element's string representation
                # in GAP has an exclamation mark in it.
                return self(rational.Rational(s.replace('!','')))
            except TypeError:
                # There is no 'E(...)' in the string representation. But it may
                # be that 'E(...)' was overwritten in GAP. We can only hope that
                # by coincidence the name in GAP is the same as the name in self
                return self._coerce_from_str(s.replace('!',''))
        j = i + s[i:].find(')')
        n = int(s[i+2:j])
        if n == self.zeta_order():
            K = self
        else:
            K = CyclotomicField(n)
        zeta = K.gen()
        zeta_name = K.variable_name()
        while zeta_name in s: # could be that gap uses the generator name for a different purpose
            zeta_name = zeta_name+'_'
        s = s.replace('E(%s)'%n,zeta_name)
        s = sage.misc.all.sage_eval(s, locals={zeta_name:zeta})
        if K is self:
            return s
        else:
            return self(s)

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from the cyclotomic field self to
        the number field codomain.

        The cat option is currently ignored.

        EXAMPLES:

        This function is implicitly called by the Hom method or
        function.

        ::

            sage: K.<a> = NumberField(x^2 + 3); K
            Number Field in a with defining polynomial x^2 + 3
            sage: CyclotomicField(3).Hom(K) # indirect doctest
            Set of field embeddings from Cyclotomic Field of order 3 and degree 2 to Number Field in a with defining polynomial x^2 + 3
            sage: End(CyclotomicField(21))
            Automorphism group of Cyclotomic Field of order 21 and degree 12
        """
        if is_NumberFieldHomsetCodomain(codomain):
            import morphism
            return morphism.CyclotomicFieldHomset(self, codomain)
        else:
            raise TypeError

    def is_galois(self):
        """
        Return True since all cyclotomic fields are automatically Galois.

        EXAMPLES::

            sage: CyclotomicField(29).is_galois()
            True
        """
        return True

    def is_isomorphic(self, other):
        """
        Return True if the cyclotomic field self is isomorphic as a number
        field to other.

        EXAMPLES::

            sage: CyclotomicField(11).is_isomorphic(CyclotomicField(22))
            True
            sage: CyclotomicField(11).is_isomorphic(CyclotomicField(23))
            False
            sage: CyclotomicField(3).is_isomorphic(NumberField(x^2 + x +1, 'a'))
            True
            sage: CyclotomicField(18).is_isomorphic(CyclotomicField(9))
            True
            sage: CyclotomicField(10).is_isomorphic(NumberField(x^4 - x^3 + x^2 - x + 1, 'b'))
            True

        Check :trac:`14300`::

            sage: K = CyclotomicField(4)
            sage: N = K.extension(x^2-5, 'z')
            sage: K.is_isomorphic(N)
            False
            sage: K.is_isomorphic(CyclotomicField(8))
            False
        """
        if is_CyclotomicField(other):
            return self.zeta_order() == other.zeta_order()
        return NumberField_generic.is_isomorphic(self, other)

    def complex_embedding(self, prec=53):
        r"""
        Return the embedding of this cyclotomic field into the approximate
        complex field with precision prec obtained by sending the generator
        `\zeta` of self to exp(2\*pi\*i/n), where `n` is
        the multiplicative order of `\zeta`.

        EXAMPLES::

            sage: C = CyclotomicField(4)
            sage: C.complex_embedding()
            Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Field with 53 bits of precision
              Defn: zeta4 |--> 6.12323399573677e-17 + 1.00000000000000*I

        Note in the example above that the way zeta is computed (using sin
        and cosine in MPFR) means that only the prec bits of the number
        after the decimal point are valid.

        ::

            sage: K = CyclotomicField(3)
            sage: phi = K.complex_embedding(10)
            sage: phi(K.0)
            -0.50 + 0.87*I
            sage: phi(K.0^3)
            1.0
            sage: phi(K.0^3 - 1)
            0.00
            sage: phi(K.0^3 + 7)
            8.0
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        return self.hom([CC.zeta(self._n())], check=False)

    def complex_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the approximate
        complex field with precision prec.

        If you want 53-bit double precision, which is faster but less
        reliable, then do ``self.embeddings(CDF)``.

        EXAMPLES::

            sage: CyclotomicField(5).complex_embeddings()
            [
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Field with 53 bits of precision
              Defn: zeta5 |--> 0.309016994374947 + 0.951056516295154*I,
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Field with 53 bits of precision
              Defn: zeta5 |--> -0.809016994374947 + 0.587785252292473*I,
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Field with 53 bits of precision
              Defn: zeta5 |--> -0.809016994374947 - 0.587785252292473*I,
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Field with 53 bits of precision
              Defn: zeta5 |--> 0.309016994374947 - 0.951056516295154*I
            ]
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        try:
            return self.__embeddings[CC]
        except AttributeError:
            self.__embeddings = {}
        except KeyError:
            pass
        n = self._n()
        z = CC.zeta(n)
        X = [m for m in range(n) if sage.rings.arith.gcd(m,n) == 1]
        v = [self.hom([z**n], check=False) for n in X]
        self.__embeddings[CC] = Sequence(v, cr=True, immutable=True,
                                         check=False, universe=self.Hom(CC))
        return self.__embeddings[CC]

    def real_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the approximate
        real field with precision prec.

        Mostly, of course, there are no such embeddings.

        EXAMPLES::

            sage: CyclotomicField(4).real_embeddings()
            []
            sage: CyclotomicField(2).real_embeddings()
            [
            Ring morphism:
              From: Cyclotomic Field of order 2 and degree 1
              To:   Real Field with 53 bits of precision
              Defn: -1 |--> -1.00000000000000
            ]
        """
        K = sage.rings.real_mpfr.RealField(prec)
        n = self._n()
        if n > 2:
            return Sequence([], cr=False, immutable=True,
                                         check=False, universe=self.Hom(K))
        else:
            return self.embeddings(K)

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this cyclotomic field,
        respectively.

        Trivial since, apart from QQ, cyclotomic fields are totally
        complex.

        EXAMPLES::

            sage: CyclotomicField(5).signature()
            (0, 2)
            sage: CyclotomicField(2).signature()
            (1, 0)
        """
        m = ZZ(self.degree())
        if m == 1:
            return (ZZ(1), ZZ(0))
        else:
            return (ZZ(0), ZZ(m/2))

    def different(self):
        """
        Returns the different ideal of the cyclotomic field self.

        EXAMPLES::

            sage: C20 = CyclotomicField(20)
            sage: C20.different()
            Fractional ideal (10, 2*zeta20^6 - 4*zeta20^4 - 4*zeta20^2 + 2)
            sage: C18 = CyclotomicField(18)
            sage: D = C18.different().norm()
            sage: D == C18.discriminant().abs()
            True
        """
        try:
            return self.__different

        except AttributeError:

            z = self.gen()
            n = self._n()
            D = self.ideal(1)
            factors = n.factor()
            for f in factors:
                p = f[0]
                r = f[1]
                e = (r*p - r - 1)*p**(r-1)
                D *= self.ideal(z**(n/p**r) - 1)**e
            self.__different = D
            return self.__different

    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the cyclotomic
        field self, or if v is specified, the determinant of the trace
        pairing on the elements of the list v.

        Uses the formula for the discriminant of a prime power cyclotomic
        field and Hilbert Theorem 88 on the discriminant of composita.

        INPUT:


        -  ``v (optional)`` - list of element of this number
           field


        OUTPUT: Integer if v is omitted, and Rational otherwise.

        EXAMPLES::

            sage: CyclotomicField(20).discriminant()
            4000000
            sage: CyclotomicField(18).discriminant()
            -19683
        """
        if v == None:
            try:
                return self.__disc
            except AttributeError:
                n = self._n()
                deg = self.degree()
                d = ZZ(1) # so that CyclotomicField(1).disc() has the right type
                factors = n.factor()
                for f in factors:
                    p = f[0]
                    r = f[1]
                    e = (r*p - r - 1)*deg/(p-1)
                    d *= p**e
                sign = 1
                if len(factors) == 1 and (n == 4 or factors[0][0].mod(4) == 3):
                    sign = -1
                elif len(factors) == 2 and factors[0] == (2, 1) and factors[1][0].mod(4) == 3:
                    sign = -1
                self.__disc = sign*d
                return self.__disc
        else:
            return NumberField_generic.discriminant(self, v)


    def next_split_prime(self, p=2):
        """
        Return the next prime integer `p` that splits completely in
        this cyclotomic field (and does not ramify).

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: K.next_split_prime(7)
            13
        """
        n = self._n()
        while True:
            p = sage.rings.arith.next_prime(p)
            if p % n == 1:
                return p

    def _pari_integral_basis(self, v=None, important=True):
        """
        Internal function returning an integral basis of this number field in
        PARI format.

        This field is cyclomotic, so this is a trivial computation, since
        the power basis on the generator is an integral basis. Thus the ``v``
        and ``important`` parameters are ignored.

        EXAMPLES::

            sage: CyclotomicField(5)._pari_integral_basis()
            [1, y, y^2, y^3]
            sage: len(CyclotomicField(137)._pari_integral_basis())
            136
        """
        try:
            return self._integral_basis_dict[tuple()]
        except KeyError:
            z = pari(self.gen())
            a = pari(1)
            B = []
            for n in xrange(self.degree()):
                B.append(a.lift())
                a *= z
            self._integral_basis_dict[tuple()] = pari(B)
            return B


    def zeta_order(self):
        """
        Return the order of the maximal root of unity contained in this
        cyclotomic field.

        EXAMPLES::

            sage: CyclotomicField(1).zeta_order()
            2
            sage: CyclotomicField(4).zeta_order()
            4
            sage: CyclotomicField(5).zeta_order()
            10
            sage: CyclotomicField(5)._n()
            5
            sage: CyclotomicField(389).zeta_order()
            778
        """
        return self.__zeta_order

    def _multiplicative_order_table(self):
        """
        Return a dictionary that maps powers of zeta to their order. This
        makes computing the orders of the elements of finite order in this
        field faster.

        EXAMPLES::

            sage: v = CyclotomicField(6)._multiplicative_order_table()
            sage: w = v.items(); w.sort(); w
            [(-1, 2), (1, 1), (-x, 3), (-x + 1, 6), (x - 1, 3), (x, 6)]
        """
        try:
            return self.__multiplicative_order_table
        except AttributeError:
            t = {}
            x = self(1)
            n = self.zeta_order()
            m = 0
            zeta = self.zeta(n)
            # todo: this desperately needs to be optimized!!!
            for i in range(n):
                t[x.polynomial()] = n//arith.GCD(m,n)   # multiplicative_order of (zeta_n)**m
                x *= zeta
                m += 1
            self.__multiplicative_order_table = t
            return t

    def zeta(self, n=None, all=False):
        """
        Returns an element of multiplicative order `n` in this this
        cyclotomic field, if there is one. Raises a ValueError if there is
        not.

        INPUT:


        -  ``n`` - integer (default: None, returns element of
           maximal order)

        -  ``all`` - bool (default: False) - whether to return
           a list of all n-th roots.


        OUTPUT: root of unity or list

        EXAMPLES::

            sage: k = CyclotomicField(7)
            sage: k.zeta()
            zeta7
            sage: k.zeta().multiplicative_order()
            7
            sage: k = CyclotomicField(49)
            sage: k.zeta().multiplicative_order()
            49
            sage: k.zeta(7).multiplicative_order()
            7
            sage: k.zeta()
            zeta49
            sage: k.zeta(7)
            zeta49^7

        ::

            sage: K.<a> = CyclotomicField(7)
            sage: K.zeta(14, all=True)
            [-a^4, -a^5, a^5 + a^4 + a^3 + a^2 + a + 1, -a, -a^2, -a^3]
            sage: K.<a> = CyclotomicField(10)
            sage: K.zeta(20, all=True)
            Traceback (most recent call last):
            ...
            ValueError: n (=20) does not divide order of generator

        ::

            sage: K.<a> = CyclotomicField(5)
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: n (=4) does not divide order of generator
            sage: v = K.zeta(5, all=True); v
            [a, a^2, a^3, -a^3 - a^2 - a - 1]
            sage: [b^5 for b in v]
            [1, 1, 1, 1]
        """
        if n is None:
            return self.gen()
        else:
            n = integer.Integer(n)
            z = self.gen()
            m = z.multiplicative_order()
            if n % 2 == 0 and m % 2 == 1:
                # In the n-th cyclotomic field, n odd, there are
                # actually 2*n-th roots of unity, so we include them.
                z = -z**((m+1)//2) # -z
                m = 2*m
            if m % n != 0:
                raise ValueError, "n (=%s) does not divide order of generator"%n
                # use generic method (factor cyclotomic polynomial)
                #  -- this is potentially really slow, so don't do it.
                #return field.Field.zeta(self, n, all=all)
            a = z**(m//n)
            if all:
                v = [a]
                b = a*a
                for i in range(2,n):
                    if sage.rings.arith.gcd(i, n) == 1:
                        v.append(b)
                    b = b * a
                return v
            else:
                return a

    def number_of_roots_of_unity(self):
        """
        Return number of roots of unity in this cyclotomic field.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(21)
            sage: K.number_of_roots_of_unity()
            42
        """
        n = self._n()
        if n%2:
            n *= 2
        return n

    def roots_of_unity(self):
        """
        Return all the roots of unity in this cyclotomic field, primitive
        or not.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(3)
            sage: zs = K.roots_of_unity(); zs
            [1, a, -a - 1, -1, -a, a + 1]
            sage: [ z**K.number_of_roots_of_unity() for z in zs ]
            [1, 1, 1, 1, 1, 1]
        """
        z = self.gen()
        n = self._n()
        v = [z**k for k in range(n)]
        if n%2:
            v += [-x for x in v]
        return v


class NumberField_quadratic(NumberField_absolute):
    r"""
    Create a quadratic extension of the rational field.

    The command ``QuadraticField(a)`` creates the field `\QQ(\sqrt{a})`.

    EXAMPLES::

        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3
        sage: QuadraticField(-4, 'b')
        Number Field in b with defining polynomial x^2 + 4
    """
    def __init__(self, polynomial, name=None, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None):
        """
        Create a quadratic number field.

        EXAMPLES::

            sage: k.<a> = QuadraticField(5, check=False); k
            Number Field in a with defining polynomial x^2 - 5

        Don't do this::

            sage: k.<a> = QuadraticField(4, check=False); k
            Number Field in a with defining polynomial x^2 - 4

        TESTS::

            sage: k.<a> = QuadraticField(7)
            sage: type(k.zero_element())
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: type(k.one_element())
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>

            sage: TestSuite(k).run()
        """
        NumberField_absolute.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
        self._standard_embedding = True
        self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic
        c, b, a = [rational.Rational(t) for t in self.defining_polynomial().list()]
        # set the generator
        Dpoly = b*b - 4*a*c
        D = (Dpoly.numer() * Dpoly.denom()).squarefree_part(bound=10000)
        self._D = D
        parts = -b/(2*a), (Dpoly/D).sqrt()/(2*a)
        self._NumberField_generic__gen = self._element_class(self, parts)

        # we must set the flag _standard_embedding *before* any element creation
        # Note that in the following code, no element is built.
        emb = self.coerce_embedding()
        if emb is not None:
            rootD = number_field_element_quadratic.NumberFieldElement_quadratic(self, (QQ(0),QQ(1)))
            if D > 0:
                from sage.rings.real_double import RDF
                self._standard_embedding = RDF.has_coerce_map_from(self) and RDF(rootD) > 0
            else:
                from sage.rings.complex_double import CDF
                self._standard_embedding = CDF.has_coerce_map_from(self) and CDF(rootD).imag() > 0

        # we reset _NumberField_generic__gen has the flag standard_embedding
        # might be modified
        self._NumberField_generic__gen = self._element_class(self, parts)


        # NumberField_absolute.__init__(...) set _zero_element and
        # _one_element to NumberFieldElement_absolute values, which is
        # wrong (and dangerous; such elements can actually be used to
        # crash Sage: see #5316).  Overwrite them with correct values.
        self._zero_element = self._element_class(self, (QQ(0), QQ(0)))
        self._one_element  = self._element_class(self, (QQ(1), QQ(0)))

    def _coerce_map_from_(self, K):
        """
        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: f = K.coerce_map_from(QQ); f # indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3
            sage: f(3/5)
            3/5
            sage: parent(f(3/5)) is K
            True

            sage: g = K.coerce_map_from(ZZ); g # indirect doctest
            Natural morphism:
              From: Integer Ring
              To:   Number Field in a with defining polynomial x^2 + 3
            sage: g(1)
            1
            sage: parent(g(1)) is K
            True
        """
        if K is ZZ:
            return number_field_element_quadratic.Z_to_quadratic_field_element(self)
        if K is QQ:
            return number_field_element_quadratic.Q_to_quadratic_field_element(self)
        return NumberField_absolute._coerce_map_from_(self, K)

    def _latex_(self):
        """
        Return the latex representation of this quadratic field.

        EXAMPLES::

            sage: Z = QuadraticField(7)
            sage: latex(Z) # indirect doctest
            \Bold{Q}(\sqrt{7})

            sage: Z = QuadraticField(7, latex_name='x')
            sage: latex(Z) # indirect doctest
            \Bold{Q}[x]/(x^{2} - 7)
        """
        v = self.latex_variable_name()
        if v.startswith('\\sqrt'):
            return "%s(%s)"%(latex(QQ), v)
        else:
            return NumberField_generic._latex_(self)

    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the number
        field, or if v is specified, the determinant of the trace pairing
        on the elements of the list v.

        INPUT:


        -  ``v (optional)`` - list of element of this number
           field


        OUTPUT: Integer if v is omitted, and Rational otherwise.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: K.discriminant()
            -4
            sage: K.<a> = NumberField(x^2+5)
            sage: K.discriminant()
            -20
            sage: K.<a> = NumberField(x^2-5)
            sage: K.discriminant()
            5
        """
        if v is None:
            try:
                return self.__disc
            except AttributeError:
                d = self._D.squarefree_part()
                if d % 4 != 1:
                    d *= 4
                self.__disc = d
                return self.__disc
        else:
            return NumberField_generic.discriminant(self, v)

    def __reduce__(self):
        """
        This is used in pickling quadratic number fields.

        TESTS::

            sage: K.<z7> = QuadraticField(7)
            sage: L = loads(dumps(K))
            sage: L
            Number Field in z7 with defining polynomial x^2 - 7
            sage: L == K
            True
        """
        return NumberField_quadratic_v1, (self.polynomial(), self.variable_name(), self.gen_embedding())


    def is_galois(self):
        """
        Return True since all quadratic fields are automatically Galois.

        EXAMPLES::

            sage: QuadraticField(1234,'d').is_galois()
            True
        """
        return True

    def class_number(self, proof=None):
        r"""
        Return the size of the class group of self.

        If proof = False (*not* the default!) and the discriminant of the
        field is negative, then the following warning from the PARI manual
        applies:

        .. warning::

            For `D<0`, this function may give incorrect results when
            the class group has a low exponent (has many cyclic
            factors), because implementing Shank's method in full
            generality slows it down immensely.

        EXAMPLES::

            sage: QuadraticField(-23,'a').class_number()
            3

        These are all the primes so that the class number of
        `\QQ(\sqrt{-p})` is `1`::

            sage: [d for d in prime_range(2,300) if not is_square(d) and QuadraticField(-d,'a').class_number() == 1]
            [2, 3, 7, 11, 19, 43, 67, 163]

        It is an open problem to *prove* that there are infinity many
        positive square-free `d` such that
        `\QQ(\sqrt{d})` has class number `1`:

        ::

            sage: len([d for d in range(2,200) if not is_square(d) and QuadraticField(d,'a').class_number() == 1])
            121

        TESTS::

            sage: type(QuadraticField(-23,'a').class_number())
            <type 'sage.rings.integer.Integer'>
            sage: type(NumberField(x^3 + 23, 'a').class_number())
            <type 'sage.rings.integer.Integer'>
            sage: type(NumberField(x^3 + 23, 'a').extension(x^2 + 5, 'b').class_number())
            <type 'sage.rings.integer.Integer'>
            sage: type(CyclotomicField(10).class_number())
            <type 'sage.rings.integer.Integer'>
        """
        proof = proof_flag(proof)
        try:
            return self.__class_number
        except AttributeError:
            D = self.discriminant()
            if D < 0 and proof:
                self.__class_number = ZZ(pari("qfbclassno(%s,1)"%D))
            else:
                self.__class_number = ZZ(pari("qfbclassno(%s)"%D))
            return self.__class_number

    def hilbert_class_field_defining_polynomial(self, name='x'):
        r"""
        Returns a polynomial over `\QQ` whose roots generate the
        Hilbert class field of this quadratic field as an extension of
        this quadratic field.

        .. note::

           Computed using PARI via Schertz's method. This
           implementation is quite fast.

        EXAMPLES::

            sage: K.<b> = QuadraticField(-23)
            sage: K.hilbert_class_field_defining_polynomial()
            x^3 - x^2 + 1

        Note that this polynomial is not the actual Hilbert class
        polynomial: see ``hilbert_class_polynomial``::

            sage: K.hilbert_class_polynomial()
            x^3 + 3491750*x^2 - 5151296875*x + 12771880859375

        ::

            sage: K.<a> = QuadraticField(-431)
            sage: K.class_number()
            21
            sage: K.hilbert_class_field_defining_polynomial(name='z')
            z^21 + 6*z^20 + 9*z^19 - 4*z^18 + 33*z^17 + 140*z^16 + 220*z^15 + 243*z^14 + 297*z^13 + 461*z^12 + 658*z^11 + 743*z^10 + 722*z^9 + 681*z^8 + 619*z^7 + 522*z^6 + 405*z^5 + 261*z^4 + 119*z^3 + 35*z^2 + 7*z + 1
        """
        f = pari(self.discriminant()).quadhilbert()
        return QQ[name](f)

    def hilbert_class_field(self, names):
        r"""
        Returns the Hilbert class field of this quadratic field as a
        relative extension of this field.

        .. note::

           For the polynomial that defines this field as a relative
           extension, see the ``hilbert_class_field_defining_polynomial``
           command, which is vastly faster than this command, since it doesn't
           construct a relative extension.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: L = K.hilbert_class_field('b'); L
            Number Field in b with defining polynomial x^3 - x^2 + 1 over its base field
            sage: L.absolute_field('c')
            Number Field in c with defining polynomial x^6 - 2*x^5 + 70*x^4 - 90*x^3 + 1631*x^2 - 1196*x + 12743
            sage: K.hilbert_class_field_defining_polynomial()
            x^3 - x^2 + 1
        """
        f = self.hilbert_class_field_defining_polynomial()
        return self.extension(f, names)

    def hilbert_class_polynomial(self, name='x'):
        r"""
        Compute the Hilbert class polynomial of this quadratic field.

        Right now, this is only implemented for imaginary quadratic
        fields.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: K.hilbert_class_polynomial()
            x

            sage: K.<a> = QuadraticField(-31)
            sage: K.hilbert_class_polynomial(name='z')
            z^3 + 39491307*z^2 - 58682638134*z + 1566028350940383
        """
        D = self.discriminant()

        if D > 0:
            raise NotImplementedError, "Hilbert class polynomial is not implemented for real quadratic fields."

        from sage.schemes.elliptic_curves.all import hilbert_class_polynomial as HCP
        return QQ[name](HCP(D))

def is_fundamental_discriminant(D):
    r"""
    Return True if the integer `D` is a fundamental
    discriminant, i.e., if `D \cong 0,1\pmod{4}`, and
    `D\neq 0, 1` and either (1) `D` is square free or
    (2) we have `D\cong 0\pmod{4}` with
    `D/4 \cong 2,3\pmod{4}` and `D/4` square free. These
    are exactly the discriminants of quadratic fields.

    EXAMPLES::

        sage: [D for D in range(-15,15) if is_fundamental_discriminant(D)]
        [-15, -11, -8, -7, -4, -3, 5, 8, 12, 13]
        sage: [D for D in range(-15,15) if not is_square(D) and QuadraticField(D,'a').disc() == D]
        [-15, -11, -8, -7, -4, -3, 5, 8, 12, 13]
    """
    d = D % 4
    if not (d in [0,1]):
        return False
    return D != 1 and  D != 0 and \
           (arith.is_squarefree(D) or \
            (d == 0 and (D//4)%4 in [2,3] and arith.is_squarefree(D//4)))


###################
# For pickling
###################


def NumberField_absolute_v1(poly, name, latex_name, canonical_embedding=None):
    """
    This is used in pickling generic number fields.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import NumberField_absolute_v1
        sage: R.<x> = QQ[]
        sage: NumberField_absolute_v1(x^2 + 1, 'i', 'i')
        Number Field in i with defining polynomial x^2 + 1
    """
    return NumberField_absolute(poly, name, latex_name, check=False, embedding=canonical_embedding)

NumberField_generic_v1 = NumberField_absolute_v1  # for historical reasons only (so old objects unpickle)

def NumberField_cyclotomic_v1(zeta_order, name, canonical_embedding=None):
    """
    This is used in pickling cyclotomic fields.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import NumberField_cyclotomic_v1
        sage: NumberField_cyclotomic_v1(5,'a')
        Cyclotomic Field of order 5 and degree 4
        sage: NumberField_cyclotomic_v1(5,'a').variable_name()
        'a'
    """
    return NumberField_cyclotomic(zeta_order, name, embedding=canonical_embedding)

def NumberField_quadratic_v1(poly, name, canonical_embedding=None):
    """
    This is used in pickling quadratic fields.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import NumberField_quadratic_v1
        sage: R.<x> = QQ[]
        sage: NumberField_quadratic_v1(x^2 - 2, 'd')
        Number Field in d with defining polynomial x^2 - 2
    """
    return NumberField_quadratic(poly, name, check=False, embedding=canonical_embedding)


def put_natural_embedding_first(v):
    """
    Helper function for embeddings() functions for number fields.

    INPUT: a list of embeddings of a number field

    OUTPUT: None. The
    list is altered in-place, so that, if possible, the first embedding
    has been switched with one of the others, so that if there is an
    embedding which preserves the generator names then it appears
    first.

    EXAMPLES::

        sage: K.<a> = CyclotomicField(7)
        sage: embs = K.embeddings(K)
        sage: [e(a) for e in embs] # random - there is no natural sort order
        [a, a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1]
        sage: id = [ e for e in embs if e(a) == a ][0]; id
        Ring endomorphism of Cyclotomic Field of order 7 and degree 6
          Defn: a |--> a
        sage: permuted_embs = list(embs); permuted_embs.remove(id); permuted_embs.append(id)
        sage: [e(a) for e in permuted_embs] # random - but natural map is not first
        [a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1, a]
        sage: permuted_embs[0] != a
        True
        sage: from sage.rings.number_field.number_field import put_natural_embedding_first
        sage: put_natural_embedding_first(permuted_embs)
        sage: [e(a) for e in permuted_embs] # random - but natural map is first
        [a, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1, a^2]
        sage: permuted_embs[0] == id
        True
    """
    for i in range(len(v)):
        phi = v[i]
        a = str(list(phi.domain().gens()))
        b = str(list(phi.im_gens()))
        if a == b:
            v[i] = v[0]
            v[0] = phi
            return



def refine_embedding(e, prec=None):
    r"""
    Given an embedding from a number field to either `\RR` or
    `\CC`, returns an equivalent embedding with higher precision.

    INPUT:


    -  ``e`` - an embedding of a number field into either
       RR or CC (with some precision)

    - ``prec`` - (default None) the desired precision; if None,
       current precision is doubled; if Infinity, the equivalent
       embedding into either ``QQbar`` or ``AA`` is returned.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import refine_embedding
        sage: K = CyclotomicField(3)
        sage: e10 = K.complex_embedding(10)
        sage: e10.codomain().precision()
        10
        sage: e25 = refine_embedding(e10, prec=25)
        sage: e25.codomain().precision()
        25

    An example where we extend a real embedding into ``AA``::

        sage: K.<a> = NumberField(x^3-2)
        sage: K.signature()
        (1, 1)
        sage: e = K.embeddings(RR)[0]; e
        Ring morphism:
        From: Number Field in a with defining polynomial x^3 - 2
        To:   Real Field with 53 bits of precision
        Defn: a |--> 1.25992104989487
        sage: e = refine_embedding(e,Infinity); e
        Ring morphism:
        From: Number Field in a with defining polynomial x^3 - 2
        To:   Algebraic Real Field
        Defn: a |--> 1.259921049894873?

    Now we can obtain arbitrary precision values with no trouble::

        sage: RealField(150)(e(a))
        1.2599210498948731647672106072782283505702515
        sage: _^3
        2.0000000000000000000000000000000000000000000
        sage: RealField(200)(e(a^2-3*a+7))
        4.8076379022835799804500738174376232086807389337953290695624

    Complex embeddings can be extended into ``QQbar``::

        sage: e = K.embeddings(CC)[0]; e
        Ring morphism:
        From: Number Field in a with defining polynomial x^3 - 2
        To:   Complex Field with 53 bits of precision
        Defn: a |--> -0.62996052494743... - 1.09112363597172*I
        sage: e = refine_embedding(e,Infinity); e
        Ring morphism:
        From: Number Field in a with defining polynomial x^3 - 2
        To:   Algebraic Field
        Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I
        sage: ComplexField(200)(e(a))
        -0.62996052494743658238360530363911417528512573235075399004099 - 1.0911236359717214035600726141898088813258733387403009407036*I
        sage: e(a)^3
        2

    Embeddings into lazy fields work::

        sage: L = CyclotomicField(7)
        sage: x = L.specified_complex_embedding(); x
        Generic morphism:
          From: Cyclotomic Field of order 7 and degree 6
          To:   Complex Lazy Field
          Defn: zeta7 -> 0.623489801858734? + 0.781831482468030?*I
        sage: refine_embedding(x, 300)
        Ring morphism:
          From: Cyclotomic Field of order 7 and degree 6
          To:   Complex Field with 300 bits of precision
          Defn: zeta7 |--> 0.623489801858733530525004884004239810632274730896402105365549439096853652456487284575942507 + 0.781831482468029808708444526674057750232334518708687528980634958045091731633936441700868007*I
        sage: refine_embedding(x, infinity)
        Ring morphism:
          From: Cyclotomic Field of order 7 and degree 6
          To:   Algebraic Field
          Defn: zeta7 |--> 0.6234898018587335? + 0.7818314824680299?*I
    """
    K = e.domain()
    RC = e.codomain()
    if RC in (sage.rings.qqbar.AA, sage.rings.qqbar.QQbar):
        return e
    if RC in (RLF, CLF):
        prec_old = e.gen_image().approx().prec()
        old_root = e(K.gen()).approx()
    else:
        prec_old = RC.precision()
        old_root = e(K.gen())

    if prec is None:
        prec = 2*prec_old
    elif prec_old >= prec:
        return e

    # We first compute all the embeddings at the new precision:
    if sage.rings.real_mpfr.is_RealField(RC) or RC is RDF:
        if prec==sage.rings.infinity.Infinity:
            elist = K.embeddings(sage.rings.qqbar.AA)
        else:
            elist = K.real_embeddings(prec)
    else:
        if prec==sage.rings.infinity.Infinity:
            elist = K.embeddings(sage.rings.qqbar.QQbar)
        else:
            elist = K.complex_embeddings(prec)

    # Now we determine which is an extension of the old one; this
    # relies on the fact that coercing a high-precision root into a
    # field with lower precision will equal the lower-precision root!
    diffs = [(RC(ee(K.gen()))-old_root).abs() for ee in elist]
    return elist[min(izip(diffs,count()))[1]]
