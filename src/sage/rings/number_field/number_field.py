# sage.doctest: needs sage.libs.linbox
r"""
Number fields

We define a quartic number field and its quadratic extension::

    sage: x = polygen(ZZ, 'x')
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

.. WARNING::

   Doing arithmetic in towers of relative fields that depends on
   canonical coercions is currently VERY SLOW. It is much better to
   explicitly coerce all elements into a common field, then do
   arithmetic with them there (which is quite fast).

AUTHORS:

- William Stein (2004, 2005): initial version
- Steven Sivek (2006-05-12): added support for relative extensions
- William Stein (2007-09-04): major rewrite and documentation
- Robert Bradshaw (2008-10): specified embeddings into ambient fields
- Simon King (2010-05): improved coercion from GAP
- Jeroen Demeyer (2010-07, 2011-04): upgraded PARI (:issue:`9343`, :issue:`10430`, :issue:`11130`)
- Robert Harron (2012-08): added is_CM(), complex_conjugation(), and
  maximal_totally_real_subfield()
- Christian Stump (2012-11): added conversion to universal cyclotomic field
- Julian Rueth (2014-04-03): absolute number fields are unique parents
- Vincent Delecroix (2015-02): comparisons/floor/ceil using embeddings
- Kiran Kedlaya (2016-05): relative number fields hash based on relative polynomials
- Peter Bruin (2016-06): made number fields fully satisfy unique representation
- John Jones (2017-07): improved check for is_galois(), add is_abelian(), building on work in patch by Chris Wuthrich
- Anna Haensch (2018-03): added :meth:`quadratic_defect`
- Michael Daub, Chris Wuthrich (2020-09-01): added Dirichlet characters for abelian fields
"""
# ****************************************************************************
#       Copyright (C) 2004-2007 William Stein <wstein@gmail.com>
#                     2014-2022 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from itertools import count
from collections import Counter

from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecation

import sage.libs.ntl.all as ntl
import sage.rings.abc
import sage.rings.complex_mpfr
from sage.rings.polynomial.polynomial_element import Polynomial
import sage.rings.real_mpfr
import sage.rings.real_mpfi
import sage.rings.complex_double
import sage.rings.real_double
import sage.rings.real_lazy

from sage.arith.misc import euler_phi, factor, gcd, next_prime
from sage.misc.fast_methods import WithEqualityById
from sage.misc.functional import is_odd, lift
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.sage_eval import sage_eval
from sage.rings.infinity import Infinity
from sage.rings.finite_rings.integer_mod import mod
from sage.categories.number_fields import NumberFields

from sage.misc.latex import latex_variable_name

from .unit_group import UnitGroup
from .class_group import ClassGroup
from .class_group import SClassGroup

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.sequence import Sequence
from sage.structure.factorization import Factorization
from sage.structure.category_object import normalize_names
import sage.structure.parent_gens
import sage.structure.coerce_exceptions

from sage.structure.proof.proof import get_flag
from . import maps
from . import structure
from . import number_field_morphisms

from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.real_mpfr import RR

from sage.interfaces.abc import GapElement
from sage.rings.number_field.morphism import RelativeNumberFieldHomomorphism_from_abs

lazy_import('sage.libs.gap.element', 'GapElement', as_='LibGapElement')
lazy_import('sage.rings.universal_cyclotomic_field', 'UniversalCyclotomicFieldElement')


_NumberFields = NumberFields()


def is_NumberFieldHomsetCodomain(codomain):
    """
    Return whether ``codomain`` is a valid codomain for a number
    field homset.

    This is used by NumberField._Hom_ to determine
    whether the created homsets should be a
    :class:`sage.rings.number_field.homset.NumberFieldHomset`.

    EXAMPLES:

    This currently accepts any parent (CC, RR, ...) in :class:`Fields`::

        sage: from sage.rings.number_field.number_field import is_NumberFieldHomsetCodomain
        sage: is_NumberFieldHomsetCodomain(QQ)
        True
        sage: x = polygen(ZZ, 'x')
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

        sage: is_NumberFieldHomsetCodomain(gap.Rationals)                               # needs sage.libs.gap
        False
    """
    from sage.categories.fields import Fields
    return codomain in Fields()


def proof_flag(t):
    """
    Used for easily determining the correct proof flag to use.

    Return ``t`` if ``t`` is not ``None``, otherwise return the system-wide
    proof-flag for number fields (default: ``True``).

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


from sage.misc.latex import latex

import sage.rings.infinity as infinity
from sage.rings.rational import Rational
from sage.rings.integer import Integer
import sage.rings.polynomial.polynomial_element as polynomial_element
import sage.groups.abelian_gps.abelian_group
import sage.rings.complex_interval_field

from sage.structure.factory import UniqueFactory
from . import number_field_element
from . import number_field_element_quadratic
from .number_field_ideal import NumberFieldIdeal, NumberFieldFractionalIdeal
from sage.libs.pari.all import pari, pari_gen

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfi import RIF
from sage.rings.cif import CIF
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.real_lazy import RLF, CLF
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing


def NumberField(polynomial, name=None, check=True, names=None, embedding=None,
                latex_name=None, assume_disc_small=False, maximize_at_primes=None, structure=None,
                *, latex_names=None, **kwds):
    r"""
    Return *the* number field (or tower of number fields) defined by the
    irreducible ``polynomial``.

    INPUT:

    - ``polynomial`` -- a polynomial over `\QQ` or a number field, or a list
      of such polynomials
    - ``names`` (or ``name``) -- string or list of strings, the names of
      the generators
    - ``check`` -- boolean (default: ``True``); do type checking and
      irreducibility checking
    - ``embedding`` -- ``None``, an element, or a list of elements, the
      images of the generators in an ambient field (default: ``None``)
    - ``latex_names`` (or ``latex_name``) -- ``None``, a string, or a
      list of strings (default: ``None``); how the generators are printed
      for latex output
    - ``assume_disc_small`` -- boolean (default: ``False``); if ``True``,
      assume that no square of a prime greater than PARI's primelimit
      (which should be 500000). Only applies for absolute fields at
      present.
    - ``maximize_at_primes`` -- ``None`` or a list of primes (default:
      ``None``); if not ``None``, then the maximal order is computed by
      maximizing only at the primes in this list, which completely avoids
      having to factor the discriminant, but of course can lead to wrong
      results; only applies for absolute fields at present.
    - ``structure`` -- ``None``, a list or an instance of
      :class:`structure.NumberFieldStructure` (default: ``None``),
      internally used to pass in additional structural information, e.g.,
      about the field from which this field is created as a subfield.

    We accept ``implementation`` and ``prec`` attributes for compatibility
    with :class:`~sage.categories.pushout.AlgebraicExtensionFunctor`
    but we ignore them as they are not used.

    EXAMPLES::

        sage: z = QQ['z'].0
        sage: K = NumberField(z^2 - 2, 's'); K
        Number Field in s with defining polynomial z^2 - 2
        sage: s = K.0; s
        s
        sage: s*s
        2
        sage: s^2
        2

    Constructing a relative number field::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 2)
        sage: R.<t> = K[]
        sage: L.<b> = K.extension(t^3 + t + a); L
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

        sage: K.<a> = NumberField(x^3 - 5)
        sage: a^3
        5
        sage: L.<a> = NumberField(x^3 - 5)
        sage: K is L
        True

    Equality of number fields depends on the variable name of the
    defining polynomial::

        sage: x = polygen(QQ, 'x'); y = polygen(QQ, 'y')
        sage: k.<a> = NumberField(x^2 + 3)
        sage: m.<a> = NumberField(y^2 + 3)
        sage: k
        Number Field in a with defining polynomial x^2 + 3
        sage: m
        Number Field in a with defining polynomial y^2 + 3
        sage: k == m
        False

    In case of conflict of the generator name with the name given by the preparser, the name given by the preparser takes precedence::

        sage: K.<b> = NumberField(x^2 + 5, 'a'); K
        Number Field in b with defining polynomial x^2 + 5

    One can also define number fields with specified embeddings, may be used
    for arithmetic and deduce relations with other number fields which would
    not be valid for an abstract number field. ::

        sage: K.<a> = NumberField(x^3 - 2, embedding=1.2)
        sage: RR.coerce_map_from(K)
        Composite map:
          From: Number Field in a with defining polynomial x^3 - 2 with a = 1.259921049894873?
          To:   Real Field with 53 bits of precision
          Defn:   Generic morphism:
                  From: Number Field in a with defining polynomial x^3 - 2
                        with a = 1.259921049894873?
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
        sage: L.<b> = NumberField(x^6 - 2, embedding=1.1)
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

        sage: K.<a> = NumberField(x^3 - 2, embedding=CC.gen() - 0.6)
        sage: CC(a)
        -0.629960524947436 + 1.09112363597172*I

        sage: # needs sage.rings.padics
        sage: L = Qp(5)
        sage: f = polygen(L)^3 - 2
        sage: K.<a> = NumberField(x^3 - 2, embedding=f.roots()[0][0])
        sage: a + L(1)
        4 + 2*5^2 + 2*5^3 + 3*5^4 + 5^5 + 4*5^6 + 2*5^8 + 3*5^9 + 4*5^12
         + 4*5^14 + 4*5^15 + 3*5^16 + 5^17 + 5^18 + 2*5^19 + O(5^20)
        sage: L.<b> = NumberField(x^6 - x^2 + 1/10, embedding=1)
        sage: K.<a> = NumberField(x^3 - x + 1/10, embedding=b^2)
        sage: a + b
        b^2 + b
        sage: CC(a) == CC(b)^2
        True
        sage: K.coerce_embedding()
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10 with a = b^2
          To:   Number Field in b with defining polynomial x^6 - x^2 + 1/10
                with b = 0.9724449978911874?
          Defn: a -> b^2

    The ``QuadraticField`` and ``CyclotomicField`` constructors
    create an embedding by default unless otherwise specified::

        sage: K.<zeta> = CyclotomicField(15)
        sage: CC(zeta)
        0.913545457642601 + 0.406736643075800*I
        sage: L.<sqrtn3> = QuadraticField(-3)
        sage: K(sqrtn3)
        2*zeta^5 + 1
        sage: sqrtn3 + zeta
        2*zeta^5 + zeta + 1

    Comparison depends on the (real) embedding specified (or the one selected by default).
    Note that the codomain of the embedding must be ``QQbar`` or ``AA`` for this to work
    (see :issue:`20184`)::

        sage: N.<g> = NumberField(x^3 + 2, embedding=1)
        sage: 1 < g
        False
        sage: g > 1
        False
        sage: RR(g)
        -1.25992104989487

    If no embedding is specified or is complex, the comparison is not
    returning something meaningful. ::

        sage: N.<g> = NumberField(x^3 + 2)
        sage: 1 < g
        False
        sage: g > 1
        True

    Since SageMath 6.9, number fields may be defined by polynomials
    that are not necessarily integral or monic.  The only notable
    practical point is that in the PARI interface, a monic integral
    polynomial defining the same number field is computed and used::

        sage: K.<a> = NumberField(2*x^3 + x + 1)
        sage: K.pari_polynomial()
        x^3 - x^2 - 2

    Elements and ideals may be converted to and from PARI as follows::

        sage: pari(a)
        Mod(-1/2*y^2 + 1/2*y, y^3 - y^2 - 2)
        sage: K(pari(a))
        a
        sage: I = K.ideal(a); I
        Fractional ideal (a)
        sage: I.pari_hnf()
        [1, 0, 0; 0, 1, 0; 0, 0, 1/2]
        sage: K.ideal(I.pari_hnf())
        Fractional ideal (a)

    Here is an example where the field has non-trivial class group::

        sage: L.<b> = NumberField(3*x^2 - 1/5)
        sage: L.pari_polynomial()
        x^2 - 15
        sage: J = L.primes_above(2)[0]; J
        Fractional ideal (2, 15*b + 1)
        sage: J.pari_hnf()
        [2, 1; 0, 1]
        sage: L.ideal(J.pari_hnf())
        Fractional ideal (2, 15*b + 1)

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

        sage: W1 = NumberField(x^2 + 1,'a')
        sage: K.<x> = CyclotomicField(5)[]
        sage: W.<a> = NumberField(x^2 + 1); W
        Number Field in a with defining polynomial x^2 + 1 over its base field

    The following has been fixed in :issue:`8800`::

        sage: P.<x> = QQ[]
        sage: K.<a> = NumberField(x^3 - 5, embedding=0)
        sage: L.<b> = K.extension(x^2 + a)
        sage: F, R = L.construction()
        sage: F(R) == L    # indirect doctest
        True

    Check that :issue:`11670` has been fixed::

        sage: K.<a> = NumberField(x^2 - x - 1)
        sage: loads(dumps(K)) is K
        True
        sage: K.<a> = NumberField(x^3 - x - 1)
        sage: loads(dumps(K)) is K
        True
        sage: K.<a> = CyclotomicField(7)
        sage: loads(dumps(K)) is K
        True

    Another problem that was found while working on :issue:`11670`,
    ``maximize_at_primes`` and ``assume_disc_small`` were lost when pickling::

        sage: # needs sage.symbolic
        sage: K.<a> = NumberField(x^3 - 2, assume_disc_small=True,
        ....:                     maximize_at_primes=[2],
        ....:                     latex_name='\\alpha', embedding=2^(1/3))
        sage: L = loads(dumps(K))
        sage: L._assume_disc_small
        True
        sage: L._maximize_at_primes
        (2,)

    It is an error not to specify the generator::

        sage: K = NumberField(x^2 - 2)
        Traceback (most recent call last):
        ...
        TypeError: You must specify the name of the generator.

    Check that we can construct morphisms to matrix space (:issue:`23418`)::

        sage: t = polygen(QQ)
        sage: K = NumberField(t^4 - 2, 'a')
        sage: K.hom([K.gen().matrix()])
        Ring morphism:
          From: Number Field in a with defining polynomial x^4 - 2
          To:   Full MatrixSpace of 4 by 4 dense matrices over Rational Field
          Defn: a |--> [0 1 0 0]
                       [0 0 1 0]
                       [0 0 0 1]
                       [2 0 0 0]
    """
    if names is not None:
        name = names
    if latex_names is not None:
        latex_name = latex_names
    for key, val in kwds.items():
        if key not in ['implementation', 'prec']:
            raise TypeError("NumberField() got an unexpected keyword argument '%s'" % key)
        if not (val is None or isinstance(val, list) and all(c is None for c in val)):
            raise NotImplementedError("Number field with prescribed %s is not implemented" % key)
    if isinstance(polynomial, (list, tuple)):
        return NumberFieldTower(polynomial, names=name, check=check, embeddings=embedding, latex_names=latex_name, assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structures=structure)

    return NumberField_version2(polynomial=polynomial, name=name, check=check, embedding=embedding, latex_name=latex_name, assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structure=structure)


class NumberFieldFactory(UniqueFactory):
    r"""
    Factory for number fields.

    This should usually not be called directly, use :meth:`NumberField`
    instead.

    INPUT:

    - ``polynomial`` -- a polynomial over `\QQ` or a number field
    - ``name`` -- string (default: ``'a'``); the name of the generator
    - ``check`` -- boolean (default: ``True``); do type checking and
      irreducibility checking
    - ``embedding`` -- ``None`` or an element, the images of the generator
      in an ambient field (default: ``None``)
    - ``latex_name`` -- ``None`` or string (default: ``None``); how the
      generator is printed for latex output
    - ``assume_disc_small`` -- boolean (default: ``False``); if ``True``,
      assume that no square of a prime greater than PARI's primelimit
      (which should be 500000). Only applies for absolute fields at
      present.
    - ``maximize_at_primes`` -- ``None`` or a list of primes (default:
      ``None``); if not ``None``, then the maximal order is computed by
      maximizing only at the primes in this list, which completely avoids
      having to factor the discriminant, but of course can lead to wrong
      results; only applies for absolute fields at present.
    - ``structure`` -- ``None`` or an instance of
      :class:`structure.NumberFieldStructure` (default: ``None``),
      internally used to pass in additional structural information, e.g.,
      about the field from which this field is created as a subfield.

    TESTS::

        sage: from sage.rings.number_field.number_field import NumberFieldFactory
        sage: nff = NumberFieldFactory("number_field_factory")
        sage: R.<x> = QQ[]
        sage: nff(x^2 + 1, name='a', check=False, embedding=None, latex_name=None,
        ....:     assume_disc_small=False, maximize_at_primes=None, structure=None)
        Number Field in a with defining polynomial x^2 + 1

    Pickling preserves the ``structure()`` of a number field::

        sage: K.<a> = QuadraticField(2)
        sage: L.<b> = K.change_names()
        sage: M = loads(dumps(L))
        sage: M.structure()
        (Isomorphism given by variable name change map:
           From: Number Field in b with defining polynomial x^2 - 2
           To:   Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?,
         Isomorphism given by variable name change map:
           From: Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?
           To:   Number Field in b with defining polynomial x^2 - 2)
    """
    def create_key_and_extra_args(self, polynomial, name, check, embedding, latex_name, assume_disc_small, maximize_at_primes, structure):
        r"""
        Create a unique key for the number field specified by the parameters.

        TESTS::

            sage: from sage.rings.number_field.number_field import NumberFieldFactory
            sage: nff = NumberFieldFactory("number_field_factory")
            sage: R.<x> = QQ[]
            sage: nff.create_key_and_extra_args(x^2+1, name='a', check=False, embedding=None,
            ....:     latex_name=None, assume_disc_small=False, maximize_at_primes=None, structure=None)
            ((Rational Field, x^2 + 1, ('a',), None, 'a', None, False, None),
             {'check': False})
        """
        if name is None:
            raise TypeError("You must specify the name of the generator.")
        name = normalize_names(1, name)

        if not isinstance(polynomial, Polynomial):
            try:
                polynomial = polynomial.polynomial(QQ)
            except (AttributeError, TypeError):
                raise TypeError("polynomial (=%s) must be a polynomial." % polynomial)

        # convert polynomial to a polynomial over a field
        polynomial = polynomial.change_ring(polynomial.base_ring().fraction_field())

        # normalize embedding
        if isinstance(embedding, (list, tuple)):
            if len(embedding) != 1:
                raise TypeError("embedding must be a list of length 1")
            embedding = embedding[0]
        if embedding is not None:
            x = number_field_morphisms.root_from_approx(polynomial, embedding)
            embedding = (x.parent(), x)

        # normalize latex_name
        if isinstance(latex_name, (list, tuple)):
            if len(latex_name) != 1:
                raise TypeError("latex_name must be a list of length 1")
            latex_name = latex_name[0]

        if latex_name is None:
            latex_name = latex_variable_name(name[0])

        if maximize_at_primes is not None:
            maximize_at_primes = tuple(maximize_at_primes)

        # normalize structure
        if isinstance(structure, (list, tuple)):
            if len(structure) != 1:
                raise TypeError("structure must be a list of length 1")
            structure = structure[0]

        return (polynomial.base_ring(), polynomial, name, embedding, latex_name, maximize_at_primes, assume_disc_small, structure), {"check": check}

    def create_object(self, version, key, check):
        r"""
        Create the unique number field defined by ``key``.

        TESTS::

            sage: from sage.rings.number_field.number_field import NumberFieldFactory
            sage: nff = NumberFieldFactory("number_field_factory")
            sage: R.<x> = QQ[]
            sage: nff.create_object(None, (QQ, x^2 + 1, ('a',), None, None, None, False, None), check=False)
            Number Field in a with defining polynomial x^2 + 1
        """
        base, polynomial, name, embedding, latex_name, maximize_at_primes, assume_disc_small, structure = key

        if isinstance(base, NumberField_generic):
            from sage.rings.number_field.number_field_rel import NumberField_relative
            # Relative number fields do not support embeddings.
            return NumberField_relative(base, polynomial, name[0], latex_name,
                                        check=check, embedding=None,
                                        structure=structure)
        if polynomial.degree() == 2:
            return NumberField_quadratic(polynomial, name, latex_name, check, embedding, assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structure=structure)
        else:
            return NumberField_absolute(polynomial, name, latex_name, check, embedding, assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structure=structure)


NumberField_version2 = NumberFieldFactory("sage.rings.number_field.number_field.NumberField_version2")


def NumberFieldTower(polynomials, names, check=True, embeddings=None, latex_names=None, assume_disc_small=False, maximize_at_primes=None, structures=None):
    r"""
    Create the tower of number fields defined by the polynomials in the list
    ``polynomials``.

    INPUT:

    - ``polynomials`` -- list of polynomials. Each entry must be polynomial
      which is irreducible over the number field generated by the roots of the
      following entries.
    - ``names`` -- list of strings or a string, the names of the generators of
      the relative number fields. If a single string, then names are generated
      from that string.
    - ``check`` -- boolean (default: ``True``); whether to check that the
      polynomials are irreducible
    - ``embeddings`` -- list of elements or ``None`` (default: ``None``);
      embeddings of the relative number fields in an ambient field
    - ``latex_names`` -- list of strings or ``None`` (default: ``None``); names
      used to print the generators for latex output
    - ``assume_disc_small`` -- boolean (default: ``False``); if ``True``,
      assume that no square of a prime greater than PARI's ``primelimit``
      (which should be 500000). Only applies for absolute fields at
      present.
    - ``maximize_at_primes`` -- ``None`` or a list of primes (default:
      ``None``); if not ``None``, then the maximal order is computed by
      maximizing only at the primes in this list, which completely avoids
      having to factor the discriminant, but of course can lead to wrong
      results; only applies for absolute fields at present.
    - ``structures`` -- ``None`` or a list (default: ``None``), internally used
      to provide additional information about the number field such as the
      field from which it was created.

    OUTPUT:

    The relative number field generated by a root of the first entry of
    ``polynomials`` over the relative number field generated by root of the
    second entry of ``polynomials`` ... over the number field over which the
    last entry of ``polynomials`` is defined.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: k.<a,b,c> = NumberField([x^2 + 1, x^2 + 3, x^2 + 5]); k  # indirect doctest
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

        sage: k.absolute_field(names='c').galois_group()                                # needs sage.groups
        Galois group 8T3 (2[x]2[x]2) with order 8 of x^8 + 36*x^6 + 302*x^4 + 564*x^2 + 121

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
        sage: m = NumberField([y^3 - 3, x^2 + x + 1, y^3 + 2], 'beta'); m
        Number Field in beta0 with defining polynomial y^3 - 3 over its base field
        sage: m.base_field ()
        Number Field in beta1 with defining polynomial x^2 + x + 1 over its base field

    A tower of quadratic fields::

        sage: K.<a> = NumberField([x^2 + 3, x^2 + 2, x^2 + 1]); K
        Number Field in a0 with defining polynomial x^2 + 3 over its base field
        sage: K.base_field()
        Number Field in a1 with defining polynomial x^2 + 2 over its base field
        sage: K.base_field().base_field()
        Number Field in a2 with defining polynomial x^2 + 1

    LaTeX versions of generator names can be specified either as::

        sage: K = NumberField([x^3 - 2, x^3 - 3, x^3 - 5], names=['a', 'b', 'c'],
        ....:                 latex_names=[r'\alpha', r'\beta', r'\gamma'])
        sage: K.inject_variables(verbose=False)
        sage: latex(a + b + c)
        \alpha + \beta + \gamma

    or as::

        sage: K = NumberField([x^3 - 2, x^3 - 3, x^3 - 5], names='a', latex_names=r'\alpha')
        sage: K.inject_variables()
        Defining a0, a1, a2
        sage: latex(a0 + a1 + a2)
        \alpha_{0} + \alpha_{1} + \alpha_{2}

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
    try:
        names = normalize_names(len(polynomials), names)
    except IndexError:
        names = normalize_names(1, names)
        if len(polynomials) > 1:
            names = ['%s%s' % (names[0], i) for i in range(len(polynomials))]

    if embeddings is None:
        embeddings = [None] * len(polynomials)
    if latex_names is None:
        latex_names = [None] * len(polynomials)
    elif isinstance(latex_names, str):
        latex_names = ['%s_{%s}' % (latex_names, i) for i in range(len(polynomials))]
    if structures is None:
        structures = [None] * len(polynomials)

    if not isinstance(polynomials, (list, tuple)):
        raise TypeError("polynomials must be a list or tuple")

    if len(polynomials) == 0:
        return QQ
    if len(polynomials) == 1:
        return NumberField(polynomials[0], names=names, check=check, embedding=embeddings[0], latex_name=latex_names[0], assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structure=structures[0])

    # create the relative number field defined by f over the tower defined by polynomials[1:]
    f = polynomials[0]
    name = names[0]
    w = NumberFieldTower(polynomials[1:], names=names[1:], check=check, embeddings=embeddings[1:], latex_names=latex_names[1:], assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structures=structures[1:])
    var = f.variable_name() if isinstance(f, Polynomial) else 'x'

    R = w[var]  # polynomial ring
    return w.extension(R(f), name, check=check, embedding=embeddings[0], structure=structures[0], latex_name=latex_names[0])  # currently, extension does not accept assume_disc_small, or maximize_at_primes


def QuadraticField(D, name='a', check=True, embedding=True, latex_name='sqrt', **args):
    r"""
    Return a quadratic field obtained by adjoining a square root of
    `D` to the rational numbers, where `D` is not a
    perfect square.

    INPUT:

    - ``D`` -- a rational number

    - ``name`` -- variable name (default: ``'a'``)

    - ``check`` -- boolean (default: ``True``)

    - ``embedding`` -- boolean or square root of `D` in an ambient field
      (default: ``True``)

    - ``latex_name`` -- latex variable name (default: `\sqrt{D}`)

    OUTPUT: a number field defined by a quadratic polynomial. Unless
    otherwise specified, it has an embedding into `\RR` or
    `\CC` by sending the generator to the positive
    or upper-half-plane root.

    EXAMPLES::

        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
        sage: K.<theta> = QuadraticField(3); K
        Number Field in theta with defining polynomial x^2 - 3 with theta = 1.732050807568878?
        sage: RR(theta)
        1.73205080756888
        sage: QuadraticField(9, 'a')
        Traceback (most recent call last):
        ...
        ValueError: D must not be a perfect square.
        sage: QuadraticField(9, 'a', check=False)
        Number Field in a with defining polynomial x^2 - 9 with a = 3

    Quadratic number fields derive from general number fields.

    ::

        sage: from sage.rings.number_field.number_field_base import NumberField
        sage: type(K)
        <class 'sage.rings.number_field.number_field.NumberField_quadratic_with_category'>
        sage: isinstance(K, NumberField)
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
        sage: list(K.latex_variable_names())
        ['\\sqrt{-7}']

    We can provide our own name as well::

        sage: K.<a> = QuadraticField(next_prime(10^10), latex_name=r'\sqrt{D}')
        sage: 1 + a
        a + 1
        sage: latex(1 + a)
        \sqrt{D} + 1
        sage: latex(QuadraticField(-1, 'a', latex_name=None).gen())
        a

    The name of the generator does not interfere with Sage preparser, see :issue:`1135`::

        sage: K1 = QuadraticField(5, 'x')
        sage: K2.<x> = QuadraticField(5)
        sage: K3.<x> = QuadraticField(5, 'x')
        sage: K1 is K2
        True
        sage: K1 is K3
        True
        sage: K1
        Number Field in x with defining polynomial x^2 - 5 with x = 2.236067977499790?


    Note that, in presence of two different names for the generator,
    the name given by the preparser takes precedence::

        sage: K4.<y> = QuadraticField(5, 'x'); K4
        Number Field in y with defining polynomial x^2 - 5 with y = 2.236067977499790?
        sage: K1 == K4
        False

    TESTS::

        sage: QuadraticField(-11, 'a') is QuadraticField(-11, 'a', latex_name='Z')
        False
        sage: QuadraticField(-11, 'a') is QuadraticField(-11, 'a', latex_name=None)
        False

    Check that :issue:`23459` is fixed::

        sage: QuadraticField(4**1000+1)
        Number Field ...

    Check quadratic fields without embedding (:issue:`28932`)::

        sage: QuadraticField(3, embedding=False)
        Number Field in a with defining polynomial x^2 - 3
    """
    D = QQ(D)
    if check:
        if D.is_square():
            raise ValueError("D must not be a perfect square.")
    R = QQ['x']
    f = R([-D, 0, 1])
    if embedding is True:
        if D > 0:
            embedding = RLF(D).sqrt()
        else:
            embedding = CLF(D).sqrt()
    elif embedding is False:
        embedding = None
    if latex_name == 'sqrt':
        latex_name = r'\sqrt{%s}' % D
    return NumberField(f, name, check=False, embedding=embedding, latex_name=latex_name, **args)


def GaussianField():
    r"""
    The field `\QQ[i]`.

    TESTS::

        sage: from sage.rings.number_field.number_field import GaussianField
        sage: QQi = GaussianField()
        sage: QQi.coerce_embedding()
        Generic morphism:
          From: Number Field in I with defining polynomial x^2 + 1 with I = 1*I
          To:   Complex Lazy Field
          Defn: I -> 1*I
        sage: (I + 1/2).parent() is GaussianField()
        True
    """
    return QuadraticField(-1, 'I', latex_name='i')


def is_AbsoluteNumberField(x):
    r"""
    Return ``True`` if ``x`` is an absolute number field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import is_AbsoluteNumberField
        sage: x = polygen(ZZ, 'x')
        sage: is_AbsoluteNumberField(NumberField(x^2 + 1, 'a'))
        doctest:warning...
        DeprecationWarning: The function is_AbsoluteNumberField is deprecated; use 'isinstance(..., NumberField_absolute)' instead.
        See https://github.com/sagemath/sage/issues/38124 for details.
        True
        sage: is_AbsoluteNumberField(NumberField([x^3 + 17, x^2 + 1], 'a'))
        False

    The rationals are a number field, but they're not of the absolute
    number field class.

    ::

        sage: is_AbsoluteNumberField(QQ)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38124,
                "The function is_AbsoluteNumberField is deprecated; "
                "use 'isinstance(..., NumberField_absolute)' instead.")
    return isinstance(x, NumberField_absolute)


class CyclotomicFieldFactory(UniqueFactory):
    r"""
    Return the `n`-th cyclotomic field, where `n` is a positive integer,
    or the universal cyclotomic field if `n=0`.

    For the documentation of the universal cyclotomic field, see
    :class:`~sage.rings.universal_cyclotomic_field.UniversalCyclotomicField`.

    INPUT:

    - ``n`` -- nonnegative integer (default: `0`)

    - ``names`` -- name of generator (default: ``zetan``)

    - ``bracket`` -- defines the brackets in the case of `n=0`, and
      is ignored otherwise. Can be any even length string, with ``'()'`` being
      the default.

    - ``embedding`` -- boolean or `n`-th root of unity in an ambient field
      (default: ``True``)

    EXAMPLES:

    If called without a parameter, we get the :class:`universal cyclotomic
    field<sage.rings.universal_cyclotomic_field.UniversalCyclotomicField>`::

        sage: CyclotomicField()                                                         # needs sage.libs.gap
        Universal Cyclotomic Field

    We create the `7`\th cyclotomic field
    `\QQ(\zeta_7)` with the default generator name.

    ::

        sage: k = CyclotomicField(7); k
        Cyclotomic Field of order 7 and degree 6
        sage: k.gen()
        zeta7

    The default embedding sends the generator to the complex primitive
    `n`-th root of unity of least argument.

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
    def create_key(self, n=0, names=None, embedding=True):
        r"""
        Create the unique key for the cyclotomic field specified by the
        parameters.

        TESTS::

            sage: CyclotomicField.create_key()
            (0, None, True)
        """
        n = ZZ(n)
        if n < 0:
            raise ValueError("n (=%s) must be a positive integer" % n)
        if n > 0:
            if embedding is True:
                embedding = (CLF, (2 * CLF.pi() * CLF.gen() / n).exp())
            elif embedding is not None:
                x = number_field_morphisms.root_from_approx(QQ['x'].cyclotomic_polynomial(n), embedding)
                embedding = (x.parent(), x)
            if names is None:
                names = "zeta%s" % n
            names = normalize_names(1, names)

        return n, names, embedding

    def create_object(self, version, key, **extra_args):
        r"""
        Create the unique cyclotomic field defined by ``key``.

        TESTS::

            sage: CyclotomicField.create_object(None, (0, None, True))                  # needs sage.libs.gap
            Universal Cyclotomic Field
        """
        n, names, embedding = key
        if n == 0:
            from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
            return UniversalCyclotomicField()
        else:
            return NumberField_cyclotomic(n, names, embedding=embedding)


CyclotomicField = CyclotomicFieldFactory("sage.rings.number_field.number_field.CyclotomicField")


from . import number_field_base


is_NumberField = number_field_base.is_NumberField


class NumberField_generic(WithEqualityById, number_field_base.NumberField):
    r"""
    Generic class for number fields defined by an irreducible
    polynomial over `\QQ`.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^3 - 2); K
        Number Field in a with defining polynomial x^3 - 2
        sage: TestSuite(K).run()

    TESTS::

        sage: k.<a> = NumberField(x^3 + 2); m.<b> = NumberField(x^3 + 2)
        sage: k == QQ
        False
        sage: k.<a> = NumberField(x^3 + 2); m.<a> = NumberField(x^3 + 2)
        sage: k is m
        True
        sage: loads(dumps(k)) is k
        True

        sage: x = QQ['x'].gen()
        sage: y = ZZ['y'].gen()
        sage: K = NumberField(x^3 + x + 3, 'a'); K
        Number Field in a with defining polynomial x^3 + x + 3
        sage: K.defining_polynomial().parent()
        Univariate Polynomial Ring in x over Rational Field

        sage: L = NumberField(y^3 + y + 3, 'a'); L
        Number Field in a with defining polynomial y^3 + y + 3
        sage: L.defining_polynomial().parent()
        Univariate Polynomial Ring in y over Rational Field
        sage: L == K
        False

        sage: NumberField(ZZ['x'].0^4 + 23, 'a') == NumberField(ZZ['y'].0^4 + 23, 'a')
        False
        sage: NumberField(ZZ['x'].0^4 + 23, 'a') == NumberField(QQ['y'].0^4 + 23, 'a')
        False
        sage: NumberField(QQ['x'].0^4 + 23, 'a') == NumberField(QQ['y'].0^4 + 23, 'a')
        False

        sage: x = polygen(QQ); y = ZZ['y'].gen()
        sage: NumberField(x^3 + x + 5, 'a') == NumberField(y^3 + y + 5, 'a')
        False
        sage: NumberField(x^3 + x + 5, 'a') == NumberField(y^4 + y + 5, 'a')
        False
        sage: NumberField(x^3 + x + 5, 'a') == NumberField(x^3 + x + 5, 'b')
        False
        sage: QuadraticField(2, 'a', embedding=2) == QuadraticField(2, 'a', embedding=-2)
        False

        sage: K.<a> = QuadraticField(2)
        sage: R.<x> = K[]
        sage: L.<b> = K.extension(x^2 + 1)
        sage: M.<b> = L.absolute_field()
        sage: M == L
        False
        sage: M['x'] == L['x']
        False

        sage: R.<x> = QQ[]
        sage: R.<y> = QQ[]
        sage: K.<a> = NumberField(x^2 + 1)
        sage: L.<a> = NumberField(y^2 + 1)
        sage: K == L
        False
        sage: hash(K) == hash(L)
        False

    Two relative number fields which are isomorphic as absolute
    fields, but which are not presented the same way, are not
    considered equal (see :issue:`18942`)::

        sage: F.<omega> = NumberField(x^2 + x + 1)
        sage: y = polygen(F)
        sage: K = F.extension(y^3 + 3*omega + 2, 'alpha')
        sage: L = F.extension(y^3 - 3*omega - 1, 'alpha')
        sage: K == L
        False
        sage: K.is_isomorphic(L)
        True
        sage: hash(K) == hash(L)
        False

    This example illustrates the issue resolved in :issue:`18942`::

        sage: F.<omega> = NumberField(x^2 + x + 1)
        sage: xx = polygen(F)
        sage: ps = [p for p, _ in F(7).factor()]
        sage: for mu in ps:
        ....:     K = F.extension(xx^3 - mu, 'alpha')
        ....:     print(K.defining_polynomial().roots(K))
        [(alpha, 1), ((-omega - 1)*alpha, 1), (omega*alpha, 1)]
        [(alpha, 1), (omega*alpha, 1), ((-omega - 1)*alpha, 1)]
        sage: for mu in ps:
        ....:     K = F.extension(xx^3 - mu, 'alpha')
        ....:     print(K.defining_polynomial().roots(K))
        [(alpha, 1), ((-omega - 1)*alpha, 1), (omega*alpha, 1)]
        [(alpha, 1), (omega*alpha, 1), ((-omega - 1)*alpha, 1)]

    This example was suggested on sage-nt; see :issue:`18942`::

        sage: G = DirichletGroup(80)                                                    # needs sage.modular
        sage: for chi in G:                     # long time                             # needs sage.modular
        ....:     D = ModularSymbols(chi, 2, -1).cuspidal_subspace().new_subspace().decomposition()
        ....:     for f in D:
        ....:         elt = f.q_eigenform(10, 'alpha')[3]
        ....:         assert elt.is_integral()
    """
    def __init__(self, polynomial, name, latex_name,
                 check=True, embedding=None, category=None,
                 assume_disc_small=False, maximize_at_primes=None, structure=None):
        """
        Create a number field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: NumberField(x^97 - 19, 'a')
            Number Field in a with defining polynomial x^97 - 19

        The defining polynomial must be irreducible::

            sage: K.<a> = NumberField(x^2 - 1)
            Traceback (most recent call last):
            ...
            ValueError: defining polynomial (x^2 - 1) must be irreducible

        If you use ``check=False``, you avoid checking irreducibility of the
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
            sage: NumberField(GF(7)['x'].0^4 + 23, 'a')                                 # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: polynomial must be defined over rational field
        """
        self._assume_disc_small = assume_disc_small
        self._maximize_at_primes = maximize_at_primes
        self._structure = structure
        default_category = _NumberFields
        if category is None:
            category = default_category
        else:
            assert category.is_subcategory(default_category), "%s is not a subcategory of %s" % (category, default_category)

        Parent.__init__(self, base=QQ, names=name, category=category)
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError("polynomial (=%s) must be a polynomial" % repr(polynomial))

        if check:
            if not polynomial.parent().base_ring() == QQ:
                raise TypeError("polynomial must be defined over rational field")
            if not polynomial.is_irreducible():
                raise ValueError("defining polynomial (%s) must be irreducible" % polynomial)

        self._assign_names(name)
        self._latex_names = (latex_name,)
        self.__polynomial = polynomial
        self._pari_bnf_certified = False
        self._integral_basis_dict = {}
        if embedding is not None:
            # Since Issue #20827, an embedding is specified as a pair
            # (parent, x) with x the image of the distinguished
            # generator (previously, it was just given as x).  This
            # allows the UniqueFactory to distinguish embeddings into
            # different fields with images of the generator that
            # compare equal.
            # We allow both formats to support old pickles.
            if isinstance(embedding, tuple):
                parent, x = embedding
            else:
                parent, x = embedding.parent(), embedding
            embedding = number_field_morphisms.NumberFieldEmbedding(self, parent, x)
        self._populate_coercion_lists_(embedding=embedding, convert_method_name='_number_field_')

    def _convert_map_from_(self, other):
        r"""
        Additional conversion maps from ``other`` may be defined by
        :meth:`structure`.

        .. SEEALSO::

            :meth:`structure.NumberFieldStructure.create_structure`

        TESTS::

            sage: K.<i> = QuadraticField(-1)
            sage: L.<j> = K.change_names()
            sage: L(i)
            j
            sage: K(j)
            i

        This also works for relative number fields and their absolute fields::

            sage: K.<a> = QuadraticField(2)
            sage: x = polygen(ZZ, 'x')
            sage: L.<i> = K.extension(x^2 + 1)
            sage: M.<b> = L.absolute_field()
            sage: M(i)
            1/6*b^3 + 1/6*b
            sage: L(b)
            i - a
        """
        from sage.categories.map import Map
        if self._structure is not None:
            structure = self.structure()
            if len(structure) >= 2:
                to_self = structure[1]
                if isinstance(to_self, Map) and to_self.domain() is other:
                    return to_self
        if isinstance(other, NumberField_generic) and other._structure is not None:
            structure = other.structure()
            if len(structure) >= 1:
                from_other = structure[0]
                if isinstance(from_other, Map) and from_other.codomain() is self:
                    return from_other

    @cached_method
    def _magma_polynomial_(self, magma):
        """
        Return Magma version of the defining polynomial of this number field.

        EXAMPLES::

            sage: # optional - magma
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K._magma_polynomial_(magma)
            x^3 + 2
            sage: magma2 = Magma()
            sage: K._magma_polynomial_(magma2)
            x^3 + 2
            sage: K._magma_polynomial_(magma) is K._magma_polynomial_(magma)
            True
            sage: K._magma_polynomial_(magma) is K._magma_polynomial_(magma2)
            False
        """
        # NB f must not be garbage-collected, otherwise the
        # return value of this function is invalid
        return magma(self.defining_polynomial())

    def _magma_init_(self, magma):
        """
        Return a Magma version of this number field.

        EXAMPLES::

            sage: # optional - magma
            sage: R.<t> = QQ[]
            sage: K.<a> = NumberField(t^2 + 1)
            sage: K._magma_init_(magma)
            'SageCreateWithNames(NumberField(_sage_[...]),["a"])'
            sage: L = magma(K)
            sage: L
            Number Field with defining polynomial t^2 + 1 over the Rational Field
            sage: L.sage()
            Number Field in a with defining polynomial t^2 + 1
            sage: L.sage() is K
            True
            sage: L.1
            a
            sage: L.1^2
            -1
            sage: m = magma(a+1/2); m
            1/2*(2*a + 1)
            sage: m.sage()
            a + 1/2

        A relative number field::

            sage: # optional - magma
            sage: S.<u> = K[]
            sage: M.<b> = NumberField(u^3+u+a)
            sage: L = magma(M)
            sage: L
            Number Field with defining polynomial u^3 + u + a over its ground field
            sage: L.sage() is M
            True
        """
        # Get magma version of defining polynomial of this number field
        f = self._magma_polynomial_(magma)
        s = 'NumberField(%s)' % f.name()
        return magma._with_names(s, self.variable_names())

    def construction(self):
        r"""
        Construction of ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 + x^2 + 1, embedding=CC.gen())
            sage: F, R = K.construction()
            sage: F
            AlgebraicExtensionFunctor
            sage: R
            Rational Field

        The construction functor respects distinguished embeddings::

            sage: F(R) is K
            True
            sage: F.embeddings
            [0.2327856159383841? + 0.7925519925154479?*I]

        TESTS::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: R.<t> = ZZ[]
            sage: a + t     # indirect doctest
            t + a
            sage: (a + t).parent()
            Univariate Polynomial Ring in t over Number Field in a with defining polynomial x^3 + x + 1

        The construction works for non-absolute number fields as well::

            sage: K.<a,b,c>=NumberField([x^3+x^2+1,x^2+1,x^7+x+1])
            sage: F,R = K.construction()
            sage: F(R) == K
            True

        ::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-5, embedding=0)
            sage: L.<b> = K.extension(x^2+a)
            sage: a*b
            a*b

        The construction preserves latex variable names::

            sage: K.<a,b,c> = NumberField([x^3+x^2+1, x^2+1, x^7+x+1], latex_names=['alpha', 'beta', 'gamma'])
            sage: F, R = K.construction()
            sage: F(R) == K
            True
        """
        from sage.categories.pushout import AlgebraicExtensionFunctor
        from sage.rings.rational_field import QQ
        names = self.variable_names()
        polys = []
        embeddings = []
        structures = []
        latex_names = []
        K = self
        while K is not QQ:
            polys.append(K.relative_polynomial())
            embeddings.append(None if K.coerce_embedding() is None else K.coerce_embedding()(K.gen()))
            structures.append(K._structure)
            latex_names.append(K.latex_variable_names()[0])
            K = K.base_field()
        return (AlgebraicExtensionFunctor(polys, names, embeddings, structures,
                                          latex_names=latex_names), QQ)

    def _element_constructor_(self, x, check=True):
        r"""
        Convert ``x`` into an element of this number field.

        INPUT:

        - ``x`` -- Sage (or Python) object

        OUTPUT:

        A :class:`~number_field_element.NumberFieldElement`
        constructed from ``x``.

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 + 17)
            sage: K(a) is a  # indirect doctest
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

        We can create number field elements from PARI::

            sage: K.<a> = NumberField(x^3 - 17)
            sage: K(pari(42))                                                           # needs sage.libs.pari
            42
            sage: K(pari("5/3"))                                                        # needs sage.libs.pari
            5/3
            sage: K(pari("[3/2, -5, 0]~"))    # Uses Z-basis                            # needs sage.libs.pari
            -5/3*a^2 + 5/3*a - 1/6

        From a PARI polynomial or ``POLMOD``, note that the variable
        name does not matter::

            sage: K(pari("-5/3*q^2 + 5/3*q - 1/6"))                                     # needs sage.libs.pari
            -5/3*a^2 + 5/3*a - 1/6
            sage: K(pari("Mod(-5/3*q^2 + 5/3*q - 1/6, q^3 - 17)"))                      # needs sage.libs.pari
            -5/3*a^2 + 5/3*a - 1/6
            sage: K(pari("x^5/17"))                                                     # needs sage.libs.pari
            a^2

        An error is raised when a PARI element with an incorrect
        modulus is given::

            sage: K(pari("Mod(-5/3*q^2 + 5/3*q - 1/6, q^3 - 999)"))                     # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            TypeError: cannot convert PARI element Mod(-5/3*q^2 + 5/3*q - 1/6, q^3 - 999)
            into Number Field in a with defining polynomial x^3 - 17

        Test round-trip conversion to PARI and back::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 - 1/2*x + 1/3)
            sage: b = K.random_element()
            sage: K(pari(b)) == b                                                       # needs sage.libs.pari
            True

            sage: F.<c> = NumberField(2*x^3 + x + 1)
            sage: d = F.random_element()
            sage: F(F.pari_nf().nfalgtobasis(d)) == d                                   # needs sage.libs.pari
            True

        If the PARI polynomial is different from the Sage polynomial,
        a warning is printed unless ``check=False`` is specified::

            sage: b = pari(a); b                                                        # needs sage.libs.pari
            Mod(-1/12*y^2 - 1/12*y + 1/6, y^3 - 3*y - 22)
            sage: K(b.lift())                                                           # needs sage.libs.pari
            doctest:...: UserWarning: interpreting PARI polynomial -1/12*y^2 - 1/12*y + 1/6
            relative to the defining polynomial x^3 - 3*x - 22 of the PARI number field
            a
            sage: K(b.lift(), check=False)                                              # needs sage.libs.pari
            a

        Using a GAP element may be tricky, as it may contain
        an exclamation mark::

            sage: L.<tau> = NumberField(x^3 - 2)
            sage: gap(tau^3)                                                            # needs sage.libs.gap
            2
            sage: gap(tau)^3                                                            # needs sage.libs.gap
            !2
            sage: L(gap(tau)^3)     # indirect doctest                                  # needs sage.libs.gap
            2

        Check that :issue:`22202` and :issue:`27765` are fixed::

            sage: y = QQ['y'].gen()
            sage: R = QQ.extension(y^2 - 1/2, 'a')['x']
            sage: R("a*x").factor()
            (a) * x

        Check that :issue:`30961` is fixed::

            sage: # needs sage.symbolic
            sage: QQi = i.parent()
            sage: x = SR.var('x')
            sage: QQi((x, x))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a rational
            sage: QQi(("1", "2"))
            2*I + 1
            sage: QQi((RR(1), RR(2)))
            2*I + 1
            sage: QQi(vector((RR(1), RR(2))))
            2*I + 1

        Check that :issue:`34059` is fixed::

            sage: K.<a> = NumberField(x)
            sage: K([1]).parent()
            Number Field in a with defining polynomial x
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            K = x.parent()
            if K is self:
                return x
            elif isinstance(x, (number_field_element.OrderElement_absolute,
                                number_field_element.OrderElement_relative,
                                number_field_element_quadratic.OrderElement_quadratic)):
                L = K.number_field()
                if L is self:
                    return self._element_class(self, x)
                x = L(x)
            return self._coerce_from_other_number_field(x)
        elif isinstance(x, pari_gen):
            if x.type() == "t_POLMOD":
                modulus = x.mod()
                if check and modulus != self.pari_polynomial(modulus.variable()):
                    raise TypeError("cannot convert PARI element %s into %s" % (x, self))
                x = x.lift()
                check = False
            elif x.type() == "t_COL":
                x = self.pari_nf().nfbasistoalg_lift(x)
                check = False
            if x.type() in ["t_INT", "t_FRAC"]:
                pass
            elif x.type() == "t_POL":
                # We consider x as a polynomial in the standard
                # generator of the PARI number field, and convert it
                # to a polynomial in the Sage generator.
                if any(x.poldegree(v) > 0 for v in x.variables()):
                    var = self.absolute_polynomial().variable_name()
                    if check and self.pari_polynomial(var) != self.absolute_polynomial().monic():
                        from warnings import warn
                        warn("interpreting PARI polynomial %s relative to the defining polynomial %s of the PARI number field"
                             % (x, self.pari_polynomial()))
                    beta = self._pari_absolute_structure()[2]
                    x = x(beta).lift()
                else:  # constant polynomial
                    x = x[0]
            else:
                raise TypeError("%s has unsupported PARI type %s" % (x, x.type()))
            x = self.absolute_polynomial().parent()(x)
            return self._element_class(self, x)
        elif isinstance(x, GapElement):
            s = x._sage_repr()
            if self.variable_name() in s:
                return self._convert_from_str(s)
            return self._convert_from_str(s.replace('!', ''))
        elif isinstance(x, str):
            return self._convert_from_str(x)
        elif isinstance(x, (tuple, list,
                            sage.modules.free_module_element.FreeModuleElement)):
            if len(x) != self.relative_degree():
                raise ValueError("Length must be equal to the degree of this number field")
            base = self.base_ring()
            return sum(base(c) * g for c, g in zip(x, self.gen(0).powers(len(x))))
        return self._convert_non_number_field_element(x)

    def _convert_non_number_field_element(self, x):
        """
        Convert a non-number field element ``x`` into this number field.

        INPUT:

        - ``x`` -- a non number field element, e.g., a list, integer,
          rational, or polynomial

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 + 2/3)
            sage: K._convert_non_number_field_element(-7/8)
            -7/8
            sage: K._convert_non_number_field_element([1,2,3])
            3*a^2 + 2*a + 1

        The list is just turned into a polynomial in the generator::

            sage: K._convert_non_number_field_element([0,0,0,1,1])
            -2/3*a - 2/3

        Any polynomial whose coefficients can be converted to rationals
        will convert to the number field, e.g., this one in
        characteristic 7::

            sage: f = GF(7)['y']([1,2,3]); f
            3*y^2 + 2*y + 1
            sage: K._convert_non_number_field_element(f)
            3*a^2 + 2*a + 1

        But not this one over a field of order 27::

            sage: F27.<g> = GF(27)                                                      # needs sage.rings.finite_rings
            sage: f = F27['z']([g^2, 2*g, 1]); f                                        # needs sage.rings.finite_rings
            z^2 + 2*g*z + g^2
            sage: K._convert_non_number_field_element(f)                                # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: unable to convert g^2 to a rational

        One can also convert an element of the polynomial quotient ring
        that is isomorphic to the number field::

            sage: K.<a> = NumberField(x^3 + 17)
            sage: b = K.polynomial_quotient_ring().random_element()
            sage: b == K(b)
            True

        We can convert symbolic expressions::

            sage: # needs sage.symbolic
            sage: I = sqrt(-1); parent(I)
            Symbolic Ring
            sage: GaussianIntegers()(2 + I)
            I + 2
            sage: K1 = QuadraticField(3)
            sage: K2 = QuadraticField(5)
            sage: (K,) = K1.composite_fields(K2, preserve_embedding=True)
            sage: K(sqrt(3) + sqrt(5))
            -1/2*a0^3 + 8*a0
            sage: K(sqrt(-3)*I)
            1/4*a0^3 - 7/2*a0
        """
        if isinstance(x, (int, Rational, Integer, pari_gen, list)):
            return self._element_class(self, x)

        if isinstance(x, sage.rings.polynomial.polynomial_quotient_ring_element.PolynomialQuotientRingElement)\
           and (x in self.polynomial_quotient_ring()):
            y = self.polynomial_ring().gen()
            return x.lift().subs({y: self.gen()})

        if isinstance(x, (sage.rings.qqbar.AlgebraicNumber, sage.rings.qqbar.AlgebraicReal)):
            return self._convert_from_qqbar(x)

        if isinstance(x, polynomial_element.Polynomial):
            return self._element_class(self, x)

        # Try converting via QQ.
        try:
            y = QQ(x)
        except (TypeError, ValueError):
            pass
        else:
            return self._element_class(self, y)

        # Final attempt: convert via QQbar. This deals in particular
        # with symbolic expressions like sqrt(-5).
        try:
            y = sage.rings.qqbar.QQbar(x)
        except (TypeError, ValueError):
            pass
        else:
            return self._convert_from_qqbar(y)

        raise TypeError("unable to convert %r to %s" % (x, self))

    def _convert_from_qqbar(self, x):
        """
        Convert an element of ``QQbar`` or ``AA`` to this number field,
        if possible.

        This requires that the given number field is equipped with an
        embedding.

        INPUT:

        - ``x`` -- an algebraic number

        EXAMPLES::

            sage: K.<a> = QuadraticField(3)
            sage: K._convert_from_qqbar(7 + 2*AA(3).sqrt())
            2*a + 7
            sage: GaussianIntegers()(QQbar(I))
            I
            sage: CyclotomicField(15)(QQbar.zeta(5))
            zeta15^3
            sage: CyclotomicField(12)(QQbar.zeta(5))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 0.3090169943749474? + 0.9510565162951536?*I
            to Cyclotomic Field of order 12 and degree 4
        """
        # We use the diagram
        #
        # self
        #  ↑  ↘
        #  F → QQbar
        #
        # Where F is the smallest number field containing x.
        #
        # y is the pre-image such that x = F(y)
        F, y, F_to_QQbar = x.as_number_field_element(minimal=True)

        # Try all embeddings from F into self
        from sage.rings.qqbar import QQbar
        for F_to_self in F.embeddings(self):
            z = F_to_self(y)
            # Check whether the diagram commutes
            if QQbar(z) == x:
                return z

        raise TypeError("unable to convert %r to %s" % (x, self))

    def _convert_from_str(self, x):
        """
        Coerce a string representation of an element of this
        number field into this number field.

        INPUT:

        - ``x`` -- string

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<theta25> = NumberField(x^3 + (2/3)*x + 1)
            sage: k._convert_from_str('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1

        This function is called by the coerce method when it gets a string
        as input::

            sage: k('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1
        """
        w = sage_eval(x, locals=self.gens_dict())
        if not (isinstance(w, Element) and w.parent() is self):
            return self(w)
        else:
            return w

    def _Hom_(self, codomain, category=None):
        """
        Return homset of homomorphisms from ``self`` to the number field codomain.

        EXAMPLES:

        This method is implicitly called by :meth:`Hom` and
        :meth:`sage.categories.homset.Hom`::

            sage: x = polygen(ZZ, 'x')
            sage: K.<i> = NumberField(x^2 + 1); K
            Number Field in i with defining polynomial x^2 + 1
            sage: K.Hom(K)  # indirect doctest
            Automorphism group of Number Field in i with defining polynomial x^2 + 1
            sage: Hom(K, QuadraticField(-1, 'b'))
            Set of field embeddings
             from Number Field in i with defining polynomial x^2 + 1
               to Number Field in b with defining polynomial x^2 + 1 with b = 1*I

        CHECKME: handling of the case where codomain is not a number field?

        ::

            sage: Hom(K, VectorSpace(QQ,3))
            Set of Morphisms
             from Number Field in i with defining polynomial x^2 + 1
               to Vector space of dimension 3 over Rational Field
               in Category of commutative additive groups

        TESTS:

        Verify that :issue:`22001` has been resolved::

            sage: R.<x> = QQ[]
            sage: K.<a> = QQ.extension(x^2 + 1)
            sage: K.hom([a]).category_for()
            Category of number fields

        ::

            sage: H = End(K)
            sage: loads(dumps(H)) is H
            True
        """
        if not is_NumberFieldHomsetCodomain(codomain):
            # Using LazyFormat fixes #28036 - infinite loop
            from sage.misc.lazy_format import LazyFormat
            raise TypeError(LazyFormat("%s is not suitable as codomain for homomorphisms from %s") % (codomain, self))
        from sage.rings.number_field.homset import NumberFieldHomset
        return NumberFieldHomset(self, codomain, category)

    @cached_method
    def structure(self):
        """
        Return fixed isomorphism or embedding structure on ``self``.

        This is used to record various isomorphisms or embeddings that
        arise naturally in other constructions.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
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
        if self._structure is None:
            f = self.hom(self)
            return f, f
        else:
            return self._structure.create_structure(self)

    def completion(self, p, prec, extras={}):
        """
        Return the completion of ``self`` at `p` to the specified precision.

        Only implemented at archimedean places, and then only if
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
            raise ValueError("No embedding into the complex numbers has been specified.")
        else:
            raise NotImplementedError

    def primitive_element(self):
        r"""
        Return a primitive element for this field, i.e., an element that
        generates it over `\QQ`.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.primitive_element()
            a
            sage: K.<a,b,c> = NumberField([x^2 - 2, x^2 - 3, x^2 - 5])
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

        - ``num_bound`` -- bound on numerator of the coefficients of
          the resulting element

        - ``den_bound`` -- bound on denominators of the coefficients
          of the resulting element

        - ``integral_coefficients`` -- boolean (default: ``False``); if
          ``True``, then the resulting element will have integral coefficients.
          This option overrides any value of ``den_bound``.

        - ``distribution`` -- distribution to use for the coefficients
          of the resulting element

        OUTPUT: element of this number field

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<j> = NumberField(x^8 + 1)
            sage: K.random_element().parent() is K
            True

            sage: while K.random_element().list()[0] != 0:
            ....:     pass
            sage: while K.random_element().list()[0] == 0:
            ....:     pass
            sage: while K.random_element().is_prime():
            ....:     pass
            sage: while not K.random_element().is_prime():
            ....:     pass

            sage: K.<a,b,c> = NumberField([x^2 - 2, x^2 - 3, x^2 - 5])
            sage: K.random_element().parent() is K
            True

            sage: while K.random_element().is_prime():
            ....:     pass
            sage: while not K.random_element().is_prime():  # long time
            ....:     pass

            sage: K.<a> = NumberField(x^5 - 2)
            sage: p = K.random_element(integral_coefficients=True)
            sage: p.is_integral()
            True
            sage: while K.random_element().is_integral():
            ....:     pass

        TESTS::

            sage: K.<a> = NumberField(x^5 - 2)
            sage: K.random_element(-1)
            Traceback (most recent call last):
            ...
            TypeError: x must be < y
            sage: K.random_element(5, 0)
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
        is a relative extension `M/L`) and a map from `K` to ``self`` that
        sends the generator of `K` to ``alpha``.

        INPUT:

        - ``alpha`` -- an element of ``self``, or something that
          coerces to an element of ``self``

        OUTPUT:

        - ``K`` -- a number field
        - ``from_K`` -- a homomorphism from `K` to ``self`` that
          sends the generator of `K` to ``alpha``

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^4 - 3); K
            Number Field in a with defining polynomial x^4 - 3
            sage: H.<b>, from_H = K.subfield(a^2)
            sage: H
            Number Field in b with defining polynomial x^2 - 3 with b = a^2
            sage: from_H(b)
            a^2
            sage: from_H
            Ring morphism:
              From: Number Field in b with defining polynomial x^2 - 3 with b = a^2
              To:   Number Field in a with defining polynomial x^4 - 3
              Defn: b |--> a^2

        A relative example. Note that the result returned is the subfield generated
        by `\alpha` over ``self.base_field()``, not over `\QQ` (see :issue:`5392`)::

            sage: L.<a> = NumberField(x^2 - 3)
            sage: M.<b> = L.extension(x^4 + 1)
            sage: K, phi = M.subfield(b^2)
            sage: K.base_field() is L
            True

        Subfields inherit embeddings::

            sage: K.<z> = CyclotomicField(5)
            sage: L, K_from_L = K.subfield(z - z^2 - z^3 + z^4)
            sage: L
            Number Field in z0 with defining polynomial x^2 - 5 with z0 = 2.236067977499790?
            sage: CLF_from_K = K.coerce_embedding(); CLF_from_K
            Generic morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Lazy Field
              Defn: z -> 0.309016994374948? + 0.951056516295154?*I
            sage: CLF_from_L = L.coerce_embedding(); CLF_from_L
            Generic morphism:
              From: Number Field in z0 with defining polynomial x^2 - 5
                    with z0 = 2.236067977499790?
              To:   Complex Lazy Field
              Defn: z0 -> 2.236067977499790?

        Check transitivity::

            sage: CLF_from_L(L.gen())
            2.236067977499790?
            sage: CLF_from_K(K_from_L(L.gen()))
            2.23606797749979? + 0.?e-14*I

        If ``self`` has no specified embedding, then `K` comes with an
        embedding in ``self``::

            sage: K.<a> = NumberField(x^6 - 6*x^4 + 8*x^2 - 1)
            sage: L.<b>, from_L = K.subfield(a^2)
            sage: L
            Number Field in b with defining polynomial x^3 - 6*x^2 + 8*x - 1 with b = a^2
            sage: L.gen_embedding()
            a^2

        You can also view a number field as having a different generator by
        just choosing the input to generate the whole field; for that it is
        better to use :meth:`change_generator`, which gives
        isomorphisms in both directions.
        """
        if names is not None:
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
        Given the number field ``self``, construct another isomorphic number
        field `K` generated by the element ``alpha`` of ``self``, along
        with isomorphisms from `K` to ``self`` and from ``self`` to
        `K`.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<i> = NumberField(x^2 + 1); L
            Number Field in i with defining polynomial x^2 + 1
            sage: K, from_K, to_K = L.change_generator(i/2 + 3)
            sage: K
            Number Field in i0 with defining polynomial x^2 - 6*x + 37/4 with i0 = 1/2*i + 3
            sage: from_K
            Ring morphism:
              From: Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
                    with i0 = 1/2*i + 3
              To:   Number Field in i with defining polynomial x^2 + 1
              Defn: i0 |--> 1/2*i + 3
            sage: to_K
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
                    with i0 = 1/2*i + 3
              Defn: i |--> 2*i0 - 6

        We can also do

        ::

            sage: K.<c>, from_K, to_K = L.change_generator(i/2 + 3); K
            Number Field in c with defining polynomial x^2 - 6*x + 37/4 with c = 1/2*i + 3


        We compute the image of the generator `\sqrt{-1}` of `L`.

        ::

            sage: to_K(i)
            2*c - 6

        Note that the image is indeed a square root of `-1`.

        ::

            sage: to_K(i)^2
            -1
            sage: from_K(to_K(i))
            i
            sage: to_K(from_K(c))
            c
        """
        if names is not None:
            name = names
        alpha = self(alpha)
        K, from_K = self.subfield(alpha, name=name)
        if K.degree() != self.degree():
            raise ValueError("alpha must generate a field of degree %s, but alpha generates a subfield of degree %s" % (self.degree(), K.degree()))
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

    def subfield_from_elements(self, alpha, name=None, polred=True, threshold=None):
        r"""
        Return the subfield generated by the elements ``alpha``.

        If the generated subfield by the elements ``alpha`` is either the
        rational field or the complete number field, the field returned is
        respectively ``QQ`` or ``self``.

        INPUT:

        - ``alpha`` -- list of elements in this number field

        - ``name`` -- a name for the generator of the new number field

        - ``polred`` -- boolean (default: ``True``); whether to optimize the
          generator of the newly created field

        - ``threshold`` -- positive number (default: ``None``) threshold to be
          passed to the ``do_polred`` function

        OUTPUT: a triple ``(field, beta, hom)`` where

        - ``field`` -- a subfield of this number field

        - ``beta`` -- list of elements of ``field`` corresponding to ``alpha``

        - ``hom`` -- inclusion homomorphism from ``field`` to ``self``

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: poly = x^4 - 4*x^2 + 1
            sage: emb = AA.polynomial_root(poly, RIF(0.51, 0.52))
            sage: K.<a> = NumberField(poly, embedding=emb)
            sage: sqrt2 = -a^3 + 3*a
            sage: sqrt3 = -a^2 + 2
            sage: assert sqrt2 ** 2 == 2 and sqrt3 ** 2 == 3
            sage: L, elts, phi = K.subfield_from_elements([sqrt2, 1 - sqrt2/3])
            sage: L
            Number Field in a0 with defining polynomial x^2 - 2 with a0 = 1.414213562373095?
            sage: elts
            [a0, -1/3*a0 + 1]
            sage: phi
            Ring morphism:
              From: Number Field in a0 with defining polynomial x^2 - 2
                    with a0 = 1.414213562373095?
              To:   Number Field in a with defining polynomial x^4 - 4*x^2 + 1
                    with a = 0.5176380902050415?
              Defn: a0 |--> -a^3 + 3*a
            sage: assert phi(elts[0]) == sqrt2
            sage: assert phi(elts[1]) == 1 - sqrt2/3


            sage: L, elts, phi = K.subfield_from_elements([1, sqrt3])
            sage: assert phi(elts[0]) == 1
            sage: assert phi(elts[1]) == sqrt3

            sage: L, elts, phi = K.subfield_from_elements([sqrt2, sqrt3])
            sage: phi
            Identity endomorphism of Number Field in a with defining polynomial
             x^4 - 4*x^2 + 1 with a = 0.5176380902050415?

        TESTS::

            sage: x = polygen(QQ)

            sage: p = x^8 - 12*x^6 + 23*x^4 - 12*x^2 + 1
            sage: K.<a> = NumberField(p)
            sage: sqrt2 = 6/7*a^7 - 71/7*a^5 + 125/7*a^3 - 43/7*a
            sage: sqrt3 = 3/7*a^6 - 32/7*a^4 + 24/7*a^2 + 10/7
            sage: sqrt5 = 8/7*a^6 - 90/7*a^4 + 120/7*a^2 - 27/7
            sage: assert sqrt2**2 == 2 and sqrt3**2 == 3 and sqrt5**2 == 5
            sage: L, elts, phi = K.subfield_from_elements([sqrt2, sqrt3], name='b')
            sage: assert phi(elts[0]) == sqrt2
            sage: assert phi(elts[1]) == sqrt3
            sage: L, elts, phi = K.subfield_from_elements([sqrt2, sqrt3], name='b', polred=False)
            sage: assert phi(elts[0]) == sqrt2
            sage: assert phi(elts[1]) == sqrt3
            sage: L, elts, phi = K.subfield_from_elements([sqrt2, sqrt5], name='b')
            sage: assert phi(elts[0]) == sqrt2
            sage: assert phi(elts[1]) == sqrt5
            sage: L, elts, phi = K.subfield_from_elements([sqrt3, sqrt5], name='b')
            sage: assert phi(elts[0]) == sqrt3
            sage: assert phi(elts[1]) == sqrt5
            sage: L, elts, phi = K.subfield_from_elements([-149582/214245 + 21423/5581*sqrt2], name='b')
            sage: assert L.polynomial() == x^2 - 2
            sage: L, elts, phi = K.subfield_from_elements([131490/777 - 1359/22*sqrt3], name='b')
            sage: assert L.polynomial() == x^2 - 3
            sage: L, elts, phi = K.subfield_from_elements([12241829/177 - 321121/22459 * sqrt5], name='b')
            sage: assert L.polynomial() == x^2 - x - 1

            sage: from sage.rings.qqbar import number_field_elements_from_algebraics
            sage: R.<x> = QQ[]
            sage: p1 = x^3 - x - 1
            sage: roots1 = p1.roots(QQbar, False)
            sage: for _ in range(10):
            ....:     p2 = R.random_element(degree=2)
            ....:     while not p2.is_irreducible(): p2 = R.random_element(degree=2)
            ....:     roots2 = p2.roots(QQbar, False)
            ....:     K, (a1,b1,c1,a2,b2), phi = number_field_elements_from_algebraics(roots1 + roots2)
            ....:     u1 = 1 - a1/17 + 3/7*a1**2
            ....:     u2 = 2 + 33/35 * a1
            ....:     L, (v1,v2), phi = K.subfield_from_elements([u1, u2], threshold=100)
            ....:     assert L.polynomial() == p1
            ....:     assert phi(v1) == u1 and phi(v2) == u2
        """
        alpha = [self(a) for a in alpha]

        # Rational case
        if all(a.is_rational() for a in alpha):
            return (QQ, [QQ(a) for a in alpha], self.coerce_map_from(QQ))

        # Saturate with multiplication
        from sage.modules.free_module import VectorSpace
        V = VectorSpace(QQ, self.degree())
        vecs = [a.vector() for a in alpha]
        U = V.subspace(vecs)
        modified = True
        while modified:
            modified = False
            d = U.dimension()
            if d == self.degree():
                # the given elements generate the whole number field
                return (self, alpha, self.hom(self))
            B = U.basis()
            for i in range(d):
                for j in range(i, d):
                    v = (self(B[i]) * self(B[j])).vector()
                    if v not in U:
                        U += V.subspace([v])
                        modified = True

        # strict subfield, find a generator
        vgen = None
        for b in U.basis():
            if self(b).minpoly().degree() == d:
                vgen = b
                break
        if vgen is None:
            s = 1
            while True:
                vgen = U.random_element(proba=0.5, x=-s, y=s)
                if self(vgen).minpoly().degree() == d:
                    break
                s *= 2

        # minimize the generator via PARI polred
        gen = self(vgen)
        p = gen.minpoly()
        if polred:
            from sage.rings.qqbar import do_polred
            if threshold:
                fwd, _, q = do_polred(p, threshold)
            else:
                fwd, _, q = do_polred(p)
        else:
            q = p
            fwd = self.polynomial_ring().gen()

        new_gen = fwd(gen)
        assert new_gen.minpoly() == q
        K, hom = self.subfield(new_gen, name=name)

        # express the elements in the basis 1, new_gen, new_gen^2, ..., new_gen^(deg-1)
        from sage.matrix.constructor import matrix
        M = matrix(QQ, [(new_gen**i).vector() for i in range(d)])
        new_alpha = [K(M.solve_left(elt.vector())) for elt in alpha]

        return (K, new_alpha, hom)

    def is_absolute(self):
        """
        Return ``True`` if ``self`` is an absolute field.

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

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^10 - 2)
            sage: K.is_absolute()
            True
            sage: K.is_relative()
            False
        """
        return not self.is_absolute()

    def quadratic_defect(self, a, p, check=True):
        r"""
        Return the valuation of the quadratic defect of `a` at `p`.

        INPUT:

        - ``a`` -- an element of ``self``

        - ``p`` -- a prime ideal

        - ``check`` -- boolean (default: ``True``); check if `p` is prime

        ALGORITHM:

        This is an implementation of Algorithm 3.1.3 from [Kir2016]_.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 2)
            sage: p = K.primes_above(2)[0]
            sage: K.quadratic_defect(5, p)
            4
            sage: K.quadratic_defect(0, p)
            +Infinity
            sage: K.quadratic_defect(a, p)
            1
            sage: K.<a> = CyclotomicField(5)
            sage: p = K.primes_above(2)[0]
            sage: K.quadratic_defect(5, p)
            +Infinity

        TESTS::

            sage: x = polygen(QQ, 'x')
            sage: K.<t> = NumberField(x^10 - x^8 - 2*x^7 - x^6 + 2*x^5 + 2*x^4 - x^2 + 1)
            sage: p = K.prime_factors(2)[0]
            sage: pi = K.uniformizer(p)
            sage: a = 1 + pi^3
            sage: K.quadratic_defect(a, p)
            3
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if a not in self:
            raise TypeError(str(a) + " must be an element of " + str(self))
        if not self == QQ and not p.parent() == self.ideal_monoid():
            raise TypeError(str(p) + " is not a prime ideal in "
                            + str(self.ideal_monoid()))
        if check and not p.is_prime():
            raise ValueError(str(p) + " must be prime")
        if a.is_zero():
            return Infinity
        v = self(a).valuation(p)
        if v % 2:
            return v
        pi = self.uniformizer(p)
        a = self(a) / pi**v
        F = p.residue_field()
        q = F.reduction_map()
        # The non-dyadic case
        if self(2).valuation(p) == 0:
            if q(a).is_square():
                return Infinity
            return v
        # The dyadic case
        s = self(F.lift((1 / F(a)).sqrt()))
        a = self(s**2) * a
        u = self(4).valuation(p)
        w = (a - 1).valuation(p)
        x = PolynomialRing(F, 'x').gen()
        f = x**2 + x
        while w < u and not w % 2:
            s = F.lift(q((a - 1) / pi**w).sqrt())
            a = a / (1 + s*(pi**(w/2)))**2
            w = (a - 1).valuation(p)
        if w < u and w % 2:
            return v + w
        if w == u and (f + F((a-1) / 4)).is_irreducible():
            return v + w
        return Infinity

    def absolute_field(self, names):
        """
        Return ``self`` as an absolute number field.

        INPUT:

        - ``names`` -- string; name of generator of the absolute field

        OUTPUT:

        - ``K`` -- this number field (since it is already absolute)

        Also, ``K.structure()`` returns ``from_K`` and ``to_K``, where
        ``from_K`` is an isomorphism from `K` to ``self`` and ``to_K``
        is an isomorphism from ``self`` to `K`.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.absolute_field('a')
            Number Field in a with defining polynomial x^4 + x^3 + x^2 + x + 1
        """
        return NumberField(self.defining_polynomial(), names, check=False, structure=structure.NameChange(self))

    def is_isomorphic(self, other, isomorphism_maps=False) -> bool:
        """
        Return ``True`` if ``self`` is isomorphic as a number field to ``other``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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

            sage: k = NumberField(x^2 - x - 1, 'b')
            sage: l = NumberField(x^2 - 7, 'a')
            sage: k.is_isomorphic(l, True)
            (False, [])

            sage: k = NumberField(x^2 - x - 1, 'b')
            sage: ky.<y> = k[]
            sage: l = NumberField(y, 'a')
            sage: k.is_isomorphic(l, True)
            (True, [-x, x + 1])

        TESTS:

        See :issue:`26239`::

            sage: k.<a> = NumberField(x)
            sage: k.is_isomorphic(k)
            True
        """
        if not isinstance(other, NumberField_generic):
            raise ValueError("other must be a generic number field.")
        t = self.pari_polynomial().nfisisom(other.pari_polynomial())
        # NB t==0 returns True when t is [0]
        if t.length() == 0:
            t = []
            res = False
        else:
            res = True

        if isomorphism_maps:
            R = self.polynomial().parent()
            return res, [R(ti) for ti in t]
        else:
            return res

    def is_totally_real(self):
        """
        Return ``True`` if ``self`` is totally real, and ``False`` otherwise.

        Totally real means that every isomorphic embedding of ``self`` into the
        complex numbers has image contained in the real numbers.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 2, 'alpha').is_totally_real()
            False
            sage: NumberField(x^2 - 2, 'alpha').is_totally_real()
            True
            sage: NumberField(x^4 - 2, 'alpha').is_totally_real()
            False
        """
        return self.signature()[1] == 0

    def is_totally_imaginary(self):
        """
        Return ``True`` if ``self`` is totally imaginary, and ``False`` otherwise.

        Totally imaginary means that no isomorphic embedding of ``self`` into
        the complex numbers has image contained in the real numbers.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 2, 'alpha').is_totally_imaginary()
            True
            sage: NumberField(x^2 - 2, 'alpha').is_totally_imaginary()
            False
            sage: NumberField(x^4 - 2, 'alpha').is_totally_imaginary()
            False
        """
        return self.signature()[0] == 0

    def is_CM(self):
        r"""
        Return ``True`` if ``self`` is a CM field (i.e., a totally imaginary
        quadratic extension of a totally real field).

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
            sage: F.<a> = NumberField(x^4 - x^3 - 3*x^2 + x + 1)
            sage: F.is_CM()
            False

        The following are non-CM totally imaginary fields.

        ::

            sage: F.<a> = NumberField(x^4 + x^3 - x^2 - x + 1)
            sage: F.is_totally_imaginary()
            True
            sage: F.is_CM()
            False
            sage: F2.<a> = NumberField(x^12 - 5*x^11 + 8*x^10 - 5*x^9 - x^8 + 9*x^7
            ....:                       + 7*x^6 - 3*x^5 + 5*x^4 + 7*x^3 - 4*x^2 - 7*x + 7)
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

            sage: E_0.<a> = NumberField(x^7 - 4*x^6 - 4*x^5 + 10*x^4 + 4*x^3
            ....:                       - 6*x^2 - x + 1)
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
        # Return cached answer if available
        try:
            return self.__is_CM
        except (AttributeError):
            pass

        # Then, deal with simple cases
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

        # Check for index 2 subextensions that are totally real
        possibilities = K.subfields(K.absolute_degree() / 2)
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
        Return the complex conjugation of ``self``.

        This is only well-defined for fields contained in CM fields
        (i.e. for totally real fields and CM fields). Recall that a CM
        field is a totally imaginary quadratic extension of a totally
        real field. For other fields, a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: QuadraticField(-1, 'I').complex_conjugation()
            Ring endomorphism of
             Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              Defn: I |--> -I
            sage: CyclotomicField(8).complex_conjugation()
            Ring endomorphism of Cyclotomic Field of order 8 and degree 4
              Defn: zeta8 |--> -zeta8^3
            sage: QuadraticField(5, 'a').complex_conjugation()
            Identity endomorphism of Number Field in a with defining
             polynomial x^2 - 5 with a = 2.236067977499790?
            sage: x = polygen(QQ, 'x')
            sage: F = NumberField(x^4 + x^3 - 3*x^2 - x + 1, 'a')
            sage: F.is_totally_real()
            True
            sage: F.complex_conjugation()
            Identity endomorphism of Number Field in a with defining
             polynomial x^4 + x^3 - 3*x^2 - x + 1
            sage: F.<b> = NumberField(x^2 - 2)
            sage: F.extension(x^2 + 1, 'a').complex_conjugation()
            Relative number field endomorphism of Number Field in a
             with defining polynomial x^2 + 1 over its base field
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
        # Return cached answer if available
        try:
            return self.__complex_conjugation
        except (AttributeError):
            pass

        # Then, deal with simple cases
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

        # In the remaining case, self.is_CM() should have cached __max_tot_real_sub
        try:
            F, phi = self.__max_tot_real_sub
        except (AttributeError):
            F, phi = self.maximal_totally_real_subfield()
        if self.is_absolute():
            K_rel = self.relativize(phi, self.variable_name() * 2)
            to_abs, from_abs = K_rel.structure()
            self.__complex_conjugation = K_rel.automorphisms()[1].pre_compose(
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
        Return the maximal totally real subfield of ``self`` together with an embedding of it into ``self``.

        EXAMPLES::

            sage: F.<a> = QuadraticField(11)
            sage: F.maximal_totally_real_subfield()
            [Number Field in a with defining polynomial x^2 - 11 with a = 3.316624790355400?,
             Identity endomorphism of
              Number Field in a with defining polynomial x^2 - 11 with a = 3.316624790355400?]
            sage: F.<a> = QuadraticField(-15)
            sage: F.maximal_totally_real_subfield()
            [Rational Field, Natural morphism:
               From: Rational Field
               To:   Number Field in a with defining polynomial x^2 + 15
                     with a = 3.872983346207417?*I]
            sage: F.<a> = CyclotomicField(29)
            sage: F.maximal_totally_real_subfield()
            (Number Field in a0 with defining polynomial x^14 + x^13 - 13*x^12 - 12*x^11
               + 66*x^10 + 55*x^9 - 165*x^8 - 120*x^7 + 210*x^6 + 126*x^5 - 126*x^4
               - 56*x^3 + 28*x^2 + 7*x - 1 with a0 = 1.953241111420174?,
             Ring morphism:
               From: Number Field in a0 with defining polynomial x^14 + x^13 - 13*x^12 - 12*x^11
                     + 66*x^10 + 55*x^9 - 165*x^8 - 120*x^7 + 210*x^6 + 126*x^5 - 126*x^4
                     - 56*x^3 + 28*x^2 + 7*x - 1 with a0 = 1.953241111420174?
               To:   Cyclotomic Field of order 29 and degree 28
               Defn: a0 |--> -a^27 - a^26 - a^25 - a^24 - a^23 - a^22 - a^21 - a^20 - a^19
                             - a^18 - a^17 - a^16 - a^15 - a^14 - a^13 - a^12 - a^11 - a^10
                             - a^9 - a^8 - a^7 - a^6 - a^5 - a^4 - a^3 - a^2 - 1)
            sage: x = polygen(QQ, 'x')
            sage: F.<a> = NumberField(x^3 - 2)
            sage: F.maximal_totally_real_subfield()
            [Rational Field,
             Coercion map:
               From: Rational Field
               To:   Number Field in a with defining polynomial x^3 - 2]
            sage: F.<a> = NumberField(x^4 - x^3 - x^2 + x + 1)
            sage: F.maximal_totally_real_subfield()
            [Rational Field,
             Coercion map:
               From: Rational Field
               To:   Number Field in a with defining polynomial x^4 - x^3 - x^2 + x + 1]
            sage: F.<a> = NumberField(x^4 - x^3 + 2*x^2 + x + 1)
            sage: F.maximal_totally_real_subfield()
            [Number Field in a1 with defining polynomial x^2 - x - 1,
             Ring morphism:
              From: Number Field in a1 with defining polynomial x^2 - x - 1
              To:   Number Field in a with defining polynomial x^4 - x^3 + 2*x^2 + x + 1
              Defn: a1 |--> -1/2*a^3 - 1/2]
            sage: F.<a> = NumberField(x^4 - 4*x^2 - x + 1)
            sage: F.maximal_totally_real_subfield()
            [Number Field in a with defining polynomial x^4 - 4*x^2 - x + 1,
             Identity endomorphism of
              Number Field in a with defining polynomial x^4 - 4*x^2 - x + 1]

        An example of a relative extension where the base field is not the maximal totally real subfield.

        ::

            sage: E_0.<a> = NumberField(x^2 - 4*x + 16)
            sage: y = polygen(E_0)
            sage: E.<z> = E_0.extension(y^2 - E_0.gen() / 2)
            sage: E.maximal_totally_real_subfield()
            [Number Field in z1 with defining polynomial x^2 - 2*x - 5,
             Composite map:
               From: Number Field in z1 with defining polynomial x^2 - 2*x - 5
               To:   Number Field in z with defining polynomial x^2 - 1/2*a over its base field
               Defn:   Ring morphism:
                       From: Number Field in z1 with defining polynomial x^2 - 2*x - 5
                       To:   Number Field in z with defining
                             polynomial x^4 - 2*x^3 + x^2 + 6*x + 3
                       Defn: z1 |--> -1/3*z^3 + 1/3*z^2 + z - 1
                     then
                       Isomorphism map:
                       From: Number Field in z with defining
                             polynomial x^4 - 2*x^3 + x^2 + 6*x + 3
                       To:   Number Field in z with defining
                             polynomial x^2 - 1/2*a over its base field]
        """

        try:
            return self.__max_tot_real_sub
        except (AttributeError):
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
        complex field with precision ``prec``.

        This always embeds into an MPFR based complex field.  If you
        want embeddings into the 53-bit double precision, which is
        faster, use ``self.embeddings(CDF)``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^5 + x + 17)
            sage: v = k.complex_embeddings()
            sage: ls = [phi(k.0^2) for phi in v]; ls  # random order
            [2.97572074038...,
             -2.40889943716 + 1.90254105304*I,
             -2.40889943716 - 1.90254105304*I,
             0.921039066973 + 3.07553311885*I,
             0.921039066973 - 3.07553311885*I]
            sage: K.<a> = NumberField(x^3 + 2)
            sage: ls = K.complex_embeddings(); ls  # random order
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
        CC = sage.rings.complex_mpfr.ComplexField(prec)
        return self.embeddings(CC)

    def real_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this number field into the approximate
        real field with precision ``prec``.

        If ``prec`` is 53 (the default), then the real double field is
        used; otherwise the arbitrary precision (but slow) real field
        is used.  If you want embeddings into the 53-bit double
        precision, which is faster, use ``self.embeddings(RDF)``.

        .. NOTE::

            This function uses finite precision real numbers.
            In functions that should output proven results, one
            could use ``self.embeddings(AA)`` instead.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.real_embeddings()
            [Ring morphism:
               From: Number Field in a with defining polynomial x^3 + 2
               To:   Real Field with 53 bits of precision
               Defn: a |--> -1.25992104989487]
            sage: K.real_embeddings(16)
            [Ring morphism:
               From: Number Field in a with defining polynomial x^3 + 2
               To:   Real Field with 16 bits of precision
               Defn: a |--> -1.260]
            sage: K.real_embeddings(100)
            [Ring morphism:
               From: Number Field in a with defining polynomial x^3 + 2
               To:   Real Field with 100 bits of precision
               Defn: a |--> -1.2599210498948731647672106073]

        As this is a numerical function, the number of embeddings
        may be incorrect if the precision is too low::

            sage: K = NumberField(x^2 + 2*10^1000*x + 10^2000 + 1, 'a')
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
        Return the embedding of this field into the complex numbers which has
        been specified.

        Fields created with the :func:`QuadraticField` or
        :func:`CyclotomicField` constructors come with an implicit
        embedding. To get one of these fields without the embedding, use
        the generic :class:`NumberField` constructor.

        EXAMPLES::

            sage: QuadraticField(-1, 'I').specified_complex_embedding()
            Generic morphism:
              From: Number Field in I with defining polynomial x^2 + 1 with I = 1*I
              To:   Complex Lazy Field
              Defn: I -> 1*I

        ::

            sage: QuadraticField(3, 'a').specified_complex_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^2 - 3
                    with a = 1.732050807568878?
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

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 - 2, 'a').specified_complex_embedding() is None
            True
            sage: NumberField(x^3 - x + 5, 'a').specified_complex_embedding() is None
            True
            sage: NumberField(x^3 - x + 5, 'a', embedding=2).specified_complex_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - x + 5
                    with a = -1.904160859134921?
              To:   Real Lazy Field
              Defn: a -> -1.904160859134921?
            sage: NumberField(x^3 - x + 5, 'a', embedding=CDF.0).specified_complex_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - x + 5
                    with a = 0.952080429567461? + 1.311248044077123?*I
              To:   Complex Lazy Field
              Defn: a -> 0.952080429567461? + 1.311248044077123?*I

        This function only returns complex embeddings::

            sage: # needs sage.rings.padics
            sage: K.<a> = NumberField(x^2 - 2, embedding=Qp(7)(2).sqrt())
            sage: K.specified_complex_embedding() is None
            True
            sage: K.gen_embedding()
            3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8 + 6*7^9 + 6*7^10
             + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16 + 7^17 + 4*7^18 + 6*7^19 + O(7^20)
            sage: K.coerce_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial x^2 - 2
                    with a = 3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8
                              + 6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16
                              + 7^17 + 4*7^18 + 6*7^19 + O(7^20)
              To:   7-adic Field with capped relative precision 20
              Defn: a -> 3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8
                          + 6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16
                          + 7^17 + 4*7^18 + 6*7^19 + O(7^20)
        """
        embedding = self.coerce_embedding()
        if embedding is not None and embedding.codomain()._is_numerical():
            return embedding

    @cached_method
    def gen_embedding(self):
        """
        If an embedding has been specified, return the image of the
        generator under that embedding. Otherwise return ``None``.

        EXAMPLES::

            sage: QuadraticField(-7, 'a').gen_embedding()
            2.645751311064591?*I
            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 7, 'a').gen_embedding()  # None
        """
        embedding = self.coerce_embedding()
        if embedding is None:
            return None
        else:
            return embedding(self.gen())

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self`` (which is ``QQbar``).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: K.algebraic_closure()
            Algebraic Field
            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: K.algebraic_closure()
            Algebraic Field
            sage: K = CyclotomicField(23)
            sage: K.algebraic_closure()
            Algebraic Field
        """
        from sage.rings.qqbar import QQbar
        return QQbar

    @cached_method
    def conductor(self, check_abelian=True):
        r"""
        Compute the conductor of the abelian field `K`.
        If ``check_abelian`` is set to ``False`` and the field is not an
        abelian extension of `\QQ`, the output is not meaningful.

        INPUT:

        - ``check_abelian`` -- boolean (default: ``True``); check to see
          that this is an abelian extension of `\QQ`

        OUTPUT: integer which is the conductor of the field

        EXAMPLES::

            sage: # needs sage.groups
            sage: K = CyclotomicField(27)
            sage: k = K.subfields(9)[0][0]
            sage: k.conductor()
            27
            sage: x = polygen(QQ, 'x')
            sage: K.<t> = NumberField(x^3 + x^2 - 2*x - 1)
            sage: K.conductor()
            7
            sage: K.<t> = NumberField(x^3 + x^2 - 36*x - 4)
            sage: K.conductor()
            109
            sage: K = CyclotomicField(48)
            sage: k = K.subfields(16)[0][0]
            sage: k.conductor()
            48
            sage: NumberField(x,'a').conductor()
            1
            sage: NumberField(x^8 - 8*x^6 + 19*x^4 - 12*x^2 + 1, 'a').conductor()
            40
            sage: NumberField(x^8 + 7*x^4 + 1, 'a').conductor()
            40
            sage: NumberField(x^8 - 40*x^6 + 500*x^4 - 2000*x^2 + 50, 'a').conductor()
            160

        ALGORITHM:

            For odd primes, it is easy to compute from the ramification
            index because the `p`-Sylow subgroup is cyclic.  For `p=2`, there
            are two choices for a given ramification index.  They can be
            distinguished by the parity of the exponent in the discriminant
            of a 2-adic completion.
        """
        m = 1
        if check_abelian and not self.is_abelian():
            raise ValueError("The conductor is only defined for abelian fields")

        try:
            De = self.__disc
        except AttributeError:
            De = self.polynomial().discriminant()
            A = De.numerator().prime_factors()+De.denominator().prime_factors()
        else:
            A = De.prime_factors()

        for p in A:
            R = self.maximal_order(p, assume_maximal=None)
            e = R.fractional_ideal(p).prime_factors()[0].ramification_index()
            if e != 1:
                if p == 2:
                    m *= e*2
                    c = R.discriminant().valuation(2)
                    c /= self.polynomial().degree()/e
                    if is_odd(c):
                        m *= 2
                else:
                    m *= p**(e.valuation(p)+1)
        return m

    def dirichlet_group(self):
        r"""
        Given a abelian field `K`, compute and return the
        set of all Dirichlet characters corresponding to the
        characters of the Galois group of `K/\QQ`.

        The output is random if the field is not abelian.

        OUTPUT: list of Dirichlet characters

        EXAMPLES::

            sage: # needs sage.groups sage.modular
            sage: x = polygen(QQ, 'x')
            sage: K.<t> = NumberField(x^3 + x^2 - 36*x - 4)
            sage: K.conductor()
            109
            sage: K.dirichlet_group()  # optional - gap_package_polycyclic
            [Dirichlet character modulo 109 of conductor 1 mapping 6 |--> 1,
             Dirichlet character modulo 109 of conductor 109 mapping 6 |--> zeta3,
             Dirichlet character modulo 109 of conductor 109 mapping 6 |--> -zeta3 - 1]

            sage: # needs sage.modular
            sage: K = CyclotomicField(44)
            sage: L = K.subfields(5)[0][0]
            sage: X = L.dirichlet_group(); X  # optional - gap_package_polycyclic
            [Dirichlet character modulo 11 of conductor 1 mapping 2 |--> 1,
             Dirichlet character modulo 11 of conductor 11 mapping 2 |--> zeta5,
             Dirichlet character modulo 11 of conductor 11 mapping 2 |--> zeta5^2,
             Dirichlet character modulo 11 of conductor 11 mapping 2 |--> zeta5^3,
             Dirichlet character modulo 11 of conductor 11
               mapping 2 |--> -zeta5^3 - zeta5^2 - zeta5 - 1]
            sage: X[4]^2  # optional - gap_package_polycyclic
            Dirichlet character modulo 11 of conductor 11 mapping 2 |--> zeta5^3
            sage: X[4]^2 in X  # optional - gap_package_polycyclic
            True
        """
        # todo : turn this into an abelian group rather than a list.
        from sage.modular.dirichlet import DirichletGroup

        m = self.conductor()
        d = self.degree()
        A = _splitting_classes_gens_(self, m, d)
        # d could be improved to be the exponent of the Galois group
        # rather than the degree, but I do not see how to go about that yet.
        G = DirichletGroup(m, CyclotomicField(d))
        H = [G.one()]
        for chi in G:
            if len(H) == d:
                break
            if chi not in H:
                if all(chi(a) == 1 for a in A):
                    H.append(chi)
        return H

    def _repr_(self):
        """
        Return string representation of this number field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._repr_()
            'Number Field in a with defining polynomial x^13 - 2/3*x + 3'
            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3, embedding=-1)
            sage: k._repr_()
            'Number Field in a with defining polynomial x^13 - 2/3*x + 3 with a = -1.106745229567614?'
        """
        result = "Number Field in {} with defining polynomial {}".format(self.variable_name(), self.polynomial())
        gen = self.gen_embedding()
        if gen is not None:
            result += " with {} = {}".format(self.variable_name(), gen)
        return result

    def _latex_(self):
        r"""
        Return latex representation of this number field. This is viewed as
        a polynomial quotient ring over a field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._latex_()
            '\\Bold{Q}[a]/(a^{13} - \\frac{2}{3} a + 3)'
            sage: latex(k)
            \Bold{Q}[a]/(a^{13} - \frac{2}{3} a + 3)

        Numbered variables are often correctly typeset::

            sage: k.<theta25> = NumberField(x^25+x+1)
            sage: print(k._latex_())
            \Bold{Q}[\theta_{25}]/(\theta_{25}^{25} + \theta_{25} + 1)
        """
        latex_name = self.latex_variable_names()[0]
        return "%s[%s]/(%s)" % (latex(QQ), latex_name,
                                self.polynomial()._latex_(latex_name))

    def _ideal_class_(self, n=0):
        """
        Return the Python class used in defining the zero ideal of the ring
        of integers of this number field.

        This function is required by the general ring/ideal machinery. The
        value defined here is the default value for all number fields.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        one of the same name in class :class:`NumberField_relative`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 2, 'c')._fractional_ideal_class_()
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
        """
        return sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal

    def ideal(self, *gens, **kwds):
        """
        Return a fractional ideal of the field, except for the
        zero ideal, which is not a fractional ideal.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.ideal(2)
            Fractional ideal (2)
            sage: K.ideal(2 + i)
            Fractional ideal (i + 2)
            sage: K.ideal(0)
            Ideal (0) of Number Field in i with defining polynomial x^2 + 1

        TESTS:

        Check that :issue:`25934` is fixed::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^6 - x^5 - 5*x^4 + 4*x^3 + 6*x^2 - 3*x - 1)
            sage: K.ideal(1, 1)
            Fractional ideal (1)
        """
        try:
            return self.fractional_ideal(*gens, **kwds)
        except ValueError:
            return self.zero_ideal()

    def idealchinese(self, ideals, residues):
        r"""
        Return a solution of the Chinese Remainder Theorem problem
        for ideals in a number field.

        This is a wrapper around the pari function :pari:`idealchinese`.

        INPUT:

        - ``ideals`` -- list of ideals of the number field

        - ``residues`` -- list of elements of the number field

        OUTPUT:

        Return an element `b` of the number field such that
        `b \equiv x_i \bmod I_i` for all residues `x_i` and
        respective ideals `I_i`.

        .. SEEALSO::

            - :func:`crt`

        EXAMPLES:

        This is the example from the pari page on ``idealchinese``::

            sage: # needs sage.symbolic
            sage: K.<sqrt2> = NumberField(sqrt(2).minpoly())
            sage: ideals = [K.ideal(4), K.ideal(3)]
            sage: residues = [sqrt2, 1]
            sage: r = K.idealchinese(ideals, residues); r
            -3*sqrt2 + 4
            sage: all((r - a) in I for I, a in zip(ideals, residues))
            True

        The result may be non-integral if the results are non-integral::

            sage: # needs sage.symbolic
            sage: K.<sqrt2> = NumberField(sqrt(2).minpoly())
            sage: ideals = [K.ideal(4), K.ideal(21)]
            sage: residues = [1/sqrt2, 1]
            sage: r = K.idealchinese(ideals, residues); r
            -63/2*sqrt2 - 20
            sage: all(
            ....:     (r - a).valuation(P) >= k
            ....:     for I, a in zip(ideals, residues)
            ....:     for P, k in I.factor()
            ....: )
            True
        """
        factorizations = [I.factor() for I in ideals]
        y = [a for a, f in zip(residues, factorizations) for _ in f]
        x = pari.Mat([
            pari.Col([p.pari_prime(), k])
            for f in factorizations
            for p, k in f
        ]).mattranspose()
        r = self.pari_nf().idealchinese(x, y)
        return self(r)

    def fractional_ideal(self, *gens, **kwds):
        r"""
        Return the ideal in `\mathcal{O}_K` generated by ``gens``.

        This overrides the :class:`sage.rings.ring.Field` method to
        use the :class:`sage.rings.ring.Ring` one instead, since
        we are not concerned with ideals in a field but in its ring
        of integers.

        INPUT:

        - ``gens`` -- list of generators, or a number field ideal

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
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
            Fractional ideal (-a)
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
        r"""
        Return all integral ideals of bounded norm.

        INPUT:

        - ``bound`` -- positive integer

        OUTPUT: a dict of all integral ideals `I` such that Norm(`I`) `\leq` ``bound``,
        keyed by norm.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 23)
            sage: d = K.ideals_of_bdd_norm(10)
            sage: for n in d:
            ....:     print(n)
            ....:     for I in sorted(d[n]):
            ....:         print(I)
                1
                Fractional ideal (1)
                2
                Fractional ideal (2, 1/2*a - 1/2)
                Fractional ideal (2, 1/2*a + 1/2)
                3
                Fractional ideal (3, 1/2*a - 1/2)
                Fractional ideal (3, 1/2*a + 1/2)
                4
                Fractional ideal (2)
                Fractional ideal (4, 1/2*a + 3/2)
                Fractional ideal (4, 1/2*a + 5/2)
                5
                6
                Fractional ideal (-1/2*a + 1/2)
                Fractional ideal (1/2*a + 1/2)
                Fractional ideal (6, 1/2*a + 5/2)
                Fractional ideal (6, 1/2*a + 7/2)
                7
                8
                Fractional ideal (4, a - 1)
                Fractional ideal (4, a + 1)
                Fractional ideal (-1/2*a - 3/2)
                Fractional ideal (1/2*a - 3/2)
                9
                Fractional ideal (3)
                Fractional ideal (9, 1/2*a + 7/2)
                Fractional ideal (9, 1/2*a + 11/2)
                10
            sage: [[I.norm() for I in sorted(d[n])] for n in d]
                [[1], [2, 2], [3, 3], [4, 4, 4], [], [6, 6, 6, 6], [], [8, 8, 8, 8], [9, 9, 9], []]
        """
        hnf_ideals = self.pari_nf().ideallist(bound)
        d = {}
        for i in range(bound):
            d[i+1] = [self.ideal(hnf) for hnf in hnf_ideals[i]]
        return d

    def primes_above(self, x, degree=None):
        r"""
        Return prime ideals of ``self`` lying over `x`.

        INPUT:

        - ``x`` -- usually an element or ideal of ``self``. It should be such
          that ``self.ideal(x)`` is sensible. This excludes `x=0`.

        - ``degree`` -- (default: ``None``) ``None`` or an integer.
          If ``None``, find all primes above `x` of any degree. If an
          integer, find all primes above `x` such that the resulting
          residue field has exactly this degree.

        OUTPUT: list of prime ideals of ``self`` lying over `x`. If ``degree``
        is specified and no such ideal exists, returns the empty list.
        The output is sorted by residue degree first, then by
        underlying prime (or equivalently, by norm).

        If there is a tie, the exact ordering should be assumed to be random.
        See the remark in :meth:`NumberFieldIdeal._richcmp_`.

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

        ::

            sage: P2s = F.primes_above(2)
            sage: P2s  # random
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

        The ideal `(3)` is totally ramified in `F`, so there is no degree 2
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

        It doesn't make sense to factor the ideal `(0)`, so this raises an error::

            sage: F.prime_above(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'prime_factors'...
        """
        if degree is not None:
            degree = ZZ(degree)
        facs = sorted((id.residue_class_degree(), id.absolute_norm(), id)
                      for id in self.prime_factors(x))
        if degree is None:
            return [id for d, n, id in facs]
        else:
            return [id for d, n, id in facs if d == degree]

    def prime_above(self, x, degree=None):
        r"""
        Return a prime ideal of ``self`` lying over `x`.

        INPUT:

        - ``x`` -- usually an element or ideal of ``self``. It should be such
          that ``self.ideal(x)`` is sensible. This excludes `x=0`.

        - ``degree`` -- (default: ``None``) ``None`` or an integer.
          If one, find a prime above `x` of any degree. If an integer, find a
          prime above `x` such that the resulting residue field has exactly
          this degree.

        OUTPUT: a prime ideal of ``self`` lying over `x`. If ``degree`` is specified
        and no such ideal exists, raises a :exc:`ValueError`.

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

        The ideal `(3)` is totally ramified in `F`, so there is no degree 2
        prime above `3`::

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

        It doesn't make sense to factor the ideal `(0)`::

            sage: F.prime_above(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'prime_factors'...
        """
        ids = self.primes_above(x, degree)
        if not ids:
            raise ValueError("No prime of degree %s above %s" % (degree, self.ideal(x)))
        return ids[0]

    def primes_of_bounded_norm(self, B):
        r"""
        Return a sorted list of all prime ideals with norm at most `B`.

        INPUT:

        - ``B`` -- positive integer or real; upper bound on the norms of the
          primes generated

        OUTPUT:

        A list of all prime ideals of this number field of norm at
        most `B`, sorted by norm.  Primes of the same norm are sorted
        using the comparison function for ideals, which is based on
        the Hermite Normal Form.

        .. NOTE::

            See also :meth:`primes_of_bounded_norm_iter` for an
            iterator version of this, but note that the iterator sorts
            the primes in order of underlying rational prime, not by
            norm.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: K.primes_of_bounded_norm(10)
            [Fractional ideal (i - 1),
             Fractional ideal (2*i - 1),
             Fractional ideal (-2*i - 1),
             Fractional ideal (3)]
            sage: K.primes_of_bounded_norm(1)
            []
            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P = K.primes_of_bounded_norm(30)
            sage: P
            [Fractional ideal (a),
             Fractional ideal (a + 1),
             Fractional ideal (a^2 + 1),
             Fractional ideal (a^2 + a - 1),
             Fractional ideal (2*a + 1),
             Fractional ideal (2*a^2 + a + 1),
             Fractional ideal (a^2 - 2*a - 1),
             Fractional ideal (a + 3)]
            sage: [p.norm() for p in P]
            [2, 3, 5, 11, 17, 23, 25, 29]
        """
        try:
            B = ZZ(B)
        except (TypeError, AttributeError):
            try:
                B = ZZ(B.ceil())
            except (TypeError, AttributeError):
                raise TypeError("%s is not valid bound on prime ideals" % B)
        if B < 2:
            return []

        from sage.rings.fast_arith import prime_range
        if self is QQ:
            return prime_range(B + 1, algorithm='pari_isprime')
        else:
            P = (pp for p in prime_range(B + 1, algorithm='pari_isprime')
                 for pp in self.primes_above(p))
            P = [p for p in P if p.norm() <= B]
            P.sort(key=lambda P: (P.norm(), P))
            return P

    def primes_of_bounded_norm_iter(self, B):
        r"""
        Iterator yielding all prime ideals with norm at most `B`.

        INPUT:

        - ``B`` -- positive integer or real; upper bound on the norms of the
          primes generated

        OUTPUT:

        An iterator over all prime ideals of this number field of norm
        at most `B`.

        .. NOTE::

            The output is not sorted by norm, but by size of the
            underlying rational prime.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: it = K.primes_of_bounded_norm_iter(10)
            sage: list(it)
            [Fractional ideal (i - 1),
             Fractional ideal (3),
             Fractional ideal (2*i - 1),
             Fractional ideal (-2*i - 1)]
            sage: list(K.primes_of_bounded_norm_iter(1))
            []
        """
        try:
            B = ZZ(B)
        except (TypeError, AttributeError):
            try:
                B = ZZ(B.ceil())
            except (TypeError, AttributeError):
                raise TypeError("%s is not valid bound on prime ideals" % B)

        if B < 2:
            return

        from sage.rings.fast_arith import prime_range
        if self is QQ:
            for p in prime_range(B + 1, algorithm='pari_isprime'):
                yield p
        else:
            for p in prime_range(B + 1, algorithm='pari_isprime'):
                for pp in self.primes_above(p):
                    if pp.norm() <= B:
                        yield pp

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

        - ``num_integer_primes`` -- (default: 10000) an
          integer. We try to find primes of absolute norm no greater than the
          ``num_integer_primes``-th prime number. For example, if
          ``num_integer_primes`` is 2, the largest norm found will be 3, since
          the second prime is 3.

        - ``max_iterations`` -- (default: 100) an integer. We
          test ``max_iterations`` integers to find small primes before raising
          :class:`StopIteration`.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(10)
            sage: it = K.primes_of_degree_one_iter()
            sage: Ps = [ next(it) for i in range(3) ]
            sage: Ps # random
            [Fractional ideal (z^3 + z + 1),
             Fractional ideal (3*z^3 - z^2 + z - 1),
             Fractional ideal (2*z^3 - 3*z^2 + z - 2)]
            sage: [P.norm() for P in Ps] # random
            [11, 31, 41]
            sage: [P.residue_class_degree() for P in Ps]
            [1, 1, 1]
        """
        from sage.rings.number_field.small_primes_of_degree_one import Small_primes_of_degree_one_iter
        return Small_primes_of_degree_one_iter(self, num_integer_primes, max_iterations)

    def primes_of_degree_one_list(self, n, num_integer_primes=10000, max_iterations=100):
        r"""
        Return a list of `n` prime ideals of absolute degree one and small
        norm.

        .. warning::

           It is possible that there are no primes of `K` of
           absolute degree one of small prime norm, and it is possible
           that this algorithm will not find any primes of small norm.

           See module :mod:`sage.rings.number_field.small_primes_of_degree_one`
           for details.

        INPUT:

        - ``num_integer_primes`` -- integer (default: 10000). We try to find
          primes of absolute norm no greater than the ``num_integer_primes``-th
          prime number. For example, if ``num_integer_primes`` is 2, the
          largest norm found will be 3, since the second prime is 3.

        - ``max_iterations`` -- integer (default: 100). We test ``max_iterations``
          integers to find small primes before raising :class:`StopIteration`.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(10)
            sage: Ps = K.primes_of_degree_one_list(3)
            sage: Ps  # random output
            [Fractional ideal (-z^3 - z^2 + 1),
             Fractional ideal (2*z^3 - 2*z^2 + 2*z - 3),
             Fractional ideal (2*z^3 - 3*z^2 + z - 2)]
            sage: [P.norm() for P in Ps]
            [11, 31, 41]
            sage: [P.residue_class_degree() for P in Ps]
            [1, 1, 1]
        """
        it = self.primes_of_degree_one_iter()
        return [next(it) for i in range(n)]

    def completely_split_primes(self, B=200):
        r"""
        Return a list of rational primes which split completely in the number field `K`.

        INPUT:

        - ``B`` -- positive integer bound (default: 200)

        OUTPUT: list of all primes `p < B` which split completely in ``K``

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<xi> = NumberField(x^3 - 3*x + 1)
            sage: K.completely_split_primes(100)
            [17, 19, 37, 53, 71, 73, 89]
        """
        from sage.rings.fast_arith import prime_range
        from sage.rings.finite_rings.finite_field_constructor import GF
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        split_primes = []
        for p in prime_range(B):
            Fp = GF(p)
            FpT = PolynomialRing(Fp, 'T')
            g = FpT(self.defining_polynomial())
            if len(factor(g)) == self.degree():
                split_primes.append(p)
        return split_primes

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        Return whether or not there is a homomorphism defined by the given
        images of generators.

        To do this we just check that the elements of the image of the
        given generator (``im_gens`` always has length 1) satisfies the
        relation of the defining poly of this field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        if len(im_gens) != 1:
            return False
        if base_map is None and not codomain.has_coerce_map_from(QQ):
            # We need that elements of the base ring of the polynomial
            # ring map canonically into codomain.
            return False
        f = self.defining_polynomial()
        try:
            # The current implementation won't productively use base_map
            # since the coefficients of f are in QQ, if there is a hom
            # from QQ to codomain it's probably unique and just the coercion
            if base_map is not None:
                f = f.map_coefficients(base_map)
            return codomain(f(im_gens[0])) == 0
        except (TypeError, ValueError):
            return False

    @cached_method
    def _pari_absolute_structure(self):
        r"""
        Return data relating the Sage and PARI absolute polynomials.

        OUTPUT:

        Let `L` be this number field, and let `f` be the defining
        polynomial of `K` over `\QQ`.  This method returns a triple
        ``(g, alpha, beta)``, where

        - ``g`` is the defining relative polynomial of the PARI ``nf``
          structure (see :meth:`pari_nf`);

        - ``alpha`` is the image of `x \bmod f` under some isomorphism
          `\phi\colon K[x]/(f) \to K[x]/(g)`

        - ``beta`` is the image of `x \bmod g` under the inverse
          isomorphism `\phi^{-1}\colon K[x]/(g) \to K[x]/(f)`

        EXAMPLES:

        If `f` is monic and integral, the result satisfies ``g = f``
        and ``alpha = beta = x``::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 - 2)
            sage: K._pari_absolute_structure()
            (y^2 - 2, Mod(y, y^2 - 2), Mod(y, y^2 - 2))

        An example where `f` neither monic nor integral::

            sage: K.<a> = NumberField(2*x^2 + 1/3)
            sage: K._pari_absolute_structure()
            (y^2 + 6, Mod(1/6*y, y^2 + 6), Mod(6*y, y^2 + 1/6))
        """
        f = self.absolute_polynomial()._pari_with_name('y')
        f = f * f.content().denominator()
        if f.pollead() == 1:
            g = f
            alpha = beta = g.variable().Mod(g)
        else:
            g, alpha = f.polredbest(flag=1)
            beta = alpha.modreverse()
        return g, alpha, beta

    def pari_polynomial(self, name='x'):
        """
        Return the PARI polynomial corresponding to this number field.

        INPUT:

        - ``name`` -- variable name (default: ``'x'``)

        OUTPUT:

        A monic polynomial with integral coefficients (PARI ``t_POL``)
        defining the PARI number field corresponding to ``self``.

        .. WARNING::

            This is *not* the same as simply converting the defining
            polynomial to PARI.

        EXAMPLES::

            sage: y = polygen(QQ)
            sage: k.<a> = NumberField(y^2 - 3/2*y + 5/3)
            sage: k.pari_polynomial()
            x^2 - x + 40
            sage: k.polynomial().__pari__()
            x^2 - 3/2*x + 5/3
            sage: k.pari_polynomial('a')
            a^2 - a + 40

        Some examples with relative number fields::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a, c> = NumberField([x^2 + 3, x^2 + 1])
            sage: k.pari_polynomial()
            x^4 + 8*x^2 + 4
            sage: k.pari_polynomial('a')
            a^4 + 8*a^2 + 4
            sage: k.absolute_polynomial()
            x^4 + 8*x^2 + 4
            sage: k.relative_polynomial()
            x^2 + 3

            sage: k.<a, c> = NumberField([x^2 + 1/3, x^2 + 1/4])
            sage: k.pari_polynomial()
            x^4 - x^2 + 1
            sage: k.absolute_polynomial()
            x^4 - x^2 + 1

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
        return self._pari_absolute_structure()[0].change_variable_name(name)

    def pari_nf(self, important=True):
        """
        Return the PARI number field corresponding to this field.

        INPUT:

        - ``important`` -- boolean (default: ``True``); if ``False``,
          raise a :exc:`RuntimeError` if we need to do a difficult
          discriminant factorization.  This is useful when an integral
          basis is not strictly required, such as for factoring
          polynomials over this number field.

        OUTPUT:

        The PARI number field obtained by calling the PARI function
        :pari:`nfinit` with ``self.pari_polynomial('y')`` as argument.

        .. NOTE::

            This method has the same effect as ``pari(self)``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
            [y^4 - 324*y + 2160, [0, 2], 48918708, 216, ..., [36, 36*y, y^3 + 6*y^2 - 252, -6*y^2], [1, 0, 0, 252; 0, 1, 0, 0; 0, 0, 0, 36; 0, 0, -6, 36], [1, 0, 0, 0, 0, 0, -18, -42, 0, -18, -46, 60, 0, -42, 60, -60; 0, 1, 0, 0, 1, 0, 2, 0, 0, 2, -11, 1, 0, 0, 1, 9; 0, 0, 1, 0, 0, 0, 6, -6, 1, 6, -5, 0, 0, -6, 0, 0; 0, 0, 0, 1, 0, -6, 6, -6, 0, 6, 1, 2, 1, -6, 2, 0]]
            sage: pari(k)
            [y^4 - 324*y + 2160, [0, 2], 48918708, 216, ...]
            sage: gp(k)
            [y^4 - 324*y + 2160, [0, 2], 48918708, 216, ...]

        With ``important=False``, we simply bail out if we cannot
        easily factor the discriminant::

            sage: p = next_prime(10^40); q = next_prime(10^41)
            sage: K.<a> = NumberField(x^2 - p*q)
            sage: K.pari_nf(important=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Unable to factor discriminant with trial division

        Next, we illustrate the ``maximize_at_primes`` and ``assume_disc_small``
        parameters of the :class:`NumberField` constructor. The following would take
        a very long time without the ``maximize_at_primes`` option::

            sage: K.<a> = NumberField(x^2 - p*q, maximize_at_primes=[p])
            sage: K.pari_nf()
            [y^2 - 100000000000000000000...]

        Since the discriminant is square-free, this also works::

            sage: K.<a> = NumberField(x^2 - p*q, assume_disc_small=True)
            sage: K.pari_nf()
            [y^2 - 100000000000000000000...]
        """
        try:
            return self._pari_nf
        except AttributeError:
            f = self.pari_polynomial("y")
            if f.poldegree() > 1:
                f = pari([f, self._pari_integral_basis(important=important)])
            self._pari_nf = f.nfinit()
            return self._pari_nf

    def pari_zk(self):
        r"""
        Integral basis of the PARI number field corresponding to this field.

        This is the same as ``pari_nf().getattr('zk')``, but much faster.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 - 17)
            sage: k.pari_zk()
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
            sage: k.pari_nf().getattr('zk')
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
        """
        return self.pari_nf().nf_get_zk()

    def __pari__(self):
        """
        Return the PARI number field corresponding to this field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k = NumberField(x^2 + x + 1, 'a')
            sage: k.__pari__()
            [y^2 + y + 1, [0, 1], -3, 1, ... [1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
            sage: pari(k)
            [y^2 + y + 1, [0, 1], -3, 1, ...[1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
        """
        return self.pari_nf()

    def _pari_init_(self):
        """
        Return the PARI number field corresponding to this field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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

        - ``proof`` -- if ``False``, assume GRH; if ``True``, run PARI's
          :pari:`bnfcertify` to make sure that the results are correct

        - ``units`` -- (default: ``True) if ``True``, insist on having
          fundamental units; if ``False``, the units may or may not be computed

        OUTPUT: the PARI ``bnf`` structure of this number field

        .. warning::

           Even with ``proof=True``, I wouldn't trust this to mean
           that everything computed involving this number field is
           actually correct.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: len(k.pari_bnf())
            10
            sage: k.pari_bnf()[:4]
            [[;], matrix(0,3), [;], ...]
            sage: len(k.pari_nf())
            9
            sage: k.<a> = NumberField(x^7 + 7); k
            Number Field in a with defining polynomial x^7 + 7
            sage: dummy = k.pari_bnf(proof=True)
        """
        proof = get_flag(proof, "number_field")
        # First compute bnf
        try:
            bnf = self._pari_bnf
        except AttributeError:
            f = self.pari_polynomial("y")
            _saved_rand = pari.getrand()
            # make this deterministic, it affects printing of ideals
            pari.setrand(1)
            if units:
                self._pari_bnf = f.bnfinit(1)
            else:
                self._pari_bnf = f.bnfinit()
            pari.setrand(_saved_rand)
            bnf = self._pari_bnf
        # Certify if needed
        if proof and not getattr(self, "_pari_bnf_certified", False):
            if bnf.bnfcertify() != 1:
                raise ValueError("The result is not correct according to bnfcertify")
            self._pari_bnf_certified = True
        return bnf

    def pari_rnfnorm_data(self, L, proof=True):
        """
        Return the PARI :pari:`rnfisnorminit` data corresponding to the
        extension `L` / ``self``.

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
            raise ValueError("L must be an extension of self")

        Kbnf = self.pari_bnf(proof=proof)
        return Kbnf.rnfisnorminit(L.pari_relative_polynomial())

    def _gap_init_(self):
        """
        Create a GAP object representing ``self`` and return its name.

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: z = QQ['z'].0
            sage: K.<zeta> = NumberField(z^2 - 2)
            sage: K._gap_init_()  # the following variable name $sage1 represents the F.base_ring() in gap and is somehow random
            'CallFuncList(function() local z,E; z:=Indeterminate($sage1,"z"); E:=AlgebraicExtension($sage1,z^2 - 2,"zeta"); return E; end,[])'
            sage: k = gap(K); k
            <algebraic extension over the Rationals of degree 2>
            sage: k.GeneratorsOfDivisionRing()
            [ zeta ]

        The following tests that it is possible to use a defining
        polynomial in the variable ``E``, even though by default
        ``E`` is used as a local variable in the above GAP
        ``CallFuncList``::

            sage: # needs sage.libs.gap
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
            raise NotImplementedError("Currently, only simple algebraic extensions are implemented in gap")
        from sage.interfaces.gap import gap as G

        q = self.polynomial()
        if q.variable_name() != 'E':
            return 'CallFuncList(function() local %s,E; %s:=Indeterminate(%s,"%s"); E:=AlgebraicExtension(%s,%s,"%s"); return E; end,[])' % (q.variable_name(), q.variable_name(), G(self.base_ring()).name(), q.variable_name(), G(self.base_ring()).name(), repr(self.polynomial()), str(self.gen()))

        return 'CallFuncList(function() local %s,F; %s:=Indeterminate(%s,"%s"); F:=AlgebraicExtension(%s,%s,"%s"); return F; end,[])' % (q.variable_name(), q.variable_name(), G(self.base_ring()).name(), q.variable_name(), G(self.base_ring()).name(), repr(self.polynomial()), str(self.gen()))

    def characteristic(self):
        """
        Return the characteristic of this number field, which is of course
        0.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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

        - ``proof`` -- if ``True`` (default), then compute the class group
          provably correctly; call :func:`number_field_proof` to change this
          default globally

        - ``names`` -- names of the generators of this class group

        OUTPUT: the class group of this number field

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 23)
            sage: G = K.class_group(); G
            Class group of order 3 with structure C3 of
             Number Field in a with defining polynomial x^2 + 23
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
            Class group of order 76 with structure C38 x C2 of
             Number Field in a with defining polynomial x^2 + 20072
            sage: G.0 # random
            Fractional ideal class (41, a + 10)
            sage: G.0^38
            Trivial principal fractional ideal class
            sage: G.1 # random
            Fractional ideal class (2, -1/2*a)
            sage: G.1^2
            Trivial principal fractional ideal class

        Class groups of Hecke polynomials tend to be very small::

            sage: # needs sage.modular
            sage: f = ModularForms(97, 2).T(2).charpoly()
            sage: f.factor()
            (x - 3) * (x^3 + 4*x^2 + 3*x - 1) * (x^4 - 3*x^3 - x^2 + 6*x - 1)
            sage: [NumberField(g,'a').class_group().order() for g,_ in f.factor()]
            [1, 1, 1]

        .. NOTE::

            Unlike in PARI/GP, class group computations *in Sage* do *not* by
            default assume the Generalized Riemann Hypothesis. To do class
            groups computations not provably correctly you must often pass the
            flag ``proof=False`` to functions or call the function
            ``proof.number_field(False)``. It can easily take 1000s of times
            longer to do computations with ``proof=True`` (the default).
        """
        proof = proof_flag(proof)
        try:
            return self.__class_group[proof, names]
        except KeyError:
            pass
        except AttributeError:
            self.__class_group = {}
        k = self.pari_bnf(proof)
        cycle_structure = tuple(ZZ(c) for c in k.bnf_get_cyc())

        # Gens is a list of ideals (the generators)
        gens = tuple(self.ideal(hnf) for hnf in k.bnf_get_gen())

        G = ClassGroup(cycle_structure, names, self, gens, proof=proof)
        self.__class_group[proof, names] = G
        return G

    def class_number(self, proof=None):
        """
        Return the class number of this number field, as an integer.

        INPUT:

        - ``proof`` -- boolean (default: ``True``, unless you called
          ``number_field_proof``)

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        Return the S-class group of this number field over its base field.

        INPUT:

        - ``S`` -- set of primes of the base field

        - ``proof`` -- if ``False``, assume the GRH in computing the class
          group. Default is ``True``. Call ``number_field_proof`` to change
          this default globally.

        - ``names`` -- names of the generators of this class group

        OUTPUT: the S-class group of this number field

        EXAMPLES:

        A well known example::

            sage: K.<a> = QuadraticField(-5)
            sage: K.S_class_group([])
            S-class group of order 2 with structure C2 of Number Field in a
             with defining polynomial x^2 + 5 with a = 2.236067977499790?*I

        When we include the prime `(2, a+1)`, the S-class group becomes
        trivial::

            sage: K.S_class_group([K.ideal(2, a + 1)])
            S-class group of order 1 of Number Field in a
             with defining polynomial x^2 + 5 with a = 2.236067977499790?*I

        TESTS::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2, a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S); CS
            S-class group of order 2 with structure C2 of Number Field in a
             with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
            sage: T = tuple()
            sage: CT = K.S_class_group(T); CT
            S-class group of order 4 with structure C4 of Number Field in a
             with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
            sage: K.class_group()
            Class group of order 4 with structure C4 of Number Field in a
             with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
        """
        proof = proof_flag(proof)
        if all(P.is_principal() for P in S):
            C = self.class_group(proof=proof)
            Slist = list(zip([g.ideal() for g in C.gens()], C.invariants()))
        else:
            Slist = self._S_class_group_and_units(tuple(S), proof=proof)[1]
        return SClassGroup(tuple(s[1] for s in Slist), names, self,
                           tuple(s[0] for s in Slist), tuple(S))

    def S_units(self, S, proof=True):
        """
        Return a list of generators of the S-units.

        INPUT:

        - ``S`` -- set of primes of the base field

        - ``proof`` -- if ``False``, assume the GRH in computing the class group

        OUTPUT: list of generators of the unit group

        .. NOTE::

            For more functionality see the function :func:`S_unit_group`.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: K.unit_group()
            Unit group with structure C6 of Number Field in a
             with defining polynomial x^2 + 3 with a = 1.732050807568878?*I
            sage: K.S_units([])  # random
            [1/2*a + 1/2]
            sage: K.S_units([])[0].multiplicative_order()
            6

        An example in a relative extension (see :issue:`8722`)::

            sage: x = polygen(QQ, 'x')
            sage: L.<a,b> = NumberField([x^2 + 1, x^2 - 5])
            sage: p = L.ideal((-1/2*b - 1/2)*a + 1/2*b - 1/2)
            sage: W = L.S_units([p]); [x.norm() for x in W]
            [9, 1, 1]

        Our generators should have the correct parent (:issue:`9367`)::

            sage: _.<x> = QQ[]
            sage: L.<alpha> = NumberField(x^3 + x + 1)
            sage: p = L.S_units([ L.ideal(7) ])
            sage: p[0].parent()
            Number Field in alpha with defining polynomial x^3 + x + 1

        TESTS:

        This checks that the multiple entries issue at :issue:`9341` is fixed::

            sage: _.<t> = QQ[]
            sage: K.<T> = NumberField(t - 1)
            sage: I = K.ideal(2)
            sage: K.S_units([I])
            [2, -1]
            sage: J = K.ideal(-2)
            sage: K.S_units([I, J, I])
            [2, -1]
        """
        return self._S_class_group_and_units(tuple(S), proof=proof)[0]

    @cached_method
    def _S_class_group_and_units(self, S, proof=True):
        """
        Compute S class group and units.

        INPUT:

        - ``S`` -- tuple of prime ideals of ``self``

        - ``proof`` -- if ``False``, assume the GRH in computing the class group

        OUTPUT:

        - ``units``: a list of generators of the unit group

        - ``clgp_gens`` -- a list of generators of the `S`-class group.
          Each generator is represented as a pair ``(gen, order)``,
          where ``gen`` is a fractional ideal of ``self`` and ``order`` is
          its order in the `S`-class group.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 5)
            sage: K._S_class_group_and_units(())
            ([-1], [(Fractional ideal (2, a + 1), 2)])

            sage: K.<a> = NumberField(polygen(QQ))
            sage: K._S_class_group_and_units( (K.ideal(5),) )
            ([5, -1], [])

        TESTS:

        Note for the following test that the representation of
        the units depends on the PARI version::

            sage: K.<a> = NumberField(x^3 - 381 * x + 127)
            sage: units, clpg_gens = K._S_class_group_and_units(tuple(K.primes_above(13)))
            sage: clpg_gens
            [(Fractional ideal (11, a - 2), 2), (Fractional ideal (19, a + 7), 2)]
            sage: units[:5]
            [2/13*a^2 + 1/13*a - 677/13,
             1/13*a^2 + 7/13*a - 332/13,
             -1/13*a^2 + 6/13*a + 345/13,
             -1,
             1/13*a^2 - 19/13*a - 7/13]
            sage: units[5] in (1/13*a^2 - 19/13*a - 7/13, 1/13*a^2 + 20/13*a - 7/13)
            True
            sage: len(units) == 6
            True

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(2*x^2 - 1/3)
            sage: K._S_class_group_and_units(tuple(K.primes_above(2) + K.primes_above(3)))
            ([6*a + 2, -6*a + 3, -1, -12*a - 5], [])  # 64-bit
            ([6*a + 2, -6*a - 3, -1, -12*a - 5], [])  # 32-bit
        """
        K_pari = self.pari_bnf(proof=proof)
        S_pari = [p.pari_prime() for p in sorted(set(S))]
        result = K_pari.bnfsunit(S_pari)
        units = [self(x, check=False) for x in result[0]] + self.unit_group().gens_values()
        orders = result[4][1].sage()
        gens = [self.ideal(_) for _ in result[4][2]]
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
            [0 1]
            [1 0]

        TESTS:

        Verify that :issue:`29364` is fixed::

            sage: R.<x> = QQ[]
            sage: L.<t> = NumberField(x^2 - 6058)
            sage: S = L.primes_above(2)
            sage: M = L._S_class_group_quotient_matrix(tuple(S))
            sage: M.dimensions()
            (2, 1)
            sage: CG = L.class_group()
            sage: SCG = L.S_class_group(S)
            sage: SCG(CG.0^M[0,0] * CG.1^M[1,0]) == SCG.0
            True
        """
        from sage.matrix.constructor import matrix
        S_clgp_gens = self._S_class_group_and_units(S)[1]
        a = len(S_clgp_gens)
        c = self.class_group().ngens()
        M = [u[0].ideal_class_log() for u in S_clgp_gens]
        M += [x.ideal_class_log() for x in S]
        M += list(matrix.diagonal(self.class_group().gens_orders()))
        M = matrix(ZZ, M)
        A, Q = M.hermite_form(transformation=True)
        assert A[:c] == 1 and A[c:] == 0
        return Q[:c, :a]

    def selmer_generators(self, S, m, proof=True, orders=False):
        r"""
        Compute generators of the group `K(S,m)`.

        INPUT:

        - ``S`` -- set of primes of ``self``

        - ``m`` -- positive integer

        - ``proof`` -- if ``False``, assume the GRH in computing the class group

        - ``orders`` -- boolean (default: ``False``); if ``True``, output two
          lists, the generators and their orders

        OUTPUT:

        A list of generators of `K(S,m)`, and (optionally) their
        orders as elements of `K^\times/(K^\times)^m`.  This is the
        subgroup of `K^\times/(K^\times)^m` consisting of elements `a`
        such that the valuation of `a` is divisible by `m` at all
        primes not in `S`.  It fits in an exact sequence between the
        units modulo `m`-th powers and the `m`-torsion in the
        `S`-class group:

        .. MATH::

            1                                    \longrightarrow
            O_{K,S}^\times / (O_{K,S}^\times)^m  \longrightarrow
            K(S,m)                               \longrightarrow
            \operatorname{Cl}_{K,S}[m]           \longrightarrow
            0.

        The group `K(S,m)` contains the subgroup of those `a` such
        that `K(\sqrt[m]{a})/K` is unramified at all primes of `K`
        outside of `S`, but may contain it properly when not all
        primes dividing `m` are in `S`.

        .. SEEALSO::

            :meth:`NumberField_generic.selmer_space`, which gives
            additional output when `m=p` is prime: as well as generators,
            it gives an abstract vector space over `\GF{p}` isomorphic to
            `K(S,p)` and maps implementing the isomorphism between this
            space and `K(S,p)` as a subgroup of `K^*/(K^*)^p`.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-5)
            sage: K.selmer_generators((), 2)
            [-1, 2]

        The previous example shows that the group generated by the
        output may be strictly larger than the group of
        elements giving extensions unramified outside `S`, since that
        has order just 2, generated by `-1`::

            sage: K.class_number()
            2
            sage: K.hilbert_class_field('b')
            Number Field in b with defining polynomial x^2 + 1 over its base field

        When `m` is prime all the orders are equal to `m`, but in general they are only divisors of `m`::

            sage: K.<a> = QuadraticField(-5)
            sage: P2 = K.ideal(2, -a + 1)
            sage: P3 = K.ideal(3, a + 1)
            sage: K.selmer_generators((), 2, orders=True)
            ([-1, 2], [2, 2])
            sage: K.selmer_generators((), 4, orders=True)
            ([-1, 4], [2, 2])
            sage: K.selmer_generators([P2], 2)
            [2, -1]
            sage: K.selmer_generators((P2,P3), 4)
            [2, -a - 1, -1]
            sage: K.selmer_generators((P2,P3), 4, orders=True)
            ([2, -a - 1, -1], [4, 4, 2])
            sage: K.selmer_generators([P2], 3)
            [2]
            sage: K.selmer_generators([P2, P3], 3)
            [2, -a - 1]
            sage: K.selmer_generators([P2, P3, K.ideal(a)], 3)  # random signs
            [2, a + 1, a]

        Example over `\QQ` (as a number field)::

            sage: K.<a> = NumberField(polygen(QQ))
            sage: K.selmer_generators([],5)
            []
            sage: K.selmer_generators([K.prime_above(p) for p in [2,3,5]],2)
            [2, 3, 5, -1]
            sage: K.selmer_generators([K.prime_above(p) for p in [2,3,5]],6, orders=True)
            ([2, 3, 5, -1], [6, 6, 6, 2])

        TESTS::

            sage: K.<a> = QuadraticField(-5)
            sage: P2 = K.ideal(2, -a + 1)
            sage: P3 = K.ideal(3, a + 1)
            sage: P5 = K.ideal(a)
            sage: S = K.selmer_generators([P2, P3, P5], 3)
            sage: S in ([2, a + 1, a], [2, a + 1, -a], [2, -a - 1, a], [2, -a - 1, -a]) or S
            True

        Verify that :issue:`14489` is fixed;
        the representation depends on the PARI version::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 381 * x + 127)
            sage: gens = K.selmer_generators(K.primes_above(13), 2)
            sage: len(gens) == 8
            True
            sage: gens[:5]
            [2/13*a^2 + 1/13*a - 677/13,
             1/13*a^2 + 7/13*a - 332/13,
             -1/13*a^2 + 6/13*a + 345/13,
             -1,
             1/13*a^2 - 19/13*a - 7/13]
            sage: gens[5] in (1/13*a^2 - 19/13*a - 7/13, 1/13*a^2 + 20/13*a - 7/13)
            True
            sage: gens[6] in (-1/13*a^2 + 45/13*a - 97/13, 1/13*a^2 - 45/13*a + 97/13)
            True
            sage: gens[7] in (2/13*a^2 + 40/13*a - 27/13, -2/13*a^2 - 40/13*a + 27/13)
            True

        Verify that :issue:`16708` is fixed::

            sage: K.<a> = QuadraticField(-5)
            sage: p = K.primes_above(2)[0]
            sage: S = K.selmer_generators((), 4)
            sage: all(4.divides(x.valuation(p)) for x in S)
            True
        """
        units, clgp_gens = self._S_class_group_and_units(tuple(S), proof=proof)
        gens = []
        ords = []
        for unit in units:
            order = unit.multiplicative_order()
            if order == Infinity:
                gens.append(unit)
                ords.append(m)
            else:
                m1 = order.gcd(m)
                if m1 != 1:
                    gens.append(unit)
                    ords.append(m1)
        card_S = len(S)
        if card_S != 0:
            from sage.matrix.constructor import Matrix
            H = self.class_group()
            gen_ords = [g.order() for g in H.gens()]
            pari_ords = pari(gen_ords).Col()
            Sords = [H(s).order() for s in S]
            MS = Matrix(ZZ, [H(s).exponents() for s in S]).transpose()
            pari_MS = pari(MS)
        for gen, order in clgp_gens:
            d = order.gcd(m)
            if d != 1:
                # The ideal I = gen^(order/d) has order d in Cl_S[m].
                # After multiplying by primes in S, the ideal
                # I^m = gen^(order*m/d) becomes principal.  We take
                # a generator of this ideal to get the corresponding
                # generator of the m-Selmer group.
                J = gen ** (order * m // d)
                if card_S != 0 and not J.is_principal():
                    B = H(J).exponents()
                    pari_B = (-pari(B)).Col()
                    exps = pari_MS.matsolvemod(pari_ords, pari_B).Vec().sage()
                    Spart = prod([S[i] ** (exps[i] % Sords[i]) for i in range(card_S)])
                    J *= Spart
                gens.append(self(J.gens_reduced()[0]))
                ords.append(d)
        if orders:
            return gens, ords
        else:
            return gens

    def selmer_group_iterator(self, S, m, proof=True):
        r"""
        Return an iterator through elements of the finite group `K(S,m)`.

        INPUT:

        - ``S`` -- set of primes of ``self``

        - ``m`` -- positive integer

        - ``proof`` -- if ``False``, assume the GRH in computing the class group

        OUTPUT:

        An iterator yielding the distinct elements of `K(S,m)`.  See
        the docstring for :meth:`NumberField_generic.selmer_generators` for
        more information.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-5)
            sage: list(K.selmer_group_iterator((), 2))
            [1, 2, -1, -2]
            sage: list(K.selmer_group_iterator((), 4))
            [1, 4, -1, -4]
            sage: list(K.selmer_group_iterator([K.ideal(2, -a + 1)], 2))
            [1, -1, 2, -2]
            sage: list(K.selmer_group_iterator([K.ideal(2, -a + 1), K.ideal(3, a + 1)], 2))
            [1, -1, -a - 1, a + 1, 2, -2, -2*a - 2, 2*a + 2]

        Examples over `\QQ` (as a number field)::

            sage: K.<a> = NumberField(polygen(QQ))
            sage: list(K.selmer_group_iterator([], 5))
            [1]
            sage: list(K.selmer_group_iterator([], 4))
            [1, -1]
            sage: list(K.selmer_group_iterator([K.prime_above(p) for p in [11,13]],2))
            [1, -1, 13, -13, 11, -11, 143, -143]
        """
        KSgens, ords = self.selmer_generators(S=S, m=m, proof=proof, orders=True)
        one = self.one()
        from sage.misc.mrange import cartesian_product_iterator
        for ev in cartesian_product_iterator([range(o) for o in ords]):
            yield prod([p ** e for p, e in zip(KSgens, ev)], one)

    def selmer_space(self, S, p, proof=None):
        r"""
        Compute the group `K(S,p)` as a vector space with maps to and from `K^*`.

        INPUT:

        - ``S`` -- set of primes ideals of ``self``

        - ``p`` -- a prime number

        - ``proof`` -- if ``False``, assume the GRH in computing the class group

        OUTPUT:

        (tuple) ``KSp``, ``KSp_gens``, ``from_KSp``, ``to_KSp`` where

        - ``KSp`` is an abstract vector space over `GF(p)` isomorphic to `K(S,p)`;

        - ``KSp_gens`` is a list of elements of `K^*` generating `K(S,p)`;

        - ``from_KSp`` is a function from ``KSp`` to `K^*`
          implementing the isomorphism from the abstract `K(S,p)` to
          `K(S,p)` as a subgroup of `K^*/(K^*)^p`;

        - ``to_KSP`` is a partial function from `K^*` to ``KSp``,
          defined on elements `a` whose image in `K^*/(K^*)^p` lies in
          `K(S,p)`, mapping them via the inverse isomorphism to the
          abstract vector space ``KSp``.

        The group `K(S,p)` is the finite subgroup of `K^*/(K^*)^p`
        consisting of elements whose valuation at all primes not in
        `S` is a multiple of `p`.  It contains the subgroup of those
        `a\in K^*` such that `K(\sqrt[p]{a})/K` is unramified at all
        primes of `K` outside of `S`, but may contain it properly when
        not all primes dividing `p` are in `S`.

        EXAMPLES:

        A real quadratic field with class number 2, where the fundamental
        unit is a generator, and the class group provides another
        generator when `p=2`::

            sage: K.<a> = QuadraticField(-5)
            sage: K.class_number()
            2
            sage: P2 = K.ideal(2, -a + 1)
            sage: P3 = K.ideal(3, a + 1)
            sage: P5 = K.ideal(a)
            sage: KS2, gens, fromKS2, toKS2 = K.selmer_space([P2, P3, P5], 2)
            sage: KS2
            Vector space of dimension 4 over Finite Field of size 2
            sage: gens
            [a + 1, a, 2, -1]

        Each generator must have even valuation at primes not in `S`::

            sage: [K.ideal(g).factor() for g in gens]
            [(Fractional ideal (2, a + 1)) * (Fractional ideal (3, a + 1)),
             Fractional ideal (-a), (Fractional ideal (2, a + 1))^2, 1]

            sage: toKS2(10)
            (0, 0, 1, 1)
            sage: fromKS2([0,0,1,1])
            -2
            sage: K(10/(-2)).is_square()
            True

            sage: KS3, gens, fromKS3, toKS3 = K.selmer_space([P2, P3, P5], 3)
            sage: KS3
            Vector space of dimension 3 over Finite Field of size 3
            sage: gens
            [1/2, 1/4*a + 1/4, a]

        An example to show that the group `K(S,2)` may be strictly
        larger than the group of elements giving extensions unramified
        outside `S`.  In this case, with `K` of class number `2` and
        `S` empty, there is only one quadratic extension of `K`
        unramified outside `S`, the Hilbert Class Field
        `K(\sqrt{-1})`::

            sage: K.<a> = QuadraticField(-5)
            sage: KS2, gens, fromKS2, toKS2 = K.selmer_space([], 2)
            sage: KS2
            Vector space of dimension 2 over Finite Field of size 2
            sage: gens
            [2, -1]
            sage: x = polygen(ZZ, 'x')
            sage: for v in KS2:
            ....:     if not v:
            ....:         continue
            ....:     a = fromKS2(v)
            ....:     print((a, K.extension(x^2 - a, 'roota').relative_discriminant().factor()))
            (2, (Fractional ideal (2, a + 1))^4)
            (-1, 1)
            (-2, (Fractional ideal (2, a + 1))^4)

            sage: K.hilbert_class_field('b')
            Number Field in b with defining polynomial x^2 + 1 over its base field
        """
        from sage.rings.number_field.selmer_group import pSelmerGroup
        return pSelmerGroup(self, S, p, proof)

    def composite_fields(self, other, names=None, both_maps=False, preserve_embedding=True):
        """
        Return the possible composite number fields formed from
        ``self`` and ``other``.

        INPUT:

        - ``other`` -- number field

        - ``names`` -- generator name for composite fields

        - ``both_maps`` -- boolean (default: ``False``)

        - ``preserve_embedding`` -- boolean (default: ``True``)

        OUTPUT: list of the composite fields, possibly with maps

        If ``both_maps`` is ``True``, the list consists of quadruples
        ``(F, self_into_F, other_into_F, k)`` such that
        ``self_into_F`` is an embedding of ``self`` in ``F``,
        ``other_into_F`` is an embedding of in ``F``, and ``k`` is one
        of the following:

        - an integer such that ``F.gen()`` equals
          ``other_into_F(other.gen()) + k*self_into_F(self.gen())``;

        - ``Infinity``, in which case ``F.gen()`` equals
          ``self_into_F(self.gen())``;

        - ``None`` (when ``other`` is a relative number field).

        If both ``self`` and ``other`` have embeddings into an ambient
        field, then each ``F`` will have an embedding with respect to
        which both ``self_into_F`` and ``other_into_F`` will be
        compatible with the ambient embeddings.

        If ``preserve_embedding`` is ``True`` and if ``self`` and
        ``other`` both have embeddings into the same ambient field, or
        into fields which are contained in a common field, only the
        compositum respecting both embeddings is returned.  In all
        other cases, all possible composite number fields are
        returned.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^4 - 2)
            sage: K.composite_fields(K)
            [Number Field in a with defining polynomial x^4 - 2,
             Number Field in a0 with defining polynomial x^8 + 28*x^4 + 2500]

        A particular compositum is selected, together with compatible maps
        into the compositum, if the fields are endowed with a real or
        complex embedding::

            sage: # needs sage.symbolic
            sage: K1 = NumberField(x^4 - 2, 'a', embedding=RR(2^(1/4)))
            sage: K2 = NumberField(x^4 - 2, 'a', embedding=RR(-2^(1/4)))
            sage: K1.composite_fields(K2)
            [Number Field in a with defining polynomial x^4 - 2 with a = 1.189207115002722?]
            sage: [F, f, g, k], = K1.composite_fields(K2, both_maps=True); F
            Number Field in a with defining polynomial x^4 - 2 with a = 1.189207115002722?
            sage: f(K1.0), g(K2.0)
            (a, -a)

        With ``preserve_embedding`` set to ``False``, the embeddings
        are ignored::

            sage: K1.composite_fields(K2, preserve_embedding=False)                     # needs sage.symbolic
            [Number Field in a with defining polynomial x^4 - 2 with a = 1.189207115002722?,
             Number Field in a0 with defining polynomial x^8 + 28*x^4 + 2500]

        Changing the embedding selects a different compositum::

            sage: K3 = NumberField(x^4 - 2, 'a', embedding=CC(2^(1/4)*I))               # needs sage.symbolic
            sage: [F, f, g, k], = K1.composite_fields(K3, both_maps=True); F            # needs sage.symbolic
            Number Field in a0 with defining polynomial x^8 + 28*x^4 + 2500
             with a0 = -2.378414230005443? + 1.189207115002722?*I
            sage: f(K1.0), g(K3.0)                                                      # needs sage.symbolic
            (1/240*a0^5 - 41/120*a0, 1/120*a0^5 + 19/60*a0)

        If no embeddings are specified, the maps into the compositum
        are chosen arbitrarily::

            sage: Q1.<a> = NumberField(x^4 + 10*x^2 + 1)
            sage: Q2.<b> = NumberField(x^4 + 16*x^2 + 4)
            sage: Q1.composite_fields(Q2, 'c')
            [Number Field in c with defining polynomial
              x^8 + 64*x^6 + 904*x^4 + 3840*x^2 + 3600]
            sage: F, Q1_into_F, Q2_into_F, k = Q1.composite_fields(Q2, 'c',
            ....:                                                  both_maps=True)[0]
            sage: Q1_into_F
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 + 10*x^2 + 1
              To:   Number Field in c with defining polynomial
                    x^8 + 64*x^6 + 904*x^4 + 3840*x^2 + 3600
              Defn: a |--> 19/14400*c^7 + 137/1800*c^5 + 2599/3600*c^3 + 8/15*c

        This is just one of four embeddings of ``Q1`` into ``F``::

            sage: Hom(Q1, F).order()
            4

        Note that even with ``preserve_embedding=True``, this method may fail
        to recognize that the two number fields have compatible embeddings, and
        hence return several composite number fields::

            sage: x = polygen(ZZ)
            sage: A.<a> = NumberField(x^3 - 7, embedding=CC(-0.95+1.65*I))
            sage: r = QQbar.polynomial_root(x^9 - 7, RIF(1.2, 1.3))
            sage: B.<a> = NumberField(x^9 - 7, embedding=r)
            sage: len(A.composite_fields(B, preserve_embedding=True))
            2

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
            [-0.6293842454258951, -0.7708351267200304]

        and we get the same embeddings via the compositum::

            sage: F, L1_into_F, L2_into_F, k = L1.composite_fields(L2, both_maps=True)[0]
            sage: [CDF(L1_into_F(L1.gen())), CDF(L2_into_F(L2.gen()))]
            [-0.6293842454258952, -0.7708351267200303]

        Let's check that if only one field has an embedding, the resulting
        fields do not have embeddings::

            sage: L1.composite_fields(K2)[0].coerce_embedding() is None
            True
            sage: L2.composite_fields(K1)[0].coerce_embedding() is None
            True

        We check that other can be a relative number field::

            sage: L.<a, b> = NumberField([x^3 - 5, x^2 + 3])
            sage: CyclotomicField(3, 'w').composite_fields(L, both_maps=True)
            [(Number Field in a with defining polynomial x^3 - 5 over its base field,
              Ring morphism:
                From: Cyclotomic Field of order 3 and degree 2
                To:   Number Field in a with defining polynomial x^3 - 5 over its base field
                Defn: w |--> -1/2*b - 1/2,
              Relative number field endomorphism of Number Field in a
               with defining polynomial x^3 - 5 over its base field
                Defn: a |--> a
                      b |--> b,
              None)]

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(x^2 + 1/2)
            sage: L.<b> = NumberField(3*x^2 - 1)
            sage: K.composite_fields(L)
            [Number Field in ab with defining polynomial 36*x^4 + 12*x^2 + 25]
            sage: C = K.composite_fields(L, both_maps=True); C
            [(Number Field in ab with defining polynomial 36*x^4 + 12*x^2 + 25,
              Ring morphism:
                From: Number Field in a with defining polynomial x^2 + 1/2
                To:   Number Field in ab with defining polynomial 36*x^4 + 12*x^2 + 25
                Defn: a |--> -3/5*ab^3 - 7/10*ab,
              Ring morphism:
                From: Number Field in b with defining polynomial 3*x^2 - 1
                To:   Number Field in ab with defining polynomial 36*x^4 + 12*x^2 + 25
                Defn: b |--> -3/5*ab^3 + 3/10*ab,
              -1)]
            sage: M, f, g, k = C[0]
            sage: M.gen() == g(b) + k*f(a)
            True

        This also fixes the bugs reported at :issue:`14164` and
        :issue:`18243`::

            sage: R.<x> = QQ[]
            sage: f = 6*x^5 + x^4 + x^2 + 5*x + 7
            sage: r = f.roots(QQbar, multiplicities=False)
            sage: F1 = NumberField(f.monic(), 'a', embedding=r[0])
            sage: F2 = NumberField(f.monic(), 'a', embedding=r[1])
            sage: (F, map1, map2, k) = F1.composite_fields(F2, both_maps=True)[0]
            sage: F.degree()
            20
            sage: F.gen() == map2(F2.gen()) + k*map1(F1.gen())
            True

            sage: f = x^8 - 3*x^7 + 61/3*x^6 - 9*x^5 + 298*x^4 + 458*x^3 + 1875*x^2 + 4293*x + 3099
            sage: F1 = NumberField(f, 'z', embedding=-1.18126721294295 + 3.02858651117832j)
            sage: F2 = NumberField(f, 'z', embedding=-1.18126721294295 - 3.02858651117832j)
            sage: (F, map1, map2, k) = F1.composite_fields(F2, both_maps=True)[0]
            sage: F.degree()
            32
            sage: F.gen() == map2(F2.gen()) + k*map1(F1.gen())
            True

        Check that the bugs reported at :issue:`24357` are fixed::

            sage: A.<a> = NumberField(x^9 - 7)
            sage: B.<b> = NumberField(x^3 - 7, embedding=a^3)
            sage: C.<c> = QuadraticField(-1)
            sage: B.composite_fields(C)
            [Number Field in bc with defining polynomial x^6 + 3*x^4 + 14*x^3 + 3*x^2 - 42*x + 50]

            sage: y = polygen(QQ, 'y')
            sage: A.<a> = NumberField(x^3 - 7, embedding=CC(-0.95+1.65*I))
            sage: B.<b> = NumberField(y^9 - 7, embedding=CC(-1.16+0.42*I))
            sage: A.composite_fields(B)
            [Number Field in b with defining polynomial y^9 - 7
              with b = -1.166502297945062? + 0.4245721146551276?*I]
        """
        if not isinstance(other, NumberField_generic):
            raise TypeError("other must be a number field.")

        sv = self.variable_name()
        ov = other.variable_name()
        if names is None:
            names = sv + (ov if ov != sv else "")
        name = normalize_names(1, names)[0]

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
            except sage.structure.coerce_exceptions.CoercionException:
                ambient_field = None
            if ambient_field is None:
                subfields_have_embeddings = False

        f = self.absolute_polynomial()
        g = other.absolute_polynomial().change_variable_name(f.variable_name())
        R = f.parent()
        f = f.__pari__()
        f /= f.content()
        g = g.__pari__()
        g /= g.content()

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
            C = [R(_) for _ in C]

            q = sum(1 for r in C if r.degree() != max(m, n))
            if q == 1 and name != sv and name != ov:
                names = [name]
            else:
                names = [name + str(i) for i in range(q)]

            i = 0
            rets = []
            for r in C:
                d = r.degree()
                if d == m:
                    rets.append(self)
                elif d == n:
                    rets.append(other)
                else:
                    rets.append(NumberField(r, names[i], check=False))
                    i += 1
            return rets

        # If flag = 1, polcompositum outputs a vector of 4-component vectors
        # [R, a, b, k], where R ranges through the list of all possible compositums
        # as above, and a (resp. b) expresses the root of P (resp. Q) as
        # an element of Q(X)/(R). Finally, k is a small integer such that
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
        # minimum value.
        if subfields_have_embeddings:
            poly_vals = []
            for r, _, _, k in C:
                r = R(r)
                k = ZZ(k)
                embedding = other.coerce_embedding()(b) + k*self.coerce_embedding()(a)
                poly_vals.append(r(embedding).abs())
            i = poly_vals.index(min(poly_vals))
            C = [C[i]]

        q = sum(1 for r, _, _, _ in C if r.poldegree() != max(m, n))
        if q == 1 and name != sv and name != ov:
            names = [name, '']
        else:
            names = [name + str(ii) for ii in range(q + 1)]

        if both_maps and not other.is_absolute():
            other_abs = other.absolute_field('z')
            from_other_abs, to_other_abs = other_abs.structure()

        embedding = None
        i = 0
        rets = []
        for r, a_in_F, b_in_F, k in C:
            r = R(r)
            d = r.degree()
            if d == m and not both_maps:
                rets.append(self)
            elif d == n and not both_maps:
                rets.append(other)
            else:
                k = ZZ(k)
                if subfields_have_embeddings:
                    embedding = other.coerce_embedding()(b) + k*self.coerce_embedding()(a)
                F = NumberField(r, names[i], check=False, embedding=embedding)
                i += 1
                if both_maps:
                    a_in_F = F(R(a_in_F.lift()))
                    b_in_F = F(R(b_in_F.lift()))
                    if other.is_absolute():
                        if d == m:
                            self_to_F = self.hom([self.gen()])
                            other_to_F = other.hom([(~self.hom([a_in_F]))(b_in_F)])
                            F = self
                            k = Infinity
                            i -= 1
                        elif d == n:
                            other_to_F = other.hom([other.gen()])
                            self_to_F = self.hom([(~other.hom([b_in_F]))(a_in_F)])
                            F = other
                            k = ZZ.zero()
                            i -= 1
                        else:
                            self_to_F = self.hom([a_in_F])
                            other_to_F = other.hom([b_in_F])
                    else:
                        other_abs_to_F = other_abs.hom([b_in_F])
                        other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(F), other_abs_to_F*to_other_abs)
                        if d == m:
                            self_to_F = self.hom([self.gen()])
                            other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(self), (~self.hom([a_in_F]))*other_abs_to_F*to_other_abs)
                            F = self
                            k = None
                            i -= 1
                        elif d == n:
                            other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(other), from_other_abs)
                            self_to_F = self.hom([from_other_abs((~other_abs_to_F)(a_in_F))])
                            F = other
                            k = None
                            i -= 1
                        else:
                            self_to_F = self.hom([a_in_F])
                            other_to_F = RelativeNumberFieldHomomorphism_from_abs(other.Hom(F), other_abs_to_F * to_other_abs)
                    rets.append((F, self_to_F, other_to_F, k))
                else:
                    rets.append(F)
        return rets

    def absolute_degree(self):
        r"""
        Return the degree of ``self`` over `\QQ`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').absolute_degree()
            3
            sage: NumberField(x + 1, 'a').absolute_degree()
            1
            sage: NumberField(x^997 + 17*x + 3, 'a', check=False).absolute_degree()
            997
        """
        return self.polynomial().degree()

    def degree(self):
        r"""
        Return the degree of this number field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        such that the trace of `xy` is an integer for
        all `y \in O_K`.

        The different is the integral ideal which is the inverse of
        the codifferent.

        See :wikipedia:`Different_ideal`

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: d = k.different()
            sage: d
            Fractional ideal (a)
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
        Return the discriminant of the ring of integers of the number
        field, or if ``v`` is specified, the determinant of the trace pairing
        on the elements of the list ``v``.

        INPUT:

        - ``v`` -- (optional) list of elements of this number field

        OUTPUT: integer if ``v`` is omitted, and Rational otherwise

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        if v is None:
            try:
                return self.__disc
            except AttributeError:
                self.__disc = ZZ(self.pari_polynomial().nfdisc())
                return self.__disc
        else:
            return QQ(self.trace_pairing(v).det())

    def disc(self, v=None):
        """
        Shortcut for :meth:`discriminant`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<b> = NumberField(x^2 - 123)
            sage: k.disc()
            492
        """
        return self.discriminant(v=v)

    def trace_dual_basis(self, b):
        r"""
        Compute the dual basis of a basis of ``self`` with respect to the trace pairing.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: b = [1, 2*a, 3*a^2]
            sage: T = K.trace_dual_basis(b); T
            [4/31*a^2 - 6/31*a + 13/31, -9/62*a^2 - 1/31*a - 3/31, 2/31*a^2 - 3/31*a + 4/93]
            sage: [(b[i]*T[j]).trace() for i in range(3) for j in range(3)]
            [1, 0, 0, 0, 1, 0, 0, 0, 1]
        """
        if not len(b) == self.degree():
            raise ValueError('Not a basis of the number field.')
        M = self.trace_pairing(b)
        if not M.is_invertible():
            raise ValueError('Not a basis of the number field.')
        return [sum([v[i]*b[i] for i in range(len(b))]) for v in M.inverse()]

    def elements_of_norm(self, n, proof=None) -> list:
        """
        Return a list of elements of norm `n`.

        INPUT:

        - ``n`` -- integer

        - ``proof`` -- boolean (default: ``True``, unless you called
          :meth:`proof.number_field` and set it otherwise)

        OUTPUT:

        A complete system of integral elements of norm `n`, modulo
        units of positive norm.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.elements_of_norm(3)
            []
            sage: K.elements_of_norm(50)
            [7*a - 1, 5*a - 5, -7*a - 1]

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`);
        the representation depends on the PARI version::

            sage: K.<a> = NumberField(7/9*x^3 + 7/3*x^2 - 56*x + 123)
            sage: [x] = K.elements_of_norm(7)
            sage: x in (7/225*a^2 - 7/75*a - 42/25, 28/225*a^2 + 77/75*a - 133/25)
            True
        """
        n = ZZ(n)
        proof = proof_flag(proof)
        B = self.pari_bnf(proof).bnfisintnorm(n)
        return [self(x, check=False) for x in B]

    def extension(self, poly, name=None, names=None, latex_name=None, latex_names=None, *args, **kwds):
        """
        Return the relative extension of this field by a given polynomial.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: R.<t> = K[]
            sage: L.<b> = K.extension(t^2 + a); L
            Number Field in b with defining polynomial t^2 + a over its base field

        We create another extension::

            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: y = polygen(QQ,'y')
            sage: m.<b> = k.extension(y^2 + 2); m
            Number Field in b with defining polynomial y^2 + 2 over its base field

        Note that `b` is a root of `y^2 + 2`::

            sage: b.minpoly()
            x^2 + 2
            sage: b.minpoly('z')
            z^2 + 2

        A relative extension of a relative extension::

            sage: k.<a> = NumberField([x^2 + 1, x^3 + x + 1])
            sage: R.<z> = k[]
            sage: L.<b> = NumberField(z^3 + 3 + a); L
            Number Field in b with defining polynomial z^3 + a0 + 3 over its base field

        Extension fields with given defining data are unique
        (:issue:`20791`)::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.extension(x^2 - 2, 'b') is K.extension(x^2 - 2, 'b')
            True
        """
        if not isinstance(poly, polynomial_element.Polynomial):
            try:
                poly = poly.polynomial(self)
            except (AttributeError, TypeError):
                raise TypeError("polynomial (=%s) must be a polynomial" % repr(poly))
        if poly.base_ring() is not self:
            poly = poly.change_ring(self)
        if names is not None:
            name = names
        if isinstance(name, (tuple, list)):
            name = name[0]
        if latex_names is not None:
            latex_name = latex_names
        if isinstance(latex_name, (tuple, list)):
            latex_name = latex_name[0]
        return NumberField(poly, name, latex_name=latex_name, *args, **kwds)

    def factor(self, n):
        r"""
        Ideal factorization of the principal ideal generated by `n`.

        EXAMPLES:

        Here we show how to factor Gaussian integers (up to units).
        First we form a number field defined by `x^2 + 1`::

            sage: x = polygen(QQ, 'x')
            sage: K.<I> = NumberField(x^2 + 1); K
            Number Field in I with defining polynomial x^2 + 1

        Here are the factors::

            sage: fi, fj = K.factor(17); fi,fj
            ((Fractional ideal (I + 4), 1), (Fractional ideal (I - 4), 1))

        Now we extract the reduced form of the generators::

            sage: zi = fi[0].gens_reduced()[0]; zi
            I + 4
            sage: zj = fj[0].gens_reduced()[0]; zj
            I - 4

        We recover the integer that was factored in `\ZZ[i]` (up to a unit)::

            sage: zi*zj
            -17

        One can also factor elements or ideals of the number field::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.factor(1/3)
            (Fractional ideal (3))^-1
            sage: K.factor(1+a)
            Fractional ideal (a - 1)
            sage: K.factor(1+a/5)
            (Fractional ideal (a - 1)) * (Fractional ideal (2*a - 1))^-1 * (Fractional ideal (-2*a - 1))^-1 * (Fractional ideal (3*a + 2))

        An example over a relative number field::

            sage: pari('setrand(2)')
            sage: L.<b> = K.extension(x^2 - 7)
            sage: f = L.factor(a + 1)
            sage: f                               # representation varies, not tested
            (Fractional ideal (1/2*a*b - a + 1/2)) * (Fractional ideal (-1/2*a*b - a + 1/2))
            sage: f.value() == a+1
            True

        It doesn't make sense to factor the ideal `(0)`, so this raises an error::

            sage: L.factor(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'factor'...

        AUTHORS:

        - Alex Clemesha (2006-05-20), Francis Clarke (2009-04-21): examples

        TESTS:

        We test the above doctest. The representation depends on the PARI version::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: L.<b> = K.extension(x^2 - 7)
            sage: f = L.factor(a + 1)
            sage: (fi, fj) = f[::]
            sage: (fi[1], fj[1])
            (1, 1)
            sage: fi[0] == L.fractional_ideal(-1/2*a*b - a + 1/2)
            True
            sage: fj[0] == L.fractional_ideal(1/2*a*b - a + 1/2)
            True
        """
        return self.ideal(n).factor()

    def prime_factors(self, x):
        """
        Return a list of the prime ideals of ``self`` which divide
        the ideal generated by `x`.

        OUTPUT: list of prime ideals (a new list is returned each time this
        function is called)

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<w> = NumberField(x^2 + 23)
            sage: K.prime_factors(w + 1)
            [Fractional ideal (2, 1/2*w - 1/2),
             Fractional ideal (2, 1/2*w + 1/2),
             Fractional ideal (3, 1/2*w + 1/2)]
        """
        return self.ideal(x).prime_factors()

    def decomposition_type(self, p):
        """
        Return how the given prime of the base field splits in this number field.

        INPUT:

        - ``p`` -- a prime element or ideal of the base field

        OUTPUT:

        A list of triples `(e, f, g)` where

        - `e` is the ramification index,

        - `f` is the residue class degree,

        - `g` is the number of primes above `p` with given `e` and `f`

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^20 + 3*x^18 + 15*x^16 + 28*x^14 + 237*x^12 + 579*x^10
            ....:                      + 1114*x^8 + 1470*x^6 + 2304*x^4 + 1296*x^2 + 729)
            sage: K.is_galois()                                                         # needs sage.groups
            True
            sage: K.discriminant().factor()
            2^20 * 3^10 * 53^10
            sage: K.decomposition_type(2)
            [(2, 5, 2)]
            sage: K.decomposition_type(3)
            [(2, 1, 10)]
            sage: K.decomposition_type(53)
            [(2, 2, 5)]

        This example is only ramified at 11::

            sage: K.<a> = NumberField(x^24 + 11^2*(90*x^12 - 640*x^8 + 2280*x^6
            ....:                                   - 512*x^4 + 2432/11*x^2 - 11))
            sage: K.discriminant().factor()
            -1 * 11^43
            sage: K.decomposition_type(11)
            [(1, 1, 2), (22, 1, 1)]

        Computing the decomposition type is feasible even in large degree::

            sage: K.<a> = NumberField(x^144 + 123*x^72 + 321*x^36 + 13*x^18 + 11)
            sage: K.discriminant().factor(limit=100000)
            2^144 * 3^288 * 7^18 * 11^17 * 31^18 * 157^18 * 2153^18 * 13907^18 * ...
            sage: K.decomposition_type(2)
            [(2, 4, 3), (2, 12, 2), (2, 36, 1)]
            sage: K.decomposition_type(3)
            [(9, 3, 2), (9, 10, 1)]
            sage: K.decomposition_type(7)
            [(1, 18, 1), (1, 90, 1), (2, 1, 6), (2, 3, 4)]

        It also works for relative extensions::

            sage: K.<a> = QuadraticField(-143)
            sage: M.<c> = K.extension(x^10 - 6*x^8 + (a + 12)*x^6 + (-7/2*a - 89/2)*x^4
            ....:                      + (13/2*a - 77/2)*x^2 + 25)

        There is a unique prime above `11` and above `13` in `K`, each of which is unramified in `M`::

            sage: M.decomposition_type(11)
            [(1, 2, 5)]
            sage: P11 = K.primes_above(11)[0]
            sage: len(M.primes_above(P11))
            5
            sage: M.decomposition_type(13)
            [(1, 1, 10)]
            sage: P13 = K.primes_above(13)[0]
            sage: len(M.primes_above(P13))
            10

        There are two primes above `2`, each of which ramifies in `M`::

            sage: Q0, Q1 = K.primes_above(2)
            sage: M.decomposition_type(Q0)
            [(2, 5, 1)]
            sage: q0, = M.primes_above(Q0)
            sage: q0.residue_class_degree()
            5
            sage: q0.relative_ramification_index()
            2
            sage: M.decomposition_type(Q1)
            [(2, 5, 1)]

        Check that :issue:`34514` is fixed::

            sage: K.<a> = NumberField(x^4 + 18*x^2 - 1)
            sage: R.<y> = K[]
            sage: L.<b> = K.extension(y^2 + 9*a^3 - 2*a^2 + 162*a - 38)
            sage: [L.decomposition_type(i) for i in K.primes_above(3)]
            [[(1, 1, 2)], [(1, 1, 2)], [(1, 2, 1)]]
        """
        v0 = self.base_ring().valuation(p)
        e0 = v0.value_group().gen().denominator()
        # Ideally we would compute f using the degree, but residue fields of relative extensions are currently implemented using polynomial quotient rings (this will hopefully be improved after #28485).
        C0 = v0.residue_field().cardinality()
        valuations = v0.extensions(self)
        ef = [(v.value_group().gen().denominator() // e0, v.residue_field().cardinality().exact_log(C0)) for v in valuations]
        return sorted([(e, f, g) for ((e, f), g) in Counter(ef).items()])

    def gen(self, n=0):
        """
        Return the generator for this number field.

        INPUT:

        - ``n`` -- must be 0 (the default), or an exception is raised

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
            raise IndexError("Only one generator.")
        try:
            return self.__gen
        except AttributeError:
            if self.__polynomial is not None:
                X = self.__polynomial.parent().gen()
            else:
                from sage.rings.polynomial.polynomial_ring import polygen
                X = polygen(QQ)
            self.__gen = self._element_class(self, X)
            return self.__gen

    @cached_method
    def _generator_matrix(self):
        """
        Return the matrix form of the generator of ``self``.

        .. SEEALSO::

            :meth:`~sage.rings.number_field.number_field_element.NumberFieldElement.matrix`

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: K.<v> = NumberField(x^4 + 514*x^2 + 64321)
            sage: R.<r> = NumberField(x^2 + 4*v*x + 5*v^2 + 514)
            sage: R._generator_matrix()
            [           0            1]
            [-5*v^2 - 514         -4*v]
        """
        x = self.gen()
        a = x
        d = self.relative_degree()
        v = x.list()
        for n in range(d-1):
            a *= x
            v += a.list()
        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(self.base_ring(), d)
        ret = M(v)
        ret.set_immutable()
        return ret

    def is_field(self, proof=True):
        """
        Return ``True`` since a number field is a field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^5 + x + 3, 'c').is_field()
            True
        """
        return True

    @cached_method
    def is_galois(self):
        r"""
        Return ``True`` if this number field is a Galois extension of
        `\QQ`.

        EXAMPLES::

            sage: # needs sage.groups
            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 1, 'i').is_galois()
            True
            sage: NumberField(x^3 + 2, 'a').is_galois()
            False
            sage: K = NumberField(x^15 + x^14 - 14*x^13 - 13*x^12 + 78*x^11 + 66*x^10
            ....:                  - 220*x^9 - 165*x^8 + 330*x^7 + 210*x^6 - 252*x^5
            ....:                  - 126*x^4 + 84*x^3 + 28*x^2 - 8*x - 1, 'a')
            sage: K.is_galois()
            True
            sage: K = NumberField(x^15 + x^14 - 14*x^13 - 13*x^12 + 78*x^11 + 66*x^10
            ....:                  - 220*x^9 - 165*x^8 + 330*x^7 + 210*x^6 - 252*x^5
            ....:                  - 126*x^4 + 84*x^3 + 28*x^2 - 8*x - 10, 'a')
            sage: K.is_galois()
            False
        """
        return self.galois_group().is_galois()

    @cached_method
    def is_abelian(self):
        r"""
        Return ``True`` if this number field is an abelian Galois extension of
        `\QQ`.

        EXAMPLES::

            sage: # needs sage.groups
            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 1, 'i').is_abelian()
            True
            sage: NumberField(x^3 + 2, 'a').is_abelian()
            False
            sage: NumberField(x^3 + x^2 - 2*x - 1, 'a').is_abelian()
            True
            sage: NumberField(x^6 + 40*x^3 + 1372, 'a').is_abelian()
            False
            sage: NumberField(x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1, 'a').is_abelian()
            True
        """
        if not self.is_galois():
            return False

        d = self.degree()
        dsqrt = d.isqrt()
        if d == 1 or d.is_prime() or (d == dsqrt**2 and dsqrt.is_prime()):
            return True

        if d <= 11:
            return self.galois_group().is_abelian()

        pari_pol = pari(self.polynomial())
        return pari_pol.galoisinit().galoisisabelian(1) == 1

    @cached_method
    def galois_group(self, type=None, algorithm='pari', names=None, gc_numbering=None):
        r"""
        Return the Galois group of the Galois closure of this number field.

        INPUT:

        - ``type`` -- deprecated; the different versions of Galois groups have been
          merged in :issue:`28782`

        - ``algorithm`` -- ``'pari'``, ``'gap'``, ``'kash'``, or ``'magma'``
          (default: ``'pari'``); for degrees between 12 and 15 default is
          ``'gap'``, and when the degree is >= 16 it is ``'kash'``)

        - ``names`` -- string giving a name for the generator of the Galois
          closure of ``self``, when this field is not Galois

        - ``gc_numbering`` -- if ``True``, permutations will be written
          in terms of the action on the roots of a defining polynomial
          for the Galois closure, rather than the defining polynomial for
          the original number field.  This is significantly faster;
          but not the standard way of presenting Galois groups.
          The default currently depends on the algorithm (``True`` for ``'pari'``,
          ``False`` for ``'magma'``) and may change in the future.

        The resulting group will only compute with automorphisms when necessary,
        so certain functions (such as :meth:`sage.rings.number_field.galois_group.GaloisGroup_v2.order`)
        will still be fast.  For more (important!)
        documentation, see the documentation for Galois groups of polynomials
        over `\QQ`, e.g., by typing ``K.polynomial().galois_group?``,
        where `K` is a number field.

        EXAMPLES::

            sage: # needs sage.groups
            sage: x = polygen(QQ, 'x')
            sage: k.<b> = NumberField(x^2 - 14)  # a Galois extension
            sage: G = k.galois_group(); G
            Galois group 2T1 (S2) with order 2 of x^2 - 14
            sage: G.gen(0)
            (1,2)
            sage: G.gen(0)(b)
            -b
            sage: G.artin_symbol(k.primes_above(3)[0])
            (1,2)

            sage: # needs sage.groups
            sage: k.<b> = NumberField(x^3 - x + 1)  # not Galois
            sage: G = k.galois_group(names='c'); G
            Galois group 3T2 (S3) with order 6 of x^3 - x + 1
            sage: G.gen(0)
            (1,2,3)(4,5,6)

            sage: NumberField(x^3 + 2*x + 1, 'a').galois_group(algorithm='magma')   # optional - magma, needs sage.groups
            Galois group Transitive group number 2 of degree 3
             of the Number Field in a with defining polynomial x^3 + 2*x + 1

        EXPLICIT GALOIS GROUP: We compute the Galois group as an explicit
        group of automorphisms of the Galois closure of a field.

        ::

            sage: # needs sage.groups
            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<b1> = K.galois_closure(); L
            Number Field in b1 with defining polynomial x^6 + 108
            sage: G = End(L); G
            Automorphism group of Number Field in b1 with defining polynomial x^6 + 108
            sage: G.list()
            [Ring endomorphism of Number Field in b1 with defining polynomial x^6 + 108
               Defn: b1 |--> b1,
             ...
             Ring endomorphism of Number Field in b1 with defining polynomial x^6 + 108
               Defn: b1 |--> -1/12*b1^4 - 1/2*b1]
            sage: G[2](b1)
            1/12*b1^4 + 1/2*b1

        Many examples for higher degrees may be found in the online databases
        http://galoisdb.math.upb.de/ by Jürgen Klüners and Gunter Malle and
        https://www.lmfdb.org/NumberField/ by the LMFDB collaboration,
        although these might need a lot of computing time.

        If `L/K` is a relative number field, this method will currently return `Gal(L/\QQ)`.  This behavior will
        change in the future, so it's better to explicitly call :meth:`absolute_field` if that is
        the desired behavior::

            sage: # needs sage.groups
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 + 1)
            sage: R.<t> = PolynomialRing(K)
            sage: L = K.extension(t^5 - t + a, 'b')
            sage: L.galois_group()
            ...DeprecationWarning: Use .absolute_field().galois_group()
            if you want the Galois group of the absolute field
            See https://github.com/sagemath/sage/issues/28782 for details.
            Galois group 10T22 (S(5)[x]2) with order 240 of t^5 - t + a

        TESTS:

        We check that the changes in :issue:`28782` won't break code that used v1 Galois groups::

            sage: # needs sage.groups
            sage: G = NumberField(x^3 - 2, 'a').galois_group(type='pari')
            ...DeprecationWarning: the different Galois types have been merged into one class
            See https://github.com/sagemath/sage/issues/28782 for details.
            sage: G.group()
            ...DeprecationWarning: the group method is deprecated; you can use _pol_galgp if you really need it
            See https://github.com/sagemath/sage/issues/28782 for details.
            PARI group [6, -1, 2, "S3"] of degree 3
        """
        if type is not None:
            deprecation(28782, "the different Galois types have been merged into one class")

        from .galois_group import GaloisGroup_v2
        return GaloisGroup_v2(self, algorithm=algorithm, names=names, gc_numbering=gc_numbering, _type=type)

    def _normalize_prime_list(self, v):
        """
        Internal function to convert into a tuple of primes either ``None`` or
        a single prime or a list.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K._normalize_prime_list(None)
            ()
            sage: K._normalize_prime_list(3)
            (3,)
            sage: K._normalize_prime_list([3,5])
            (3, 5)
            sage: K._normalize_prime_list([5, 3, 3])
            (3, 5)
            sage: K._normalize_prime_list([3, 5, -3])
            (3, 5)
        """
        if v is None:
            v = []
        elif not isinstance(v, (list, tuple)):
            v = [v]

        v = [ZZ(p).abs() for p in v]
        v = sorted(set(v))

        return tuple(v)

    def power_basis(self):
        r"""
        Return a power basis for this number field over its base field.

        If this number field is represented as `k[t]/f(t)`, then
        the basis returned is `1, t, t^2, \ldots, t^{d-1}` where
        `d` is the degree of this number field over its base
        field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        return [g**i for i in range(self.relative_degree())]

    def integral_basis(self, v=None):
        """
        Return a list containing a ``ZZ``-basis for the full ring of integers
        of this number field.

        INPUT:

        - ``v`` -- ``None``, a prime, or a list of primes; see the
          documentation for :meth:`maximal_order`

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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

        ALGORITHM: Uses the PARI library (via :pari:`_pari_integral_basis`).
        """
        return self.maximal_order(v=v).basis()

    def _pari_integral_basis(self, v=None, important=True):
        """
        Internal function returning an integral basis of this number field in
        PARI format.

        INPUT:

        - ``v`` -- ``None``, a prime, or a list of primes; see the
          documentation for :meth:``maximal_order``

        - ``important`` -- boolean (default: ``True``); if ``False``,
          raise a :exc:`RuntimeError` if we need to do a difficult
          discriminant factorization.  This is useful when an integral
          basis is not strictly required.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
            pass

        f = self.pari_polynomial("y")
        if v:
            # NOTE: here we make pari know about potentially big primes factors of
            # the discriminant, see
            # https://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=2257
            primelimit = pari.default("primelimit")
            primes = [p for p in v if p > primelimit]
            if primes:
                pari.addprimes(primes)
            B = f.nfbasis(fa=v)
        elif self._assume_disc_small:
            B = f.nfbasis(1)
        elif not important:
            # Trial divide the discriminant with primes up to 10^6
            m = self.pari_polynomial().poldisc().abs().factor(limit=10**6)
            # Since we only need a *squarefree* factorization for
            # primes with exponent 1, we need trial division up to D^(1/3)
            # instead of D^(1/2).
            trialdivlimit2 = pari(10**12)
            trialdivlimit3 = pari(10**18)
            if all(p < trialdivlimit2 or (e == 1 and p < trialdivlimit3) or p.isprime() for p, e in zip(m[0], m[1])):
                B = f.nfbasis(fa=m)
            else:
                raise RuntimeError("Unable to factor discriminant with trial division")
        else:
            B = f.nfbasis()

        self._integral_basis_dict[v] = B
        return B

    def reduced_basis(self, prec=None):
        r"""
        Return an LLL-reduced basis for the Minkowski-embedding
        of the maximal order of a number field.

        INPUT:

        - ``prec`` -- (default: ``None``) the precision with which to
          compute the Minkowski embedding

        OUTPUT:

        An LLL-reduced basis for the Minkowski-embedding of the
        maximal order of a number field, given by a sequence of (integral)
        elements from the field.

        .. NOTE::

            In the non-totally-real case, the LLL routine we call is
            currently PARI's :pari:`qflll`, which works with floating point
            approximations, and so the result is only as good as the
            precision promised by PARI. The matrix returned will always
            be integral; however, it may only be only "almost" LLL-reduced
            when the precision is not sufficiently high.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<t> = NumberField(x^6 - 7*x^4 - x^3 + 11*x^2 + x - 1)
            sage: F.maximal_order().basis()
            [1/2*t^5 + 1/2*t^4 + 1/2*t^2 + 1/2, t, t^2, t^3, t^4, t^5]
            sage: F.reduced_basis()
            [-1, -1/2*t^5 + 1/2*t^4 + 3*t^3 - 3/2*t^2 - 4*t - 1/2, t,
             1/2*t^5 + 1/2*t^4 - 4*t^3 - 5/2*t^2 + 7*t + 1/2,
             1/2*t^5 - 1/2*t^4 - 2*t^3 + 3/2*t^2 - 1/2,
             1/2*t^5 - 1/2*t^4 - 3*t^3 + 5/2*t^2 + 4*t - 5/2]
            sage: CyclotomicField(12).reduced_basis()
            [1, zeta12^2, zeta12, zeta12^3]

        TESTS:

        Check that the bug reported at :issue:`10017` is fixed::

            sage: x = polygen(QQ)
            sage: k.<a> = NumberField(x^6 + 2218926655879913714112*x^4 - 32507675650290949030789018433536*x^3 + 4923635504174417014460581055002374467948544*x^2 - 36066074010564497464129951249279114076897746988630016*x + 264187244046129768986806800244258952598300346857154900812365824)
            sage: new_basis = k.reduced_basis(prec=120)
            sage: [c.minpoly() for c in new_basis]
            [x - 1,
             x^2 + x + 1,
             x^6 + 3*x^5 - 102*x^4 - 103*x^3 + 10572*x^2 - 59919*x + 127657,
             x^6 + 3*x^5 - 102*x^4 - 103*x^3 + 10572*x^2 - 59919*x + 127657,
             x^3 - 171*x + 848,
             x^6 + 171*x^4 + 1696*x^3 + 29241*x^2 + 145008*x + 719104]
            sage: R = k.order(new_basis)
            sage: R.discriminant()==k.discriminant()
            True
        """
        ZK = self.integral_basis()
        d = self.absolute_degree()

        # If self is totally real, then we can use (x*y).trace() as
        # the inner product on the Minkowski embedding, which is
        # faster than computing all the conjugates, etc ...

        # flag to disable FLATTER, which is much more unstable than fplll
        flag = 1 if pari.version() >= (2,17) else 0
        if self.is_totally_real():
            from sage.matrix.constructor import matrix
            M = matrix(ZZ, d, d, [[(x*y).trace() for x in ZK] for y in ZK])
            T = pari(M).qflllgram(flag=flag)
        else:
            M = self.minkowski_embedding(ZK, prec=prec)
            T = pari(M).qflll(flag=flag)

        return [sum([ZZ(T[i][j]) * ZK[j] for j in range(d)]) for i in range(d)]

    def reduced_gram_matrix(self, prec=None):
        r"""
        Return the Gram matrix of an LLL-reduced basis for
        the Minkowski embedding of the maximal order of a number field.

        INPUT:

        - ``prec`` -- (default: ``None``) the precision with which
          to calculate the Minkowski embedding (see NOTE below)

        OUTPUT: the Gram matrix `[\langle x_i,x_j \rangle]` of an LLL reduced
        basis for the maximal order of ``self``, where the integral basis for
        ``self`` is given by `\{x_0, \dots, x_{n-1}\}`. Here `\langle , \rangle` is
        the usual inner product on `\RR^n`, and ``self`` is embedded in `\RR^n` by
        the Minkowski embedding. See the docstring for
        :meth:`NumberField_absolute.minkowski_embedding` for more information.

        .. NOTE::

           In the non-totally-real case, the LLL routine we call is
           currently PARI's :pari:`qflll`, which works with floating point
           approximations, and so the result is only as good as the
           precision promised by PARI. In particular, in this case,
           the returned matrix will *not* be integral, and may not
           have enough precision to recover the correct Gram matrix
           (which is known to be integral for theoretical
           reasons). Thus the need for the ``prec`` parameter above.

        If the following run-time error occurs: "PariError: not a definite
        matrix in lllgram (42)", try increasing the ``prec`` parameter,

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<t> = NumberField(x^6 - 7*x^4 - x^3 + 11*x^2 + x - 1)
            sage: F.reduced_gram_matrix()
            [ 6  3  0  2  0  1]
            [ 3  9  0  1  0 -2]
            [ 0  0 14  6 -2  3]
            [ 2  1  6 16 -3  3]
            [ 0  0 -2 -3 16  6]
            [ 1 -2  3  3  6 19]
            sage: Matrix(6, [(x*y).trace()
            ....:            for x in F.integral_basis() for y in F.integral_basis()])
            [2550  133  259  664 1368 3421]
            [ 133   14    3   54   30  233]
            [ 259    3   54   30  233  217]
            [ 664   54   30  233  217 1078]
            [1368   30  233  217 1078 1371]
            [3421  233  217 1078 1371 5224]

        ::

            sage: x = polygen(QQ)
            sage: F.<alpha> = NumberField(x^4 + x^2 + 712312*x + 131001238)
            sage: F.reduced_gram_matrix(prec=128)
            [   4.0000000000000000000000000000000000000   0.00000000000000000000000000000000000000   -1.9999999999999999999999999999999999037  -0.99999999999999999999999999999999383702]
            [  0.00000000000000000000000000000000000000    46721.539331563218381658483353092335550   -11488.910026551724275122749703614966768   -418.12718083977141198754424579680468382]
            [  -1.9999999999999999999999999999999999037   -11488.910026551724275122749703614966768  5.5658915310500611768713076521847709187e8  1.4179092271494070050433368847682152174e8]
            [ -0.99999999999999999999999999999999383702   -418.12718083977141198754424579680468382  1.4179092271494070050433368847682152174e8 1.3665897267919181137884111201405279175e12]
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
                                                [[(x * y).trace() for x in B]
                                                 for y in B])
        else:
            M = self.minkowski_embedding(prec=prec)
            T = matrix(d, flatten([a.vector().list()
                                   for a in self.reduced_basis(prec=prec)]))
            A = M * T.transpose()
            self.__reduced_gram_matrix = A.transpose() * A
            if prec is None:
                # this is the default choice for minkowski_embedding
                self.__reduced_gram_matrix_prec = 53
            else:
                self.__reduced_gram_matrix_prec = prec

        return self.__reduced_gram_matrix

    # ******************************************************
    #  Supplementary algorithm to enumerate lattice points
    # ******************************************************

    def _positive_integral_elements_with_trace(self, C):
        r"""
        Find all totally positive integral elements in ``self`` whose
        trace is between C[0] and C[1], inclusive.

        .. NOTE::

           This is currently only implemented in the case that ``self`` is
           totally real, since it requires exact computation of
           :meth:`.reduced_gram_matrix`.

        EXAMPLES::

            sage: K.<alpha> = NumberField(ZZ['x'].0^2 - 2)
            sage: K._positive_integral_elements_with_trace([0,5])
            [alpha + 2, -alpha + 2, 2, 1]
            sage: L.<beta> = NumberField(ZZ['x'].0^2 + 1)
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
            raise ValueError("bounds must be positive")

        if not self.is_totally_real():
            raise NotImplementedError("exact computation of LLL reduction only implemented in the totally real case")

        B = self.reduced_basis()
        T = self.reduced_gram_matrix()
        P = pari(T).qfminim((C[1]**2) * 0.5, 10**6)[2]

        S = []
        for p in P:
            theta = sum([p.list()[i] * B[i] for i in range(self.degree())])
            if theta.trace() < 0:
                theta *= -1
            if C[0] <= theta.trace() <= C[1]:
                if self(theta).is_totally_positive():
                    S.append(self(theta))
        return S

    @cached_method
    def narrow_class_group(self, proof=None):
        r"""
        Return the narrow class group of this field.

        INPUT:

        - ``proof`` -- (default: ``None``) use the global proof setting, which
          defaults to ``True``

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^3 + x + 9, 'a').narrow_class_group()
            Multiplicative Abelian group isomorphic to C2

        TESTS::

            sage: QuadraticField(3, 'a').narrow_class_group()
            Multiplicative Abelian group isomorphic to C2
        """
        proof = proof_flag(proof)
        k = self.pari_bnf(proof)
        s = k.bnfnarrow().sage()
        return sage.groups.abelian_gps.abelian_group.AbelianGroup(s[1])

    def ngens(self):
        """
        Return the number of generators of this number field (always 1).

        OUTPUT: the python integer 1

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 19,'a').order()
            +Infinity
        """
        return infinity.infinity

    def absolute_polynomial_ntl(self):
        r"""
        Alias for :meth:`~polynomial_ntl`. Mostly for internal use.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').absolute_polynomial_ntl()
            ([-27 34 51], 51)
        """
        return self.polynomial_ntl()

    def polynomial_ntl(self):
        r"""
        Return defining polynomial of this number field as a pair, an ntl
        polynomial and a denominator.

        This is used mainly to implement some internal arithmetic.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        r"""
        Return the defining polynomial of this number field.

        This is exactly the same as :meth:`defining_polynomial`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').polynomial()
            x^2 + 2/3*x - 9/17
        """
        return self.__polynomial

    def defining_polynomial(self):   # do not overload this -- overload polynomial instead
        r"""
        Return the defining polynomial of this number field.

        This is exactly the same as :meth:`polynomial`.

        EXAMPLES::

            sage: k5.<z> = CyclotomicField(5)
            sage: k5.defining_polynomial()
            x^4 + x^3 + x^2 + x + 1
            sage: y = polygen(QQ, 'y')
            sage: k.<a> = NumberField(y^9 - 3*y + 5); k
            Number Field in a with defining polynomial y^9 - 3*y + 5
            sage: k.defining_polynomial()
            y^9 - 3*y + 5
        """
        return self.polynomial()

    def polynomial_ring(self):
        r"""
        Return the polynomial ring that we view this number field as being
        a quotient of (by a principal ideal).

        EXAMPLES: An example with an absolute field::

            sage: x = polygen(QQ, 'x')
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
            Univariate Polynomial Ring in y over
             Number Field in a1 with defining polynomial y^2 + 1
        """
        return self.relative_polynomial().parent()

    def polynomial_quotient_ring(self):
        r"""
        Return the polynomial quotient ring isomorphic to this number
        field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K = NumberField(x^3 + 2*x - 5, 'alpha')
            sage: K.polynomial_quotient_ring()
            Univariate Quotient Polynomial Ring in alpha over Rational Field
             with modulus x^3 + 2*x - 5
        """
        return self.polynomial_ring().quotient(self.relative_polynomial(), self.variable_name())

    def regulator(self, proof=None):
        r"""
        Return the regulator of this number field.

        Note that PARI computes the regulator to higher precision than the
        Sage default.

        INPUT:

        - ``proof`` -- (default: ``True``) unless you set it otherwise

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 - 2, 'a').regulator()
            0.881373587019543
            sage: NumberField(x^4 + x^3 + x^2 + x + 1, 'a').regulator()
            0.962423650119207
        """
        proof = proof_flag(proof)
        try:
            return self.__regulator
        except AttributeError:
            from sage.rings.real_mpfr import RealField
            k = self.pari_bnf(proof)
            self.__regulator = RealField(53)(k.bnf_get_reg())
            return self.__regulator

    def residue_field(self, prime, names=None, check=True):
        r"""
        Return the residue field of this number field at a given prime, ie
        `O_K / p O_K`.

        INPUT:

        - ``prime`` -- a prime ideal of the maximal order in
          this number field, or an element of the field which generates a
          principal prime ideal.

        - ``names`` -- the name of the variable in the residue field

        - ``check`` -- whether or not to check the primality of ``prime``

        OUTPUT: the residue field at this prime

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4 + 3*x^2 - 17)
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
        if isinstance(prime, NumberFieldIdeal) and prime.number_field() is not self:
            raise ValueError("%s is not an ideal of %s" % (prime, self))
        # This allows principal ideals to be specified using a generator:
        try:
            prime = self.ideal(prime)
        except TypeError:
            pass

        if not isinstance(prime, NumberFieldIdeal) or prime.number_field() is not self:
            raise ValueError("%s is not an ideal of %s" % (prime, self))
        if check and not prime.is_prime():
            raise ValueError("%s is not a prime ideal" % prime)
        from sage.rings.finite_rings.residue_field import ResidueField
        return ResidueField(prime, names=names, check=False)

    def signature(self):
        r"""
        Return `(r_1, r_2)`, where `r_1` and `r_2` are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: NumberField(x^2 + 1, 'a').signature()
            (0, 1)
            sage: NumberField(x^3 - 2, 'a').signature()
            (1, 1)
        """
        r1, r2 = self.pari_nf().nf_get_sign()
        return (ZZ(r1), ZZ(r2))

    def trace_pairing(self, v):
        r"""
        Return the matrix of the trace pairing on the elements of the list
        ``v``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<zeta3> = NumberField(x^2 + 3)
            sage: K.trace_pairing([1, zeta3])
            [ 2  0]
            [ 0 -6]
        """
        import sage.matrix.matrix_space
        A = sage.matrix.matrix_space.MatrixSpace(self.base_ring(), len(v))(0)
        for i in range(len(v)):
            for j in range(i, len(v)):
                t = (self(v[i]*v[j])).trace()
                A[i, j] = t
                A[j, i] = t
        return A

    def uniformizer(self, P, others='positive'):
        """
        Return an element of ``self`` with valuation 1 at the prime ideal `P`.

        INPUT:

        - ``self`` -- a number field

        - ``P`` -- a prime ideal of ``self``

        - ``others`` -- either ``'positive'`` (default), in which case the
          element will have nonnegative valuation at all other primes of
          ``self``, or ``'negative'``, in which case the element will have
          nonpositive valuation at all other primes of ``self``

        .. NOTE::

           When `P` is principal (e.g., always when ``self`` has class number
           one) the result may or may not be a generator of `P`!

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 5); K
            Number Field in a with defining polynomial x^2 + 5
            sage: P, Q = K.ideal(3).prime_factors()
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
            sage: Plist = K.ideal(17).prime_factors()
            sage: pilist = [K.uniformizer(P) for P in Plist]
            sage: [pi.is_integral() for pi in pilist]
            [True, True, True]
            sage: [pi.valuation(P) for pi, P in zip(pilist, Plist)]
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
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(7))]  # representation varies, not tested
            [t^2 + 3*t + 1]
            sage: [K.uniformizer(P) for P,e in factor(K.ideal(67))]
            [t + 23, t + 26, t - 32, t - 18]

        ALGORITHM:

            Use PARI. More precisely, use the second component of
            :pari:`idealprimedec` in the "positive" case. Use :pari:`idealappr`
            with exponent of `-1` and invert the result in the "negative"
            case.

        TESTS:

        We test the above doctest. The representation depends on the PARI version::

            sage: K.<t> = NumberField(x^4 - x^3 - 3*x^2 - x + 1)
            sage: [x] = [K.uniformizer(P) for P,e in factor(K.ideal(7))]
            sage: x in (t^2 + 3*t +1, t^2 - 4*t +1)
            True
        """
        if not isinstance(P, NumberFieldIdeal):
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

        ALGORITHM: Uses PARI's :pari:`bnfinit` command.

        INPUT:

        - ``proof`` -- boolean (default: ``True``); flag passed to PARI

        .. NOTE::

            For more functionality see :meth:`unit_group`.

        .. SEEALSO::

            :meth:`unit_group`
            :meth:`S_unit_group`
            :meth:`S_units`

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: K.units()
            (-1/275*a^3 - 4/55*a^2 + 5/11*a - 3,)

        For big number fields, provably computing the unit group can
        take a very long time.  In this case, one can ask for the
        conjectural unit group (correct if the Generalized Riemann
        Hypothesis is true)::

            sage: K = NumberField(x^17 + 3, 'a')
            sage: K.units(proof=True)  # takes forever, not tested
            ...
            sage: K.units(proof=False)  # result not independently verified
            (a^9 + a - 1,
             -a^15 + a^12 - a^10 + a^9 + 2*a^8 - 3*a^7 - a^6 + 3*a^5 - a^4 - 4*a^3 + 3*a^2 + 2*a - 2,
             a^15 + a^14 + a^13 + a^12 + a^10 - a^7 - a^6 - a^2 - 1,
             2*a^16 - 3*a^15 + 3*a^14 - 3*a^13 + 3*a^12 - a^11 + a^9 - 3*a^8 + 4*a^7 - 5*a^6 + 6*a^5 - 4*a^4 + 3*a^3 - 2*a^2 - 2*a + 4,
             -a^16 + a^15 - a^14 + a^12 - a^11 + a^10 + a^8 - a^7 + 2*a^6 - a^4 + 3*a^3 - 2*a^2 + 2*a - 1,
             a^16 - 2*a^15 - 2*a^13 - a^12 - a^11 - 2*a^10 + a^9 - 2*a^8 + 2*a^7 - 3*a^6 - 3*a^4 - 2*a^3 - a^2 - 4*a + 2,
             -a^15 - a^14 - 2*a^11 - a^10 + a^9 - a^8 - 2*a^7 + a^5 - 2*a^3 + a^2 + 3*a - 1,
             -3*a^16 - 3*a^15 - 3*a^14 - 3*a^13 - 3*a^12 - 2*a^11 - 2*a^10 - 2*a^9 - a^8 + a^7 + 2*a^6 + 3*a^5 + 3*a^4 + 4*a^3 + 6*a^2 + 8*a + 8)

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(1/2*x^2 - 1/6)
            sage: K.units()
            (3*a + 2,)
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

        # get PARI to compute the fundamental units
        B = self.pari_bnf(proof).bnf_get_fu()
        B = tuple(self(b, check=False) for b in B)
        if proof:
            # cache the provable results and return them
            self.__units = B
            return self.__units
        else:
            # cache the conjectural results and return them
            self.__units_no_proof = B
            return self.__units_no_proof

    def unit_group(self, proof=None):
        """
        Return the unit group (including torsion) of this number field.

        ALGORITHM: Uses PARI's :pari:`bnfinit` command.

        INPUT:

        - ``proof`` -- boolean (default: ``True``); flag passed to PARI

        .. NOTE::

           The group is cached.

        .. SEEALSO::

            :meth:`units`
            :meth:`S_unit_group`
            :meth:`S_units`

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: U = K.unit_group(); U
            Unit group with structure C10 x Z of Number Field in a
             with defining polynomial x^4 - 10*x^3 + 100*x^2 - 375*x + 1375
            sage: U.gens()
            (u0, u1)
            sage: U.gens_values()  # random
            [-1/275*a^3 + 7/55*a^2 - 6/11*a + 4, 1/275*a^3 + 4/55*a^2 - 5/11*a + 3]
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
            Unit group with structure C2 x Z x Z x Z x Z x Z x Z x Z x Z of
             Number Field in a with defining polynomial x^17 + 3
            sage: U.gens()
            (u0, u1, u2, u3, u4, u5, u6, u7, u8)
            sage: U.gens_values()  # result not independently verified
            [-1,
             a^9 + a - 1,
             -a^15 + a^12 - a^10 + a^9 + 2*a^8 - 3*a^7 - a^6 + 3*a^5 - a^4 - 4*a^3 + 3*a^2 + 2*a - 2,
             a^15 + a^14 + a^13 + a^12 + a^10 - a^7 - a^6 - a^2 - 1,
             2*a^16 - 3*a^15 + 3*a^14 - 3*a^13 + 3*a^12 - a^11 + a^9 - 3*a^8 + 4*a^7 - 5*a^6 + 6*a^5 - 4*a^4 + 3*a^3 - 2*a^2 - 2*a + 4,
             -a^16 + a^15 - a^14 + a^12 - a^11 + a^10 + a^8 - a^7 + 2*a^6 - a^4 + 3*a^3 - 2*a^2 + 2*a - 1,
             a^16 - 2*a^15 - 2*a^13 - a^12 - a^11 - 2*a^10 + a^9 - 2*a^8 + 2*a^7 - 3*a^6 - 3*a^4 - 2*a^3 - a^2 - 4*a + 2,
             -a^15 - a^14 - 2*a^11 - a^10 + a^9 - a^8 - 2*a^7 + a^5 - 2*a^3 + a^2 + 3*a - 1,
             -3*a^16 - 3*a^15 - 3*a^14 - 3*a^13 - 3*a^12 - 2*a^11 - 2*a^10 - 2*a^9 - a^8 + a^7 + 2*a^6 + 3*a^5 + 3*a^4 + 4*a^3 + 6*a^2 + 8*a + 8]
        """
        proof = proof_flag(proof)

        try:
            return self._unit_group
        except AttributeError:
            pass

        if not proof:
            try:
                return self._unit_group_no_proof
            except AttributeError:
                pass

        U = UnitGroup(self, proof)
        if proof:
            self._unit_group = U
        else:
            self._unit_group_no_proof = U
        return U

    def S_unit_group(self, proof=None, S=None):
        """
        Return the `S`-unit group (including torsion) of this number field.

        ALGORITHM: Uses PARI's :pari:`bnfsunit` command.

        INPUT:

        - ``proof`` -- boolean (default: ``True``); flag passed to PARI

        - ``S`` -- list or tuple of prime ideals, or an ideal, or a single
          ideal or element from which an ideal can be constructed, in
          which case the support is used.  If ``None``, the global unit
          group is constructed; otherwise, the `S`-unit group is
          constructed.

        .. NOTE::

           The group is cached.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3)
            sage: U = K.S_unit_group(S=a); U
            S-unit group with structure C10 x Z x Z x Z of
             Number Field in a with defining polynomial x^4 - 10*x^3 + 100*x^2 - 375*x + 1375
             with S = (Fractional ideal (5, -7/275*a^3 + 1/11*a^2 - 9/11*a),
                       Fractional ideal (11, -7/275*a^3 + 1/11*a^2 - 9/11*a + 3))
            sage: U.gens()
            (u0, u1, u2, u3)
            sage: U.gens_values()  # random
            [-1/275*a^3 + 7/55*a^2 - 6/11*a + 4, 1/275*a^3 + 4/55*a^2 - 5/11*a + 3,
             1/275*a^3 + 4/55*a^2 - 5/11*a + 5, -14/275*a^3 + 21/55*a^2 - 29/11*a + 6]
            sage: U.invariants()
            (10, 0, 0, 0)
            sage: [u.multiplicative_order() for u in U.gens()]
            [10, +Infinity, +Infinity, +Infinity]
            sage: U.primes()
            (Fractional ideal (5, -7/275*a^3 + 1/11*a^2 - 9/11*a),
             Fractional ideal (11, -7/275*a^3 + 1/11*a^2 - 9/11*a + 3))

        With the default value of `S`, the S-unit group is the same as
        the global unit group::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.unit_group(proof=False)
            sage: U.is_isomorphic(K.S_unit_group(proof=False))
            True

        The value of `S` may be specified as a list of prime ideals,
        or an ideal, or an element of the field::

            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.S_unit_group(proof=False, S=K.ideal(6).prime_factors()); U
            S-unit group with structure C2 x Z x Z x Z x Z
             of Number Field in a with defining polynomial x^3 + 3
             with S = (Fractional ideal (-a^2 + a - 1),
                       Fractional ideal (a + 1),
                       Fractional ideal (a))
            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.S_unit_group(proof=False, S=K.ideal(6)); U
            S-unit group with structure C2 x Z x Z x Z x Z
             of Number Field in a with defining polynomial x^3 + 3
             with S = (Fractional ideal (-a^2 + a - 1),
                       Fractional ideal (a + 1),
                       Fractional ideal (a))
            sage: K.<a> = NumberField(x^3 + 3)
            sage: U = K.S_unit_group(proof=False, S=6); U
            S-unit group with structure C2 x Z x Z x Z x Z
             of Number Field in a with defining polynomial x^3 + 3
             with S = (Fractional ideal (-a^2 + a - 1),
                       Fractional ideal (a + 1),
                       Fractional ideal (a))
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
            if isinstance(S, list):
                S = tuple(S)
            if not isinstance(S, tuple):
                try:
                    S = tuple(self.ideal(S).prime_factors())
                except (NameError, TypeError, ValueError):
                    raise ValueError(f"Cannot make a set of primes from {S}")
            else:
                try:
                    S = tuple(self.ideal(P) for P in S)
                except (NameError, TypeError, ValueError):
                    raise ValueError(f"Cannot make a set of primes from {S}")
                if not all(P.is_prime() for P in S):
                    raise ValueError(f"Not all elements of {S} are prime ideals")

        try:
            return self._S_unit_group_cache[S]
        except AttributeError:
            self._S_unit_group_cache = {}
        except KeyError:
            pass

        if not proof:
            try:
                return self._S_unit_group_no_proof_cache[S]
            except AttributeError:
                self._S_unit_group_no_proof_cache = {}
            except KeyError:
                pass

        U = UnitGroup(self, proof, S=S)
        if proof:
            self._S_unit_group_cache[S] = U
        else:
            self._S_unit_group_no_proof_cache[S] = U
        return U

    def S_unit_solutions(self, S=[], prec=106, include_exponents=False, include_bound=False, proof=None):
        r"""
        Return all solutions to the `S`-unit equation `x + y = 1` over ``self``.

        INPUT:

        - ``S`` -- list of finite primes in this number field
        - ``prec`` -- precision used for computations in real, complex, and
          `p`-adic fields (default: 106)
        - ``include_exponents`` -- whether to include the exponent vectors in
          the returned value (default: ``True``)
        - ``include_bound`` -- whether to return the final computed bound
          (default: ``False``)
        - ``proof`` -- if ``False``, assume the GRH in computing the class group;
          default is ``True``

        OUTPUT:

        A list `[(A_1, B_1, x_1, y_1), (A_2, B_2, x_2, y_2), \dots, (A_n, B_n, x_n, y_n)]` of tuples such that:

        1. The first two entries are tuples `A_i = (a_0, a_1, \dots, a_t)` and `B_i = (b_0, b_1, \dots, b_t)` of exponents.
           These will be omitted if ``include_exponents`` is ``False``.

        2. The last two entries are `S`-units `x_i` and `y_i` in ``self`` with `x_i + y_i = 1`.

        3. If the default generators for the `S`-units of ``self`` are `(\rho_0, \rho_1, \dots, \rho_t)``,
           then these satisfy `x_i = \prod(\rho_i)^{(a_i)}` and `y_i = \prod(\rho_i)^{(b_i)}`.

        If ``include_bound`` is ``True``, will return a pair ``(sols, bound)`` where ``sols`` is as above
        and ``bound`` is the bound used for the entries in the exponent vectors.

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: x = polygen(QQ, 'x')
            sage: K.<xi> = NumberField(x^2 + x + 1)
            sage: S = K.primes_above(3)
            sage: K.S_unit_solutions(S)  # random, due to ordering
            [(xi + 2, -xi - 1), (1/3*xi + 2/3, -1/3*xi + 1/3), (-xi, xi + 1), (-xi + 1, xi)]

        You can get the exponent vectors::

            sage: # needs sage.rings.padics
            sage: K.S_unit_solutions(S, include_exponents=True)  # random, due to ordering
            [((2, 1), (4, 0), xi + 2, -xi - 1),
             ((5, -1), (4, -1), 1/3*xi + 2/3, -1/3*xi + 1/3),
             ((5, 0), (1, 0), -xi, xi + 1),
             ((1, 1), (2, 0), -xi + 1, xi)]

        And the computed bound::

            sage: # needs sage.rings.padics
            sage: solutions, bound = K.S_unit_solutions(S, prec=100, include_bound=True)
            sage: bound
            6
        """
        from .S_unit_solver import solve_S_unit_equation
        return solve_S_unit_equation(self, S, prec, include_exponents, include_bound, proof)

    def zeta(self, n=2, all=False):
        """
        Return one, or a list of all, primitive `n`-th root of unity in this field.

        INPUT:

        - ``n`` -- positive integer

        - ``all`` -- boolean; if ``False`` (default), return a primitive
          `n`-th root of unity in this field, or raise a :exc:`ValueError`
          exception if there are none.  If ``True``, return a list of all
          primitive `n`-th roots of unity in this field (possibly empty).

        .. NOTE::

            To obtain the maximal order of a root of unity in this field,
            use :meth:`number_of_roots_of_unity`.

        .. NOTE::

            We do not create the full unit group since that can be
            expensive, but we do use it if it is already known.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
            ValueError: there are no 4th roots of unity in self

        ::

            sage: r.<x> = QQ[]
            sage: K.<b> = NumberField(x^2 + 1)
            sage: K.zeta(4)
            b
            sage: K.zeta(4,all=True)
            [b, -b]
            sage: K.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: there are no 3rd roots of unity in self
            sage: K.zeta(3, all=True)
            []

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(1/2*x^2 + 1/6)
            sage: K.zeta(3)
            -3/2*a - 1/2

        TESTS::

            sage: K = NumberField(x**60+691*x**12-25,'a')
            sage: K.zeta(15,all=True)
            []
        """
        try:
            return self._unit_group.zeta(n, all)
        except AttributeError:
            pass
        try:
            return self._unit_group_no_proof.zeta(n, all)
        except AttributeError:
            pass

        K = self
        n = ZZ(n)
        if n <= 0:
            raise ValueError("n (=%s) must be positive" % n)
        if n == 1:
            if all:
                return [K.one()]
            return K.one()
        elif n == 2:
            if all:
                return [K(-1)]
            return K(-1)

        # First check if the degree of K is compatible
        # with an inclusion QQ(\zeta_n) -> K.
        if euler_phi(n).divides(K.absolute_degree()):
            w, zeta_w = self.pari_nf().nfrootsof1()
            w = w.sage()
            zeta_w = K(zeta_w)
            if not w % n:
                zeta_n = zeta_w**(w // n)
                if all:
                    return [zeta_n**i for i in n.coprime_integers(n)]
                return zeta_n
        if all:
            return []
        raise ValueError("there are no %s roots of unity in self" % n.ordinal_str())

    def zeta_order(self):
        r"""
        Return the number of roots of unity in this field.

        .. NOTE::

           We do not create the full unit group since that can be
           expensive, but we do use it if it is already known.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^22 + 3)
            sage: F.zeta_order()
            6
            sage: F.<alpha> = NumberField(x^2 - 7)
            sage: F.zeta_order()
            2

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(1/2*x^2 + 1/6)
            sage: K.zeta_order()
            6
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
        cyclotomic fields. Repeated calls of this function may return
        a different value.

        .. NOTE::

           We do not create the full unit group since that can be
           expensive, but we do use it if it is already known.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: z = K.primitive_root_of_unity(); z
            i
            sage: z.multiplicative_order()
            4

            sage: K.<a> = NumberField(x^2 + x + 1)
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

            sage: K = CyclotomicField(3)
            sage: z = K.primitive_root_of_unity(); z
            zeta3 + 1
            sage: z.multiplicative_order()
            6

        TESTS:

        Check for :issue:`15027`. We use a new variable name::

            sage: K.<f> = NumberField(x^2 + x + 1)
            sage: K.primitive_root_of_unity()
            f + 1
            sage: UK = K.unit_group()
            sage: K.primitive_root_of_unity()
            f + 1

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(3*x^2 + 1)
            sage: K.primitive_root_of_unity()
            -3/2*a + 1/2
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
        return self(z, check=False)

    def roots_of_unity(self):
        """
        Return all the roots of unity in this field, primitive or not.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<b> = NumberField(x^2 + 1)
            sage: zs = K.roots_of_unity(); zs
            [b, -1, -b, 1]
            sage: [z**K.number_of_roots_of_unity() for z in zs]
            [1, 1, 1, 1]
        """
        z = self.primitive_root_of_unity()
        n = self.zeta_order()
        return [z**k for k in range(1, n + 1)]

    def zeta_coefficients(self, n):
        """
        Compute the first `n` coefficients of the Dedekind zeta function of
        this field as a Dirichlet series.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: NumberField(x^2 + 1, 'a').zeta_coefficients(10)
            [1, 1, 0, 1, 2, 0, 0, 1, 1, 2]
        """
        return self.pari_nf().dirzetak(n)

    def solve_CRT(self, reslist, Ilist, check=True):
        r"""
        Solve a Chinese remainder problem over this number field.

        INPUT:

        - ``reslist`` -- list of residues, i.e. integral number field elements

        - ``Ilist`` -- list of integral ideals, assumed pairwise coprime

        - ``check`` -- boolean (default: ``True``); if ``True``, result is checked

        OUTPUT:

        An integral element `x` such that ``x - reslist[i]`` is in ``Ilist[i]`` for all `i`.

        .. NOTE::

           The current implementation requires the ideals to be pairwise
           coprime.  A more general version would be possible.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 - 10)
            sage: Ilist = [K.primes_above(p)[0] for p in prime_range(10)]
            sage: b = K.solve_CRT([1,2,3,4], Ilist, True)
            sage: all(b - i - 1 in Ilist[i] for i in range(4))
            True
            sage: Ilist = [K.ideal(a), K.ideal(2)]
            sage: K.solve_CRT([0,1], Ilist, True)
            Traceback (most recent call last):
            ...
            ArithmeticError: ideals in solve_CRT() must be pairwise coprime
            sage: Ilist[0] + Ilist[1]
            Fractional ideal (2, a)
        """
        n = len(reslist)
        try:
            reslist = [self(x) for x in reslist]
        except ValueError:
            raise ValueError("solve_CRT requires a list of arguments in the field")
        if n == 0:
            return self.zero()
        if n == 1:
            return reslist[0]
        if n == 2:
            try:
                r = Ilist[0].element_1_mod(Ilist[1])
            except TypeError:
                raise ArithmeticError("ideals in solve_CRT() must be pairwise coprime")
            x = ((1-r)*reslist[0]+r*reslist[1]).mod(prod(Ilist))
        else:  # n>2;, use induction / recursion
            x = self.solve_CRT([reslist[0], self.solve_CRT(reslist[1:], Ilist[1:])],
                               [Ilist[0], prod(Ilist[1:])], check=check)
        if check and not all(x - xi in Ii for xi, Ii in zip(reslist, Ilist)):
            raise RuntimeError("Error in number field solve_CRT()")
        return self(x)

    def valuation(self, prime):
        r"""
        Return the valuation on this field defined by ``prime``.

        INPUT:

        - ``prime`` -- a prime that does not split, a discrete
          (pseudo-)valuation or a fractional ideal

        EXAMPLES:

        The valuation can be specified with an integer ``prime`` that is
        completely ramified in ``R``::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.valuation(2)                                                        # needs sage.rings.padics
            2-adic valuation

        It can also be unramified in ``R``::

            sage: K.valuation(3)                                                        # needs sage.rings.padics
            3-adic valuation

        A ``prime`` that factors into pairwise distinct factors, results in an error::

            sage: K.valuation(5)                                                        # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: The valuation Gauss valuation induced by 5-adic valuation does not
            approximate a unique extension of 5-adic valuation with respect to x^2 + 1

        The valuation can also be selected by giving a valuation on the base
        ring that extends uniquely::

            sage: CyclotomicField(5).valuation(ZZ.valuation(5))                         # needs sage.rings.padics
            5-adic valuation

        When the extension is not unique, this does not work::

            sage: K.valuation(ZZ.valuation(5))                                          # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: The valuation Gauss valuation induced by 5-adic valuation does not
            approximate a unique extension of 5-adic valuation with respect to x^2 + 1

        For a number field which is of the form `K[x]/(G)`, you can specify a
        valuation by providing a discrete pseudo-valuation on `K[x]` which sends
        `G` to infinity. This lets us specify which extension of the 5-adic
        valuation we care about in the above example::

            sage: # needs sage.rings.padics
            sage: R.<x> = QQ[]
            sage: G5 = GaussValuation(R, QQ.valuation(5))
            sage: v = K.valuation(G5.augmentation(x + 2, infinity))
            sage: w = K.valuation(G5.augmentation(x + 1/2, infinity))
            sage: v == w
            False

        Note that you get the same valuation, even if you write down the
        pseudo-valuation differently::

            sage: # needs sage.rings.padics
            sage: ww = K.valuation(G5.augmentation(x + 3, infinity))
            sage: w is ww
            True

        The valuation ``prime`` does not need to send the defining polynomial `G`
        to infinity. It is sufficient if it singles out one of the valuations on
        the number field.  This is important if the prime only factors over the
        completion, i.e., if it is not possible to write down one of the factors
        within the number field::

            sage: # needs sage.rings.padics
            sage: v = G5.augmentation(x + 3, 1)
            sage: K.valuation(v)
            [ 5-adic valuation, v(x + 3) = 1 ]-adic valuation

        Finally, ``prime`` can also be a fractional ideal of a number field if it
        singles out an extension of a `p`-adic valuation of the base field::

            sage: K.valuation(K.fractional_ideal(a + 1))                                # needs sage.rings.padics
            2-adic valuation

        .. SEEALSO::

            :meth:`Order.valuation() <sage.rings.number_field.order.Order.valuation>`,
            :meth:`pAdicGeneric.valuation() <sage.rings.padics.padic_generic.pAdicGeneric.valuation>`
        """
        from sage.rings.padics.padic_valuation import pAdicValuation
        return pAdicValuation(self, prime)

    def some_elements(self):
        """
        Return a list of elements in the given number field.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: K.<a> =  QQ.extension(t^2 - 2); K
            Number Field in a with defining polynomial t^2 - 2
            sage: K.some_elements()
            [1, a, 2*a, 3*a - 4, 1/2, 1/3*a, 1/6*a, 0, 1/2*a, 2, ..., 12, -12*a + 18]

            sage: T.<u> = K[]
            sage: M.<b> = K.extension(t^3 - 5); M
            Number Field in b with defining polynomial t^3 - 5 over its base field
            sage: M.some_elements()
            [1, b, 1/2*a*b, ..., 2/5*b^2 + 2/5, 1/6*b^2 + 5/6*b + 13/6, 2]

        TESTS:

        This also works in trivial extensions::

            sage: R.<t> = QQ[]
            sage: K.<a> = QQ.extension(t); K
            Number Field in a with defining polynomial t
            sage: K.some_elements()
            [0, 1, 2, -1, 1/2, -1/2, 1/4, -2, 4]
        """
        elements = []

        polynomials = [self(f) for f in self.polynomial_ring().some_elements()]

        for numerator in polynomials:
            for denominator in polynomials:
                if denominator:
                    some_element = numerator/denominator
                    if some_element not in elements:
                        elements.append(some_element)

        return elements

    def lmfdb_page(self):
        r"""
        Open the LMFDB web page of the number field in a browser.

        See https://www.lmfdb.org

        EXAMPLES::

            sage: E = QuadraticField(-1)
            sage: E.lmfdb_page()  # optional -- webbrowser

        Even if the variable name is different it works::

            sage: R.<y>= PolynomialRing(QQ, "y")
            sage: K = NumberField(y^2 + 1 , "i")
            sage: K.lmfdb_page()  # optional -- webbrowser
        """
        import webbrowser
        from urllib.parse import quote
        lmfdb_url = 'https://www.lmfdb.org/NumberField/?jump={}'
        poly = self.absolute_polynomial()
        f = poly.parent().change_var('x')(poly)
        poly = pari(f).polredabs()
        url = lmfdb_url.format(quote(str(poly)))
        webbrowser.open(url)

    def maximal_order(self, v=None, assume_maximal='non-maximal-non-unique'):
        """
        Return the maximal order, i.e., the ring of integers, associated to
        this number field.

        INPUT:

        - ``v`` -- ``None``, a prime, or a list of integer primes (default: ``None``)

          - if ``None``, return the maximal order.

          - if a prime `p`, return an order that is `p`-maximal.

          - if a list, return an order that is maximal at each prime of these primes

        - ``assume_maximal`` -- ``True``, ``False``, ``None``, or
          ``'non-maximal-non-unique'`` (default: ``'non-maximal-non-unique'``)
          ignored when ``v`` is ``None``; otherwise, controls whether we assume
          that the order :meth:`order.is_maximal` outside of ``v``.

          - if ``True``, the order is assumed to be maximal at all primes.

          - if ``False``, the order is assumed to be non-maximal at some prime
            not in ``v``.

          - if ``None``, no assumptions are made about primes not in ``v``.

          - if ``'non-maximal-non-unique'`` (deprecated), like ``False``,
            however, the order is not a unique parent, so creating the same
            order later does typically not poison caches with the information
            that the order is not maximal.

        EXAMPLES:

        In this example, the maximal order cannot be generated by a single
        element::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + x^2 - 2*x+8)
            sage: o = k.maximal_order()
            sage: o
            Maximal Order generated by [1/2*a^2 + 1/2*a, a^2] in Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8

        We compute `p`-maximal orders for several `p`. Note
        that computing a `p`-maximal order is much faster in
        general than computing the maximal order::

            sage: p = next_prime(10^22)
            sage: q = next_prime(10^23)
            sage: K.<a> = NumberField(x^3 - p*q)

            sage: K.maximal_order([3], assume_maximal=None).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

            sage: K.maximal_order([2], assume_maximal=None).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

            sage: K.maximal_order([p], assume_maximal=None).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

            sage: K.maximal_order([q], assume_maximal=None).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

            sage: K.maximal_order([p, 3], assume_maximal=None).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

        An example with bigger discriminant::

            sage: p = next_prime(10^97)
            sage: q = next_prime(10^99)
            sage: K.<a> = NumberField(x^3 - p*q)
            sage: K.maximal_order(prime_range(10000), assume_maximal=None).basis()
            [1, a, a^2]

        An example in a relative number field::

            sage: K.<a, b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order()
            sage: OK.basis()
            [1, 1/2*a - 1/2*b, -1/2*b*a + 1/2, a]

            sage: charpoly(OK.1)
            x^2 + b*x + 1
            sage: charpoly(OK.2)
            x^2 - x + 1

            sage: O2 = K.order([3*a, 2*b])
            sage: O2.index_in(OK)
            144

        An order that is maximal at a prime. We happen to know that it is
        actually maximal and mark it as such::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.maximal_order(v=2, assume_maximal=True)
            Gaussian Integers generated by i in Number Field in i with defining polynomial x^2 + 1

        It is an error to create a maximal order and declare it non-maximal,
        however, such mistakes are only caught automatically if they evidently
        contradict previous results in this session::

            sage: K.maximal_order(v=2, assume_maximal=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot assume this order to be non-maximal
            because we already found it to be a maximal order

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(3*x^2 + 1)
            sage: K.maximal_order().basis()
            [3/2*a + 1/2, 3*a]

        The following was previously "ridiculously slow"; see :issue:`4738`::

            sage: K.<a, b> = NumberField([x^4 + 1, x^4 - 3])
            sage: K.maximal_order()
            Maximal Relative Order generated by
            [(131/4*b^3 + 9/4*b)*a^3 + (-69/4*b^2 - 521/4)*a^2 + (69/4*b^3 + 333/4*b)*a - 57/4*b^2 - 103/4,
             (21/2*b^3 - 39/4*b^2 - 25/4*b - 553/2)*a^3 + (69/4*b^3 - 21/4*b^2 + 473/2*b - 85/2)*a^2 + (37/4*b^3 - 189/2*b^2 + 55/2*b - 155/4)*a - 35*b^3 - 5*b^2 + 47/4*b - 83/4,
             (73/2*b^3 - 115/4*b^2 - 2*b - 1153/4)*a^3 + (56*b^3 - 59/4*b^2 + 237*b - 507/4)*a^2 + (17*b^3 - 353/4*b^2 + 73*b - 523/4)*a - 42*b^3 - 21/4*b^2 + 87/2*b - 117/4]
            in Number Field in a with defining polynomial x^4 + 1 over its base field

        An example with nontrivial ``v``::

            sage: L.<a, b> = NumberField([x^2 - 1000003, x^2 - 5*1000099^2])
            sage: O3 = L.maximal_order([3], assume_maximal=None)
            sage: O3.absolute_discriminant()
            400160824478095086350656915693814563600
            sage: O3.is_maximal()
            False
        """
        v = self._normalize_prime_list(v)

        if not v:
            if assume_maximal is False:
                raise ValueError("cannot assume that a maximal order is non-maximal")
            assume_maximal = True

        if assume_maximal == "non-maximal-non-unique":
            deprecation(33386, 'maximal_order(v=[primes], assume_maximal="non-maximal-non-unique") has been deprecated since it can silently produce wrong results and does not play nicely with caching. An order that is maximal at some primes should be created with assume_maximal=None instead to make no assumptions about maximality at other primes.')

        return self._maximal_order(v, assume_maximal=assume_maximal)


class NumberField_absolute(NumberField_generic):
    def __init__(self, polynomial, name, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None, structure=None):
        r"""
        Function to initialize an absolute number field.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K = NumberField(x^17 + 3, 'a'); K
            Number Field in a with defining polynomial x^17 + 3
            sage: type(K)
            <class 'sage.rings.number_field.number_field.NumberField_absolute_with_category'>
            sage: TestSuite(K).run()
        """
        NumberField_generic.__init__(self, polynomial, name, latex_name, check, embedding,
                                     assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structure=structure)
        self._element_class = number_field_element.NumberFieldElement_absolute
        self._zero_element = self._element_class(self, 0)
        self._one_element = self._element_class(self, 1)

        self._init_embedding_approx()

    def _coerce_from_other_number_field(self, x):
        r"""
        Coerce a number field element x into this number field.

        Unless `x` is in ``QQ``, this requires ``x.parent()`` and
        ``self`` to have compatible embeddings: either they both embed
        in a common field, or there is an embedding of ``x.parent()``
        into ``self`` or the other way around.  If no compatible
        embeddings are found and `x` is not in ``QQ``, then raise
        :exc:`TypeError`.  This guarantees that these conversions respect
        the field operations and conversions between several fields
        commute.

        REMARK:

        The name of this method was chosen for historical reasons.
        In fact, what it does is not a coercion but a conversion.

        INPUT:

        - ``x`` -- an element of some number field

        OUTPUT: an element of ``self`` corresponding to ``x``

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 + 2)
            sage: L.<b> = NumberField(x^2 + 1)
            sage: K._coerce_from_other_number_field(L(2/3))
            2/3
            sage: L._coerce_from_other_number_field(K(0))
            0
            sage: K._coerce_from_other_number_field(b)
            Traceback (most recent call last):
            ...
            TypeError: No compatible natural embeddings found for
            Number Field in a with defining polynomial x^3 + 2 and
            Number Field in b with defining polynomial x^2 + 1

        Two number fields both containing `i`::

            sage: K.<a> = NumberField(x^4 + 6*x^2 + 1, embedding=CC(-2.4*I))
            sage: L.<b> = NumberField(x^4 + 8*x^2 + 4, embedding=CC(2.7*I))
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
            sage: L((Ki+2)^1000)
            737533628...075020804*b^3 + 442520177...450124824*b + 793311113...453515313

        This fails if we don't specify the embeddings::

            sage: K.<a> = NumberField(x^4 + 6*x^2 + 1)
            sage: L.<b> = NumberField(x^4 + 8*x^2 + 4)
            sage: L(1/2*a^3 + 5/2*a)
            Traceback (most recent call last):
            ...
            TypeError: No compatible natural embeddings found for
            Number Field in b with defining polynomial x^4 + 8*x^2 + 4 and
            Number Field in a with defining polynomial x^4 + 6*x^2 + 1

        Embeddings can also be `p`-adic::

            sage: # needs sage.rings.padics
            sage: F = Qp(73)
            sage: K.<a> = NumberField(x^4 + 6*x^2 + 1,
            ....:                     embedding=F(1290990671961076190983179596556712119))
            sage: L.<b> = NumberField(x^4 + 8*x^2 + 4,
            ....:                     embedding=F(1773398470280167815153042237103591466))
            sage: L(2*a^3 + 10*a + 3)
            b^3 + 6*b + 3

        If we take the same non-Galois number field with two different
        embeddings, conversion fails::

            sage: K.<a> = NumberField(x^3 - 4*x + 1, embedding=0.254)
            sage: L.<b> = NumberField(x^3 - 4*x + 1, embedding=1.86)
            sage: L(a)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert a to Number Field in b with defining polynomial x^3 - 4*x + 1
            with b = 1.860805853111704? (using the specified embeddings)

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
            sage: i = (x^2 + 1).roots(ring=K)[0][0]
            sage: r2 = (x^2 - 2).roots(ring=K)[0][0]
            sage: r3 = (x^2 - 3).roots(ring=K)[0][0]
            sage: L2.<a2>, phi2 = K.subfield(r2 + i)
            sage: L3.<a3>, phi3 = K.subfield(r3 + i)
            sage: i_in_L2 = L2(i); i_in_L2
            1/6*a2^3 + 1/6*a2
            sage: i_in_L3 = L3(i); i_in_L3
            1/8*a3^3
            sage: L2(i_in_L3) == i_in_L2
            True
            sage: L3(i_in_L2) == i_in_L3
            True

        TESTS:

        The following was fixed in :issue:`8800`::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3 - 5, embedding=0)
            sage: L.<b> = K.extension(x^2 + a)
            sage: F,R = L.construction()
            sage: F(R) == L   #indirect doctest
            True

        AUTHORS:

        - Jeroen Demeyer (2011-09-30): :issue:`11869`
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
        from sage.rings.complex_mpfr import ComplexField
        from sage.rings.real_mpfr import RealField
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
            raise TypeError("No compatible natural embeddings found for %s and %s" % (KF, LF))

        # List of candidates for K(x)
        f = x.minpoly()
        ys = f.roots(ring=K, multiplicities=False)
        if not ys:
            raise ValueError("cannot convert %s to %s (regardless of embeddings)" % (x, K))

        # Define a function are_roots_equal to determine whether two
        # roots of f are equal.  A simple a == b does not suffice for
        # inexact fields because of floating-point errors.
        if F.is_exact():
            are_roots_equal = lambda a, b: a == b
        else:
            #  Compute a lower bound on the distance between the roots of f.
            #  This essentially gives the precision to work with.

            # A function
            # log2abs: F --> RR
            #          x |-> log2(abs(x))
            # This should work for all fields F with an absolute value.
            # The p-adic absolute value goes into QQ, so we need the RR().
            log2abs = lambda x: RR(F(x).abs()).log2()

            # Compute half Fujiwara's bound on the roots of f
            n = f.degree()
            log_half_root_bound = log2abs(f[0]/2)/n
            for i in range(1, n):
                bd = log2abs(f[i])/(n-i)
                log_half_root_bound = max(bd, log_half_root_bound)
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
            are_roots_equal = lambda a, b: (a-b).abs() <= eps
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
        raise ValueError("cannot convert %s to %s (using the specified embeddings)" % (x, K))

    def _coerce_map_from_(self, R):
        r"""
        Canonical coercion of a ring `R` into ``self``.

        Currently any ring coercing into the base ring canonically coerces
        into this field, as well as orders in any number field coercing into
        this field, and of course the field itself as well.

        Two embedded number fields may mutually coerce into each other, if
        the pushout of the two ambient fields exists and if it is possible
        to construct an :class:`~sage.rings.number_field.number_field_morphisms.EmbeddedNumberFieldMorphism`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: S.<y> = NumberField(x^3 + x + 1)
            sage: S.coerce(int(4))  # indirect doctest
            4
            sage: S.coerce(-Integer(2))
            -2
            sage: z = S.coerce(-7/8); z, type(z)
            (-7/8, <class 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>)
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
              To:   Number Field in a with defining polynomial x^2 + 3 with a = 1.732050807568878?*I
              Defn: zeta3 -> 1/2*a - 1/2

        Two embedded number fields with mutual coercions (testing against a
        bug that was fixed in :issue:`8800`)::

            sage: K.<r4> = NumberField(x^4 - 2)
            sage: L1.<r2_1> = NumberField(x^2 - 2, embedding=r4**2)
            sage: L2.<r2_2> = NumberField(x^2 - 2, embedding=-r4**2)
            sage: r2_1 + r2_2    # indirect doctest
            0
            sage: (r2_1 + r2_2).parent() is L1
            True
            sage: (r2_2 + r2_1).parent() is L2
            True

        Coercion of an order (testing against a bug that was fixed in
        :issue:`8800`)::

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
            sage: K.coerce(sqrt(2))                                                     # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Symbolic Ring to Number Field in a
            with defining polynomial x^2 - 2 with a = 1.414213562373095?

        TESTS::

            sage: K.<a> = NumberField(polygen(QQ)^3 - 2)
            sage: type(K.coerce_map_from(QQ))
            <class 'sage.structure.coerce_maps.DefaultConvertMap_unique'>

        Make sure we still get our optimized morphisms for special fields::

            sage: K.<a> = NumberField(polygen(QQ)^2 - 2)
            sage: type(K.coerce_map_from(QQ))
            <class 'sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element'>
        """
        if R is int:
            return self._generic_coerce_map(R)
        if R in (ZZ, QQ, self.base()):
            return self._generic_coerce_map(R)
        if isinstance(R, sage.rings.abc.Order) and self.has_coerce_map_from(R.number_field()):
            return self._generic_coerce_map(R)
        # R is not QQ by the above tests
        if isinstance(R, number_field_base.NumberField) and R.coerce_embedding() is not None:
            if self.coerce_embedding() is not None:
                try:
                    return number_field_morphisms.EmbeddedNumberFieldMorphism(R, self)
                except ValueError:  # no common embedding found
                    return None
            else:
                # R is embedded, self isn't. So, we could only have
                # the forgetful coercion. But this yields to non-commuting
                # coercions, as was pointed out at issue #8800
                return None

    def base_field(self):
        r"""
        Return the base field of ``self``, which is always ``QQ``.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.base_field()
            Rational Field
        """
        return QQ

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: it = iter(K)
            sage: [next(it) for _ in range(20)]
            [0, 1, zeta5, zeta5^2, zeta5^3, -1, zeta5 + 1, zeta5^2 + 1,
             zeta5^3 + 1, -zeta5, zeta5^2 + zeta5, zeta5^3 + zeta5, -zeta5^2,
             zeta5^3 + zeta5^2, -zeta5^3, 1/2, zeta5 - 1, zeta5^2 - 1,
             zeta5^3 - 1, -zeta5 + 1]
        """
        from sage.categories.sets_cat import cartesian_product
        M, f, g = self.free_module()
        # We iterate over the Cartesian product since the free module does
        #   not iterate using the diagonal embedding.
        for v in cartesian_product([M.base_ring()] * M.dimension()):
            yield f(M(list(v)))

    def is_absolute(self):
        r"""
        Return ``True`` since ``self`` is an absolute field.

        EXAMPLES::

            sage: K = CyclotomicField(5)
            sage: K.is_absolute()
            True
        """
        return True

    def absolute_polynomial(self):
        r"""
        Return absolute polynomial that defines this absolute field. This
        is the same as :meth:`polynomial`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 - 17)
            sage: K.absolute_generator()
            a
        """
        return self.gen()

    def optimized_representation(self, name=None, both_maps=True):
        r"""
        Return a field isomorphic to ``self`` with a better defining polynomial
        if possible, along with field isomorphisms from the new field to
        ``self`` and from ``self`` to the new field.

        EXAMPLES: We construct a compositum of 3 quadratic fields, then
        find an optimized representation and transform elements back and
        forth.

        ::

            sage: x = polygen(QQ, 'x')
            sage: K = NumberField([x^2 + p for p in [5, 3, 2]],'a').absolute_field('b'); K
            Number Field in b with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
            sage: L, from_L, to_L = K.optimized_representation()
            sage: L    # your answer may different, since algorithm is random
            Number Field in b1 with defining polynomial x^8 + 4*x^6 + 7*x^4 + 36*x^2 + 81
            sage: to_L(K.0)   # random
            4/189*b1^7 + 1/63*b1^6 + 1/27*b1^5 - 2/9*b1^4 - 5/27*b1^3 - 8/9*b1^2 + 3/7*b1 - 3/7
            sage: from_L(L.0)   # random
            1/1152*b^7 - 1/192*b^6 + 23/576*b^5 - 17/96*b^4 + 37/72*b^3 - 5/6*b^2 + 55/24*b - 3/4

        The transformation maps are mutually inverse isomorphisms.

        ::

            sage: from_L(to_L(K.0)) == K.0
            True
            sage: to_L(from_L(L.0)) == L.0
            True

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(7/9*x^3 + 7/3*x^2 - 56*x + 123)
            sage: K.optimized_representation()  # representation varies, not tested
            (Number Field in a1 with defining polynomial x^3 - 7*x - 7,
             Ring morphism:
               From: Number Field in a1 with defining polynomial x^3 - 7*x - 7
               To:   Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               Defn: a1 |--> 7/225*a^2 - 7/75*a - 42/25,
             Ring morphism:
               From: Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               To:   Number Field in a1 with defining polynomial x^3 - 7*x - 7
               Defn: a |--> -15/7*a1^2 + 9)

        TESTS:

        We test the above doctest. The representation depends on the PARI version::

            sage: K.<a> = NumberField(7/9*x^3 + 7/3*x^2 - 56*x + 123)
            sage: N, M1, M2 = K.optimized_representation(); N, M1, M2
            (Number Field in a1 with defining polynomial x^3 - 7*x - 7,
             Ring morphism:
               From: Number Field in a1 with defining polynomial x^3 - 7*x - 7
               To:   Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               Defn: a1 |--> ...,
             Ring morphism:
               From: Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               To:   Number Field in a1 with defining polynomial x^3 - 7*x - 7
               Defn: a |--> ...)
            sage: a1 = M1.domain().gens()[0]
            sage: M2(a) in (-15/7*a1^2 + 9, -60/7*a1^2 + 15*a1 + 39)
            True
            sage: M1(M2(a)) == a
            True
        """
        if name is None:
            name = self.variable_names()
        name = normalize_names(1, name)[0]

        f = self.absolute_polynomial().__pari__()

        g, alpha = f.polredbest(flag=1)
        beta = alpha.modreverse()

        b = self(QQ['x'](lift(beta)))
        h = QQ['x'](g)

        embedding = None
        if self.coerce_embedding() is not None:
            embedding = self.coerce_embedding()(b)
        # trac 7695 add a _ to prevent zeta70 etc.
        if name[-1].isdigit():
            new_name = name + '_1'
        else:
            new_name = name + '1'

        K = NumberField(h, names=new_name, embedding=embedding)
        from_K = K.hom([b])

        if both_maps:
            a = K(alpha)
            to_K = self.hom([a])

            return K, from_K, to_K

        return K, from_K

    def optimized_subfields(self, degree=0, name=None, both_maps=True):
        """
        Return optimized representations of many (but *not* necessarily
        all!) subfields of ``self`` of the given ``degree``, or of all possible degrees if
        ``degree`` is 0.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K = NumberField([x^2 + p for p in [5, 3, 2]],'a').absolute_field('b'); K
            Number Field in b with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
            sage: L = K.optimized_subfields(name='b')
            sage: L[0][0]
            Number Field in b0 with defining polynomial x
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
              To:   Number Field in a1 with defining polynomial
                    x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
              Defn: b2 |--> -5/1152*a1^7 + 1/96*a1^6 - 97/576*a1^5 + 17/48*a1^4
                            - 95/72*a1^3 + 17/12*a1^2 - 53/24*a1 - 1

        The ``to_M`` map is ``None``, since there is no map from `K` to `M`::

            sage: to_M

        We apply the from_M map to the generator of M, which gives a
        rather large element of `K`::

            sage: from_M(M.0)          # random
            -5/1152*a1^7 + 1/96*a1^6 - 97/576*a1^5 + 17/48*a1^4
             - 95/72*a1^3 + 17/12*a1^2 - 53/24*a1 - 1

        Nevertheless, that large-ish element lies in a degree 4 subfield::

            sage: from_M(M.0).minpoly()   # random
            x^4 - 5*x^2 + 25

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(2*x^4 + 6*x^2 + 1/2)
            sage: K.optimized_subfields()
            [(Number Field in a0 with defining polynomial x,
              Ring morphism:
                From: Number Field in a0 with defining polynomial x
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: 0 |--> 0,
              None),
             (Number Field in a1 with defining polynomial x^2 - 2*x + 2,
              Ring morphism:
                From: Number Field in a1 with defining polynomial x^2 - 2*x + 2
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: a1 |--> a^3 + 7/2*a + 1,
              None),
             (Number Field in a2 with defining polynomial x^2 - 2*x + 2,
              Ring morphism:
                From: Number Field in a2 with defining polynomial x^2 - 2*x + 2
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: a2 |--> -a^3 - 7/2*a + 1,
              None),
             (Number Field in a3 with defining polynomial x^2 - 2,
              Ring morphism:
                From: Number Field in a3 with defining polynomial x^2 - 2
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: a3 |--> a^2 + 3/2,
              None),
             (Number Field in a4 with defining polynomial x^2 + 1,
              Ring morphism:
                From: Number Field in a4 with defining polynomial x^2 + 1
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: a4 |--> a^3 + 7/2*a,
              None),
             (Number Field in a5 with defining polynomial x^2 + 2,
              Ring morphism:
                From: Number Field in a5 with defining polynomial x^2 + 2
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: a5 |--> 2*a^3 + 5*a,
              None),
             (Number Field in a6 with defining polynomial x^4 + 1,
              Ring morphism:
                From: Number Field in a6 with defining polynomial x^4 + 1
                To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                Defn: a6 |--> a^3 + 1/2*a^2 + 5/2*a + 3/4,
              Ring morphism:
                From: Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
                To:   Number Field in a6 with defining polynomial x^4 + 1
                Defn: a |--> -1/2*a6^3 + a6^2 - 1/2*a6)]
        """
        return self._subfields_helper(degree=degree, name=name,
                                      both_maps=both_maps, optimize=True)

    def change_names(self, names):
        r"""
        Return number field isomorphic to ``self`` but with the given generator
        name.

        INPUT:

        - ``names`` -- should be exactly one variable name

        Also, ``K.structure()`` returns ``from_K`` and ``to_K``,
        where ``from_K`` is an isomorphism from `K` to ``self`` and ``to_K`` is an
        isomorphism from ``self`` to `K`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        Return all subfields of ``self`` of the given ``degree``,
        or of all possible degrees if ``degree`` is 0.  The subfields are returned as
        absolute fields together with an embedding into ``self``.  For the case of the
        field itself, the reverse isomorphism is also provided.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField([x^3 - 2, x^2 + x + 1])
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

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: K.<a> = NumberField(2*x^4 + 6*x^2 + 1/2)
            sage: K.subfields()
            [(Number Field in a0 with defining polynomial x, Ring morphism:
              From: Number Field in a0 with defining polynomial x
              To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
              Defn: 0 |--> 0, None),
            (Number Field in a1 with defining polynomial x^2 - 2, Ring morphism:
              From: Number Field in a1 with defining polynomial x^2 - 2
              To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
              Defn: a1 |--> -a^2 - 3/2, None),
            (Number Field in a2 with defining polynomial x^2 + 4, Ring morphism:
              From: Number Field in a2 with defining polynomial x^2 + 4
              To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
              Defn: a2 |--> 2*a^3 + 7*a, None),
            (Number Field in a3 with defining polynomial x^2 + 2, Ring morphism:
              From: Number Field in a3 with defining polynomial x^2 + 2
              To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
              Defn: a3 |--> -2*a^3 - 5*a, None),
            (Number Field in a4 with defining polynomial x^4 + 1, Ring morphism:
              From: Number Field in a4 with defining polynomial x^4 + 1
              To:   Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
              Defn: a4 |--> -a^3 - 1/2*a^2 - 5/2*a - 3/4, Ring morphism:
              From: Number Field in a with defining polynomial 2*x^4 + 6*x^2 + 1/2
              To:   Number Field in a4 with defining polynomial x^4 + 1
              Defn: a |--> 1/2*a4^3 + a4^2 + 1/2*a4)]
        """
        return self._subfields_helper(degree=degree, name=name,
                                      both_maps=True, optimize=False)

    def _subfields_helper(self, degree=0, name=None, both_maps=True, optimize=False):
        """
        Internal function: common code for :meth:`optimized_subfields` and :meth:`subfields`.

        TESTS:

        Let's make sure embeddings are being respected::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^4 - 23, embedding=50)
            sage: K, CDF(a)
            (Number Field in a with defining polynomial x^4 - 23 with a = 2.189938703094843?,
             2.1899387030948425)
            sage: Ss = K.subfields(); len(Ss)  # indirect doctest
            3
            sage: diffs = [ S.coerce_embedding()(S.gen()) - CDF(S_into_K(S.gen())) for S, S_into_K, _ in Ss ]
            sage: all(abs(diff) < 1e-5 for diff in diffs)
            True

            sage: L1, _, _ = K.subfields(2)[0]; L1, CDF(L1.gen())  # indirect doctest
            (Number Field in a0 with defining polynomial x^2 - 23 with a0 = -4.795831523312720?,
             -4.795831523312719)

        If we take a different embedding of the large field, we get a
        different embedding of the degree 2 subfield::

            sage: K.<a> = NumberField(x^4 - 23, embedding=-50)
            sage: L2, _, _ = K.subfields(2)[0]; L2, CDF(L2.gen())  # indirect doctest
            (Number Field in a0 with defining polynomial x^2 - 23 with a0 = -4.795831523312720?,
             -4.795831523312719)

        Test for :issue:`7695`::

            sage: F = CyclotomicField(7)
            sage: K = F.subfields(3)[0][0]
            sage: K
            Number Field in zeta7_0 with defining polynomial x^3 + x^2 - 2*x - 1 with zeta7_0 = 1.246979603717467?
        """
        if name is None:
            name = self.variable_names()
        name = normalize_names(1, name)[0]
        try:
            return self.__subfields[name, degree, both_maps, optimize]
        except AttributeError:
            self.__subfields = {}
        except KeyError:
            pass
        f = self.pari_polynomial()
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
            a = self(elts[i], check=False)
            if self.coerce_embedding() is not None:
                embedding = self.coerce_embedding()(a)
            # trac 7695 add a _ to prevent zeta70 etc.
            if name[-1].isdigit():
                new_name = name + '_' + str(i)
            else:
                new_name = name + str(i)
            K = NumberField(f, names=new_name, embedding=embedding)

            from_K = K.hom([a])  # check=False here ??   would be safe unless there are bugs.

            if both_maps and K.degree() == self.degree():
                g = K['x'](self.polynomial())
                a = from_K(K.gen())
                for root in g.roots(multiplicities=False):
                    to_K = self.hom([root])    # check=False here ??
                    if to_K(a) == K.gen():
                        break
            else:
                to_K = None
            ans.append((K, from_K, to_K))
        ans = Sequence(ans, immutable=True, cr=ans != [])
        self.__subfields[name, degree, both_maps, optimize] = ans
        return ans

    @cached_method
    def _maximal_order(self, v=(), assume_maximal=None):
        r"""
        Return a :meth:`maximal_order` in this number field.

        TESTS::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + x^2 - 2*x+8)
            sage: k.maximal_order() is k.maximal_order()  # indirect doctest
            True
        """
        B = [self(b, check=False) for b in self._pari_integral_basis(v=v)]

        from sage.rings.number_field.order import absolute_order_from_module_generators
        return absolute_order_from_module_generators(B, check_integral=False, check_rank=False,
                                                     check_is_ring=False, is_maximal=assume_maximal,
                                                     is_maximal_at=v)

    def order(self, *args, **kwds):
        r"""
        Return the order with given ring generators in the maximal order of
        this number field.

        INPUT:

        - ``gens`` -- list of elements in this number field; if no generators
          are given, just returns the cardinality of this number field
          (`\infty`) for consistency.

        - ``check_is_integral`` -- boolean (default: ``True``); whether to
          check that each generator is integral

        - ``check_rank`` -- boolean (default: ``True``); whether to check that
          the ring generated by ``gens`` is of full rank

        - ``allow_subfield`` -- boolean (default: ``False``); if ``True`` and
          the generators do not generate an order, i.e., they generate a
          subring of smaller rank, instead of raising an error, return an order
          in a smaller number field

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<i> = NumberField(x^2 + 1)
            sage: k.order(2*i)
            Order of conductor 2 generated by 2*i in Number Field in i with defining polynomial x^2 + 1
            sage: k.order(10*i)
            Order of conductor 10 generated by 10*i in Number Field in i with defining polynomial x^2 + 1
            sage: k.order(3)
            Traceback (most recent call last):
            ...
            ValueError: the rank of the span of gens is wrong
            sage: k.order(i/2)
            Traceback (most recent call last):
            ...
            ValueError: each generator must be integral

        Alternatively, an order can be constructed by adjoining elements to
        `\ZZ`::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: ZZ[a]
            Order generated by a0 in Number Field in a0 with defining polynomial x^3 - 2 with a0 = a

        TESTS:

        We verify that :issue:`2480` is fixed::

            sage: K.<a> = NumberField(x^4 + 4*x^2 + 2)
            sage: B = K.integral_basis()
            sage: K.order(*B)
            Maximal Order generated by a in Number Field in a with defining polynomial x^4 + 4*x^2 + 2
            sage: K.order(B)
            Maximal Order generated by a in Number Field in a with defining polynomial x^4 + 4*x^2 + 2
            sage: K.order(gens=B)
            Maximal Order generated by a in Number Field in a with defining polynomial x^4 + 4*x^2 + 2
        """
        # set gens appropriately from the arguments
        gens = kwds.pop('gens', args)

        if len(gens) == 0:
            return NumberField_generic.order(self)
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens = map(self, gens)
        return self._order(tuple(gens), **kwds)

    @cached_method
    def _order(self, gens, **kwds):
        r"""
        Helper method for :meth:`order` which adds caching. See :meth:`order`
        for a description of the parameters and keyword parameters.

        TESTS:

        Test that caching works::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: K.order(a) is K.order(a)      # indirect doctest
            True

        Keywords have no influence on the caching::

            sage: K.order(a) is K.order(a,check_is_integral=True) is K.order(a,check_is_integral=False)
            True

        Even if the order lives in a different field, caching works (currently,
        however, ``allow_subfield`` is incorrect :issue:`16046`)::

            sage: K.<a> = NumberField(x**4 + 3)
            sage: o = K.order([a**2], allow_subfield=True)
            sage: o is K.order([a**2], allow_subfield=True)
            True

        Different generators for the same order::

            sage: K.order(a) is K.order(a,a^2) is K.order(a^2,a)
            True
        """
        from sage.rings.number_field.order import absolute_order_from_ring_generators
        return absolute_order_from_ring_generators(gens, **kwds)

    @cached_method(key=lambda self, base, basis, map: (base or self.base_ring(), basis, map))
    def free_module(self, base=None, basis=None, map=True):
        r"""
        Return a vector space `V` and isomorphisms ``self`` `\to` `V` and `V` `\to` ``self``.

        INPUT:

        - ``base`` -- a subfield (default: ``None``); the returned vector
          space is over this subfield `R`, which defaults to the base field of this
          function field

        - ``basis`` -- a basis for this field over the base

        - ``maps`` -- boolean (default: ``True``); whether to return
          `R`-linear maps to and from `V`

        OUTPUT:

        - ``V`` -- a vector space over the rational numbers

        - ``from_V`` -- an isomorphism from `V` to ``self`` (if requested)

        - ``to_V`` -- an isomorphism from ``self`` to `V` (if requested)

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 2)
            sage: V, from_V, to_V  = k.free_module()
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
        if base is None:
            base = QQ
        elif base is self:
            return super().free_module(base=base, basis=basis, map=map)
        if basis is not None or base is not QQ:
            raise NotImplementedError
        V = QQ**self.degree()
        if not map:
            return V
        from_V = maps.MapVectorSpaceToNumberField(V, self)
        to_V = maps.MapNumberFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def absolute_vector_space(self, *args, **kwds):
        r"""
        Return vector space over `\QQ` corresponding to this
        number field, along with maps from that space to this number field
        and in the other direction.

        For an absolute extension this is identical to
        :meth:`vector_space`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        return self.free_module(*args, **kwds)

    def _galois_closure_and_embedding(self, names=None):
        r"""
        Return number field `K` that is the Galois closure of ``self`` and an
        embedding of ``self`` into `K`.

        INPUT:

        - ``names`` -- variable name for Galois closure

        .. warning::

           This is an internal function; see :meth:`galois_closure`.

        EXAMPLES:

        For medium-sized Galois groups of fields with small discriminants,
        this computation is feasible::

            sage: # needs sage.groups
            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^6 + 4*x^2 + 2)
            sage: K.galois_group().order()
            48
            sage: L, phi = K._galois_closure_and_embedding('c')
            sage: phi.domain() is K, phi.codomain() is L
            (True, True)
            sage: L
            Number Field in c with defining polynomial x^48 + 8*x^46 - 20*x^44 - 520*x^42 + 12106*x^40 - 68344*x^38 + 463156*x^36 - 1823272*x^34 + 8984591*x^32 - 25016080*x^30 + 84949344*x^28 - 163504384*x^26 + 417511068*x^24 - 394687376*x^22 + 836352224*x^20 + 72845696*x^18 + 1884703919*x^16 + 732720520*x^14 + 3960878676*x^12 + 2507357768*x^10 + 5438373834*x^8 + 3888508744*x^6 + 4581432268*x^4 + 1765511400*x^2 + 1723993441
            sage: K.defining_polynomial()( phi(K.gen()) )
            0
        """
        if names is None:
            raise TypeError("You must specify the name of the generator.")

        try:
            # compose with variable renaming
            L = self.__galois_closure.change_names(names)
            L_to_orig, orig_to_L = L.structure()
            # "flatten" the composition by hand
            self_into_L = self.hom([(orig_to_L * self.__galois_closure_embedding)(self.gen())])
            return (L, self_into_L)
        except AttributeError:
            pass

        # Compute degree of Galois closure if possible
        deg = self.galois_group().easy_order()

        L, self_into_L = self.defining_polynomial().change_ring(self).splitting_field(names, map=True, degree_multiple=deg)
        self.__galois_closure = L
        self.__galois_closure_embedding = self_into_L
        return (self.__galois_closure, self.__galois_closure_embedding)

    def galois_closure(self, names=None, map=False):
        """
        Return number field `K` that is the Galois closure of ``self``,
        i.e., is generated by all roots of the defining polynomial of
        ``self``, and possibly an embedding of ``self`` into `K`.

        INPUT:

        - ``names`` -- variable name for Galois closure

        - ``map`` -- boolean (default: ``False``); also return an embedding of
          ``self`` into `K`

        EXAMPLES::

            sage: # needs sage.groups
            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^4 - 2)
            sage: M = K.galois_closure('b'); M
            Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
            sage: L.<a2> = K.galois_closure(); L
            Number Field in a2 with defining polynomial x^8 + 28*x^4 + 2500
            sage: K.galois_group(names=("a3")).order()
            8

        ::

            sage: # needs sage.groups
            sage: phi = K.embeddings(L)[0]
            sage: phi(K.0)
            1/120*a2^5 + 19/60*a2
            sage: phi(K.0).minpoly()
            x^4 - 2

            sage: # needs sage.groups
            sage: L, phi = K.galois_closure('b', map=True)
            sage: L
            Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
            sage: phi
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 - 2
              To:   Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
              Defn: a |--> 1/240*b^5 - 41/120*b

        A cyclotomic field is already Galois::

            sage: # needs sage.groups
            sage: K.<a> = NumberField(cyclotomic_polynomial(23))
            sage: L.<z> = K.galois_closure()
            sage: L
            Number Field in z with defining polynomial
             x^22 + x^21 + x^20 + x^19 + x^18 + x^17 + x^16 + x^15 + x^14 + x^13 + x^12
              + x^11 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1

        TESTS:

        Let's make sure we're renaming correctly::

            sage: # needs sage.groups
            sage: K.<a> = NumberField(x^4 - 2)
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
        return (L, self_into_L) if map else L

    def automorphisms(self):
        r"""
        Compute all Galois automorphisms of ``self``.

        This uses PARI's :pari:`nfgaloisconj` and is much faster than
        root finding for many fields.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 10000)
            sage: K.automorphisms()
            [Ring endomorphism of Number Field in a with defining polynomial x^2 + 10000
               Defn: a |--> a,
             Ring endomorphism of Number Field in a with defining polynomial x^2 + 10000
               Defn: a |--> -a]

        Here's a larger example, that would take some time if we found
        roots instead of using PARI's specialized machinery::

            sage: K = NumberField(x^6 - x^4 - 2*x^2 + 1, 'a')
            sage: len(K.automorphisms())
            2

        `L` is the Galois closure of `K`::

            sage: L = NumberField(x^24 - 84*x^22 + 2814*x^20 - 15880*x^18 - 409563*x^16
            ....:                  - 8543892*x^14 + 25518202*x^12 + 32831026956*x^10
            ....:                  - 672691027218*x^8 - 4985379093428*x^6 + 320854419319140*x^4
            ....:                  + 817662865724712*x^2 + 513191437605441, 'a')
            sage: len(L.automorphisms())
            24

        Number fields defined by non-monic and non-integral
        polynomials are supported (:issue:`252`)::

            sage: R.<x> = QQ[]
            sage: f = 7/9*x^3 + 7/3*x^2 - 56*x + 123
            sage: K.<a> = NumberField(f)
            sage: A = K.automorphisms(); A
            [Ring endomorphism of Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               Defn: a |--> a,
             Ring endomorphism of Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               Defn: a |--> -7/15*a^2 - 18/5*a + 96/5,
             Ring endomorphism of Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
               Defn: a |--> 7/15*a^2 + 13/5*a - 111/5]
            sage: prod(x - sigma(a) for sigma in A) == f.monic()
            True
        """
        return self.embeddings(self)

    @cached_method
    def embeddings(self, K):
        """
        Compute all field embeddings of this field into the field `K` (which need
        not even be a number field, e.g., it could be the complex numbers).
        This will return an identical result when given `K` as input again.

        If possible, the most natural embedding of this field into `K`
        is put first in the list.

        INPUT:

        - ``K`` -- a field

        EXAMPLES::

            sage: # needs sage.groups
            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<a1> = K.galois_closure(); L
            Number Field in a1 with defining polynomial x^6 + 108
            sage: K.embeddings(L)[0]
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in a1 with defining polynomial x^6 + 108
              Defn: a |--> 1/18*a1^4
            sage: K.embeddings(L) is K.embeddings(L)
            True

        We embed a quadratic field into a cyclotomic field::

            sage: L.<a> = QuadraticField(-7)
            sage: K = CyclotomicField(7)
            sage: L.embeddings(K)
            [Ring morphism:
               From: Number Field in a with defining polynomial x^2 + 7 with a = 2.645751311064591?*I
               To:   Cyclotomic Field of order 7 and degree 6
               Defn: a |--> 2*zeta7^4 + 2*zeta7^2 + 2*zeta7 + 1,
             Ring morphism:
               From: Number Field in a with defining polynomial x^2 + 7 with a = 2.645751311064591?*I
               To:   Cyclotomic Field of order 7 and degree 6
               Defn: a |--> -2*zeta7^4 - 2*zeta7^2 - 2*zeta7 - 1]

        We embed a cubic field in the complex numbers::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: K.embeddings(CC)  # abs tol 1e-12
            [Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Complex Field with 53 bits of precision
               Defn: a |--> -0.629960524947437 - 1.09112363597172*I,
             Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Complex Field with 53 bits of precision
               Defn: a |--> -0.629960524947437 + 1.09112363597172*I,
             Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Complex Field with 53 bits of precision
               Defn: a |--> 1.25992104989487]

        Some more (possible and impossible) embeddings of cyclotomic fields::

            sage: CyclotomicField(5).embeddings(QQbar)
            [Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Algebraic Field
               Defn: zeta5 |--> 0.3090169943749474? + 0.9510565162951536?*I,
             Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Algebraic Field
               Defn: zeta5 |--> -0.8090169943749474? + 0.5877852522924731?*I,
             Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Algebraic Field
               Defn: zeta5 |--> -0.8090169943749474? - 0.5877852522924731?*I,
             Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Algebraic Field
               Defn: zeta5 |--> 0.3090169943749474? - 0.9510565162951536?*I]
            sage: CyclotomicField(3).embeddings(CyclotomicField(7))
            []
            sage: CyclotomicField(3).embeddings(CyclotomicField(6))
            [Ring morphism:
               From: Cyclotomic Field of order 3 and degree 2
               To:   Cyclotomic Field of order 6 and degree 2
               Defn: zeta3 |--> zeta6 - 1,
             Ring morphism:
               From: Cyclotomic Field of order 3 and degree 2
               To:   Cyclotomic Field of order 6 and degree 2
               Defn: zeta3 |--> -zeta6]

        Test that :issue:`15053` is fixed::

            sage: K = NumberField(x^3 - 2, 'a')
            sage: K.embeddings(GF(3))
            []
        """
        if K is self:
            f = self.pari_polynomial('y')
            # Compute the conjugates of Mod(x, f).
            conj = self.pari_nf().nfgaloisconj()
            # Convert these to conjugates of self.gen().
            P = self._pari_absolute_structure()[1].lift()
            conj = sorted([self(P(g.Mod(f))) for g in conj])
            v = [self.hom([e]) for e in conj]    # check=False here?
            put_natural_embedding_first(v)
            return Sequence(v, cr=(v != []), immutable=True,
                            check=False, universe=self.Hom(self))
        elif K.characteristic() != 0:
            return Sequence([], immutable=True, check=False, universe=self.Hom(K))

        f = self.defining_polynomial()
        r = sorted(f.roots(K, multiplicities=False))
        v = [self.hom([e], check=False) for e in r]
        # If there is an embedding that preserves variable names
        # then it is most natural, so we put it first.
        put_natural_embedding_first(v)
        return Sequence(v, cr=bool(v), immutable=True,
                        check=False, universe=self.Hom(K))

    def minkowski_embedding(self, B=None, prec=None):
        r"""
        Return an `n \times n` matrix over ``RDF`` whose columns are the images of the
        basis `\{1, \alpha, \dots, \alpha^{n-1}\}` of ``self`` over
        `\QQ` (as vector spaces), where here
        `\alpha` is the generator of ``self`` over
        `\QQ`, i.e. ``self.gen(0)``. If `B` is not ``None``, return
        the images of the vectors in `B` as the columns instead. If ``prec`` is
        not ``None``, use ``RealField(prec)`` instead of ``RDF``.

        This embedding is the so-called "Minkowski embedding" of a number
        field in `\RR^n`: given the `n` embeddings
        `\sigma_1, \dots, \sigma_n` of ``self`` in
        `\CC`, write `\sigma_1, \dots, \sigma_r`
        for the real embeddings, and
        `\sigma_{r+1}, \dots, \sigma_{r+s}` for choices of one of
        each pair of complex conjugate embeddings (in our case, we simply
        choose the one where the image of `\alpha` has positive
        real part). Here `(r,s)` is the signature of ``self``. Then the
        Minkowski embedding is given by

        .. MATH::

            x \mapsto ( \sigma_1(x), \dots,
                     \sigma_r(x), \sqrt{2}\Re(\sigma_{r+1}(x)),
                     \sqrt{2}\Im(\sigma_{r+1}(x)), \dots,
                     \sqrt{2}\Re(\sigma_{r+s}(x)),
                     \sqrt{2}\Im(\sigma_{r+s}(x)))

        Equivalently, this is an embedding of ``self`` in `\RR^n` so
        that the usual norm on `\RR^n` coincides with
        `|x| = \sum_i |\sigma_i(x)|^2` on ``self``.

        .. TODO::

            This could be much improved by implementing homomorphisms
            over VectorSpaces.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^3 + 2)
            sage: F.minkowski_embedding()
            [ 1.00000000000000 -1.25992104989487  1.58740105196820]
            [ 1.41421356237... 0.8908987181... -1.12246204830...]
            [0.000000000000000  1.54308184421...  1.94416129723...]
            sage: F.minkowski_embedding([1, alpha+2, alpha^2-alpha])
            [ 1.00000000000000 0.740078950105127  2.84732210186307]
            [ 1.41421356237...  3.7193258428... -2.01336076644...]
            [0.000000000000000  1.54308184421... 0.40107945302...]
            sage: F.minkowski_embedding() * (alpha + 2).vector().column()
            [0.740078950105127]
            [ 3.7193258428...]
            [ 1.54308184421...]
        """
        n = self.degree()
        if prec is None:
            R = sage.rings.real_double.RDF
        else:
            R = sage.rings.real_mpfr.RealField(prec)
        r, s = self.signature()
        places = self.places(prec=prec)

        if B is None:
            B = [(self.gen(0))**i for i in range(n)]

        A = ZZ['x']
        f = A.gen(0)**2 - 2
        sqrt2 = f.roots(R)[1][0]

        d = {}

        for col in range(n):

            for row in range(r):
                d[(row, col)] = places[row](B[col])

            for i in range(s):
                z = places[r+i](B[col])
                d[(r+2*i, col)] = z.real()*sqrt2
                d[(r+2*i+1, col)] = z.imag()*sqrt2

        return sage.matrix.all.matrix(d)

    def logarithmic_embedding(self, prec=53):
        r"""
        Return the morphism of ``self`` under the logarithmic embedding
        in the category Set.

        The logarithmic embedding is defined as a map from the number field
        ``self`` to `\RR^n`.

        It is defined under Definition 4.9.6 in [Coh1993]_.

        INPUT:

        - ``prec`` -- desired floating point precision

        OUTPUT: the morphism of ``self`` under the logarithmic embedding in the
        category Set

        EXAMPLES::

            sage: CF.<a> = CyclotomicField(5)
            sage: f = CF.logarithmic_embedding()
            sage: f(0)
            (-1, -1)
            sage: f(7)
            (3.89182029811063, 3.89182029811063)

        ::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 + 5)
            sage: f = K.logarithmic_embedding()
            sage: f(0)
            (-1, -1)
            sage: f(7)
            (1.94591014905531, 3.89182029811063)

        ::

            sage: F.<a> = NumberField(x^4 - 8*x^2 + 3)
            sage: f = F.logarithmic_embedding()
            sage: f(0)
            (-1, -1, -1, -1)
            sage: f(7)
            (1.94591014905531, 1.94591014905531, 1.94591014905531, 1.94591014905531)
        """
        def closure_map(x, prec=53):
            """
            The function closure of the logarithmic embedding.
            """
            K = self
            K_embeddings = K.places(prec)
            r1, r2 = K.signature()
            r = r1 + r2 - 1

            from sage.rings.real_mpfr import RealField
            Reals = RealField(prec)

            if x == 0:
                return vector([-1 for _ in range(r + 1)])

            x_logs = []
            for i in range(r1):
                sigma = K_embeddings[i]
                x_logs.append(Reals(abs(sigma(x))).log())
            for i in range(r1, r + 1):
                tau = K_embeddings[i]
                x_logs.append(2 * Reals(abs(tau(x))).log())

            return vector(x_logs)

        hom = Hom(self, VectorSpace(RR, len(closure_map(self(0), prec))), Sets())
        return hom(closure_map)

    def places(self, all_complex=False, prec=None):
        r"""
        Return the collection of all infinite places of ``self``.

        By default, this returns the set of real places as
        homomorphisms into ``RIF`` first, followed by a choice of one of
        each pair of complex conjugate homomorphisms into ``CIF``.

        On the other hand, if ``prec`` is not ``None``, we simply return places
        into ``RealField(prec)`` and ``ComplexField(prec)`` (or ``RDF``, ``CDF`` if
        ``prec=53``). One can also use ``prec=infinity``, which returns embeddings
        into the field `\overline{\QQ}` of algebraic numbers (or its subfield
        `\mathbb{A}` of algebraic reals); this permits exact computation, but
        can be extremely slow.

        There is an optional flag ``all_complex``, which defaults to ``False``. If
        ``all_complex`` is ``True``, then the real embeddings are returned as
        embeddings into ``CIF`` instead of ``RIF``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^3 - 100*x + 1); F.places()
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

            sage: F.<alpha> = NumberField(x^3 + 7); F.places()
            [Ring morphism:
               From: Number Field in alpha with defining polynomial x^3 + 7
               To:   Real Field with 106 bits of precision
               Defn: alpha |--> -1.912931182772389101199116839549,
             Ring morphism:
               From: Number Field in alpha with defining polynomial x^3 + 7
               To:   Complex Field with 53 bits of precision
               Defn: alpha |--> 0.956465591386195 + 1.65664699997230*I]

        ::

            sage: F.<alpha> = NumberField(x^3 + 7) ; F.places(all_complex=True)
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

        elif prec == Infinity:
            from sage.rings.qqbar import AA, QQbar
            R = AA
            C = QQbar

        else:
            R = sage.rings.real_mpfr.RealField(prec)
            C = sage.rings.complex_mpfr.ComplexField(prec)

        # first, find the intervals with roots, and see how much
        # precision we need to approximate the roots
        #
        all_intervals = [x[0] for x in self.defining_polynomial().roots(C)]

        # first, set up the real places
        if all_complex:
            real_intervals = [x for x in all_intervals if x.imag().is_zero()]
        else:
            real_intervals = [x[0] for x in self.defining_polynomial().roots(R)]

        if prec is None:
            real_places = [self.hom([i.center()], check=False)
                           for i in real_intervals]

            complex_places = [self.hom([i.center()], check=False)
                              for i in all_intervals if i.imag() > 0]
        else:
            real_places = [self.hom([i], check=False) for i in real_intervals]

            complex_places = [self.hom([i], check=False)
                              for i in all_intervals if i.imag() > 0]

        return real_places + complex_places

    def real_places(self, prec=None):
        """
        Return all real places of ``self`` as homomorphisms into ``RIF``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<alpha> = NumberField(x^4 - 7) ; F.real_places()
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

    def abs_val(self, v, iota, prec=None):
        r"""
        Return the value `|\iota|_{v}`.

        INPUT:

        - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of ``K.places(prec)``)
        - ``iota`` -- an element of ``K``
        - ``prec`` -- (default: ``None``) the precision of the real field

        OUTPUT: the absolute value as a real number

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<xi> = NumberField(x^3 - 3)
            sage: phi_real = K.places()[0]
            sage: phi_complex = K.places()[1]
            sage: v_fin = tuple(K.primes_above(3))[0]

            sage: K.abs_val(phi_real, xi^2)
            2.08008382305190

            sage: K.abs_val(phi_complex, xi^2)
            4.32674871092223

            sage: K.abs_val(v_fin, xi^2)
            0.111111111111111

        Check that :issue:`28345` is fixed::

            sage: K.abs_val(v_fin, K.zero())
            0.000000000000000
        """
        if prec is None:
            prec = 53
        R = sage.rings.real_mpfr.RealField(prec)
        if iota == 0:
            return R.zero()
        try:
            p = v.smallest_integer()
            iota_ideal = self.ideal(self(iota))
            exponent = - v.residue_class_degree() * iota_ideal.valuation(v)
            return R(p**exponent)
        except AttributeError:
            if is_real_place(v):
                return R(v(iota).abs())
            else:
                return R(v(iota).abs()**2)

    def relativize(self, alpha, names, structure=None):
        r"""
        Given an element in ``self`` or an embedding of a subfield into ``self``,
        return a relative number field `K` isomorphic to ``self`` that is relative
        over the absolute field `\QQ(\alpha)` or the domain of `\alpha`, along
        with isomorphisms from `K` to ``self`` and from ``self`` to `K`.

        INPUT:

        - ``alpha`` -- an element of ``self``  or an embedding of a subfield into
          ``self``
        - ``names`` -- 2-tuple of names of generator for output field `K` and the
          subfield `\QQ(\alpha)`
        - ``structure`` -- an instance of
          :class:`structure.NumberFieldStructure` or ``None`` (default:
          ``None``), if ``None``, then the resulting field's :meth:`structure`
          will return isomorphisms from and to this field. Otherwise, the field
          will be equipped with ``structure``.

        OUTPUT: `K` -- relative number field

        Also, ``K.structure()`` returns ``from_K`` and ``to_K``, where
        ``from_K`` is an isomorphism from `K` to ``self`` and ``to_K`` is an isomorphism
        from ``self`` to `K`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^10 - 2)
            sage: L.<c,d> = K.relativize(a^4 + a^2 + 2); L
            Number Field in c with defining polynomial
             x^2 - 1/5*d^4 + 8/5*d^3 - 23/5*d^2 + 7*d - 18/5 over its base field
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
            Number Field in b with defining polynomial x^2 - a0 + 2 over its base field
            sage: L2 = K.relativize(tau, 'b'); L2
            Number Field in b with defining polynomial x^2 + a0 over its base field
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
            [Number Field in c with defining polynomial x^4 + 1 over its base field,
             Number Field in c with defining polynomial x^2 - a1*x + 1 over its base field,
             Number Field in c with defining polynomial x^2 - 1/2*a2 over its base field,
             Number Field in c with defining polynomial x^2 - a3*x - 1 over its base field,
             Number Field in c with defining polynomial x - a4 over its base field]

        We can relativize over a relative field::

            sage: K.<z> = CyclotomicField(16)
            sage: L, L_into_K, _ = K.subfields(4)[0]; L
            Number Field in z0 with defining polynomial x^4 + 16
             with z0 = 1.414213562373095? + 1.414213562373095?*I
            sage: F, F_into_L, _ = L.subfields(2)[0]; F
            Number Field in z0_0 with defining polynomial x^2 + 64 with z0_0 = 8*I

            sage: L_over_F = L.relativize(F_into_L, 'c'); L_over_F
            Number Field in c with defining polynomial x^2 - 1/2*z0_0 over its base field
            sage: L_over_F_into_L, _ = L_over_F.structure()

            sage: K_over_rel = K.relativize(L_into_K * L_over_F_into_L, 'a'); K_over_rel
            Number Field in a with defining polynomial x^2 - 1/2*c over its base field
            sage: K_over_rel.base_field() is L_over_F
            True
            sage: K_over_rel.structure()
            (Relative number field morphism:
              From: Number Field in a with defining polynomial x^2 - 1/2*c over its base field
              To:   Cyclotomic Field of order 16 and degree 8
              Defn: a |--> z
                    c |--> 2*z^2
                    z0_0 |--> 8*z^4, Ring morphism:
              From: Cyclotomic Field of order 16 and degree 8
              To:   Number Field in a with defining polynomial x^2 - 1/2*c over its base field
              Defn: z |--> a)

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

        Only one name is required when a morphism is given (fixing :issue:`12005`)::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: L.<b> = NumberField(x^4 - x^2 + 1)
            sage: phi = K.hom(b^3, L)
            sage: M.<r> = L.relativize(phi)
            sage: M
            Number Field in r with defining polynomial x^2 - i*x - 1 over its base field
            sage: M.base_field()
            Number Field in i with defining polynomial x^2 + 1

        See :issue:`27469`::

            sage: L.<z24> = CyclotomicField(24)
            sage: K.<z8> = L.subfield(z24^3)[0]
            sage: L.relativize(K.hom(L), L.variable_name()+'0' )
            Number Field in z2400 with defining polynomial x^2 + z240^3*x - z240^2 over its base field
        """
        # step 1: construct the abstract field generated by alpha.w
        # step 2: make a relative extension of it.
        # step 3: construct isomorphisms
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        from sage.categories.map import Map
        if isinstance(alpha, Map):
            # alpha better be a morphism with codomain self
            if alpha.codomain() != self:
                raise ValueError("Co-domain of morphism must be self")
            L = alpha.domain()
            alpha = alpha(L.gen())  # relativize over phi's domain
            if L is QQ:
                from sage.rings.polynomial.polynomial_ring import polygen
                f = polygen(QQ)
            else:
                f = L.defining_polynomial()  # = alpha.minpoly()
            names = normalize_names(len(names), names)
        else:
            # alpha must be an element coercible to self
            alpha = self(alpha)
            f = alpha.minpoly()
            names = normalize_names(2, names)
            L = NumberField(f, names[1])

        # now we do some linear algebra to find the minpoly of self.gen() over L
        L_into_self = L.hom([alpha])

        extdeg = self.absolute_degree() // L.absolute_degree()  # [ L : self ]
        a = self.gen()

        # we will find a linear relation between small powers of a over L
        basis = [a**i * b for i in range(extdeg)
                 for b in map(L_into_self, L.power_basis())]
        basis.append(a**extdeg)  # this one makes the basis no longer a basis
        mat = matrix([b.vector() for b in basis])
        soln_space = mat.left_kernel(mat.row_space()(0))
        # the solution space is one dimensional and the last entry is nonzero
        # because a satisfies no smaller linear relation
        assert soln_space.dimension() == 1
        (reln, ) = soln_space.basis()
        assert reln[-1] != 0
        reln = reln * ~reln[-1]

        # now we need to get those coeffs in L
        coeff_mat = matrix(extdeg, f.degree(), list(reln)[:-1])  # easy way to divide into the correct lengths
        coeffs_in_L = [r * vector(L.power_basis()) for r in coeff_mat.rows()]
        # f is the minimal polynomial of a over L
        f = L['x'](coeffs_in_L + [1])
        # sanity check...

        mp_in_self = self['x']([L_into_self(_) for _ in f.coefficients(sparse=False)])
        assert mp_in_self(a) == 0

        if structure is None:
            from sage.rings.number_field.structure import RelativeFromAbsolute
            structure = RelativeFromAbsolute(self, alpha)
        if L is QQ:
            return L.extension(f, names[0])
        else:
            return L.extension(f, names[0], structure=structure)

    # Synonyms so that terminology appropriate to relative number fields
    # can be applied to an absolute number field:

    def absolute_degree(self):
        """
        A synonym for :meth:`degree`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_degree()
            2
        """
        return self.degree()

    def relative_degree(self):
        """
        A synonym for :meth:`degree`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_degree()
            2
        """
        return self.degree()

    def relative_polynomial(self):
        """
        A synonym for :meth:`polynomial`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_polynomial()
            x^2 + 1
        """
        return self.polynomial()

    def relative_vector_space(self, *args, **kwds):
        """
        A synonym for :meth:`vector_space`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
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
        return self.free_module(*args, **kwds)

    def absolute_discriminant(self):
        """
        A synonym for :meth:`discriminant`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_discriminant()
            -4
        """
        return self.discriminant()

    def relative_discriminant(self):
        """
        A synonym for :meth:`discriminant`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_discriminant()
            -4
        """
        return self.discriminant()

    def absolute_different(self):
        """
        A synonym for :meth:`different`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.absolute_different()
            Fractional ideal (2)
        """
        return self.different()

    def relative_different(self):
        """
        A synonym for :meth:`different`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.relative_different()
            Fractional ideal (2)
        """
        return self.different()

    def hilbert_symbol(self, a, b, P=None):
        r"""
        Return the Hilbert symbol `(a,b)_P` for a prime `P` of ``self``
        and nonzero elements `a` and `b` of ``self``.

        If `P` is omitted, return the global Hilbert symbol `(a,b)` instead.

        INPUT:

        - ``a``, ``b`` -- elements of ``self``

        - ``P`` -- (default: ``None``) if ``None``, compute the global
          symbol.  Otherwise, `P` should be either a prime ideal of ``self``
          (which may also be given as a generator or set of generators)
          or a real or complex embedding.

        OUTPUT: if `a` or `b` is zero, returns 0

        If `a` and `b` are nonzero and `P` is specified, returns
        the Hilbert symbol `(a,b)_P`, which is `1` if the equation
        `a x^2 + b y^2 = 1` has a solution in the completion of
        ``self`` at `P`, and is `-1` otherwise.

        If `a` and `b` are nonzero and `P` is unspecified, returns `1`
        if the equation has a solution in ``self`` and `-1` otherwise.

        EXAMPLES:

        Some global examples::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 - 23)
            sage: K.hilbert_symbol(0, a + 5)
            0
            sage: K.hilbert_symbol(a, 0)
            0
            sage: K.hilbert_symbol(-a, a + 1)
            1
            sage: K.hilbert_symbol(-a, a + 2)
            -1
            sage: K.hilbert_symbol(a, a + 5)
            -1

        That the latter two are unsolvable should be visible in local
        obstructions.  For the first, this is a prime ideal above 19.
        For the second, the ramified prime above 23::

            sage: K.hilbert_symbol(-a, a + 2, a + 2)
            -1
            sage: K.hilbert_symbol(a, a + 5, K.ideal(23).factor()[0][0])
            -1

        More local examples::

            sage: K.hilbert_symbol(a, 0, K.fractional_ideal(5))
            0
            sage: K.hilbert_symbol(a, a + 5, K.fractional_ideal(5))
            1
            sage: K.hilbert_symbol(a + 1, 13, (a+6)*K)
            -1
            sage: [emb1, emb2] = K.embeddings(AA)
            sage: K.hilbert_symbol(a, -1, emb1)
            -1
            sage: K.hilbert_symbol(a, -1, emb2)
            1

        Ideals P can be given by generators::

            sage: K.<a> = NumberField(x^5 - 23)
            sage: pi = 2*a^4 + 3*a^3 + 4*a^2 + 15*a + 11
            sage: K.hilbert_symbol(a, a + 5, pi)
            1
            sage: rho = 2*a^4 + 3*a^3 + 4*a^2 + 15*a + 11
            sage: K.hilbert_symbol(a, a + 5, rho)
            1

        This also works for non-principal ideals::

            sage: K.<a> = QuadraticField(-5)
            sage: P = K.ideal(3).factor()[0][0]
            sage: P.gens_reduced()  # random, could be the other factor
            (3, a + 1)
            sage: K.hilbert_symbol(a, a + 3, P)
            1
            sage: K.hilbert_symbol(a, a + 3, [3, a+1])
            1

        Primes above 2::

            sage: K.<a> = NumberField(x^5 - 23)
            sage: p = [p[0] for p in (2*K).factor() if p[0].norm() == 16][0]
            sage: K.hilbert_symbol(a, a + 5, p)
            1
            sage: K.hilbert_symbol(a, 2, p)
            1
            sage: K.hilbert_symbol(-1, a - 2, p)
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
                    with a = 2.236067977499790?
              To:   Complex Field with 53 bits of precision
              Defn: a |--> -2.23606797749979)
            given as complex embedding in hilbert_symbol. Is it real or complex?
            sage: K.hilbert_symbol(-1, -1, K.embeddings(QQbar)[0])
            Traceback (most recent call last):
            ...
            ValueError: Possibly real place (=Ring morphism:
              From: Number Field in a with defining polynomial x^2 - 5
                    with a = 2.236067977499790?
              To:   Algebraic Field
              Defn: a |--> -2.236067977499790?)
            given as complex embedding in hilbert_symbol. Is it real or complex?
            sage: K.<b> = QuadraticField(-5)
            sage: K.hilbert_symbol(-1, -1, K.embeddings(CDF)[0])
            1
            sage: K.hilbert_symbol(-1, -1, K.embeddings(QQbar)[0])
            1

        `a` and `b` do not have to be integral or coprime::

            sage: K.<i> = QuadraticField(-1)
            sage: K.hilbert_symbol(1/2, 1/6, 3*K)
            1
            sage: p = 1 + i
            sage: K.hilbert_symbol(p, p, p)
            1
            sage: K.hilbert_symbol(p, 3*p, p)
            -1
            sage: K.hilbert_symbol(3, p, p)
            -1
            sage: K.hilbert_symbol(1/3, 1/5, 1 + i)
            1
            sage: L = QuadraticField(5, 'a')
            sage: L.hilbert_symbol(-3, -1/2, 2)
            1

        Various other examples::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: K.hilbert_symbol(-6912, 24, -a^2 - a - 2)
            1
            sage: K.<a> = NumberField(x^5 - 23)
            sage: P = K.ideal(-1105*a^4 + 1541*a^3 - 795*a^2 - 2993*a + 11853)
            sage: Q = K.ideal(-7*a^4 + 13*a^3 - 13*a^2 - 2*a + 50)
            sage: b = -a+5
            sage: K.hilbert_symbol(a, b, P)
            1
            sage: K.hilbert_symbol(a, b, Q)
            1
            sage: K.<a> = NumberField(x^5 - 23)
            sage: P = K.ideal(-1105*a^4 + 1541*a^3 - 795*a^2 - 2993*a + 11853)
            sage: K.hilbert_symbol(a, a + 5, P)
            1
            sage: K.hilbert_symbol(a, 2, P)
            1
            sage: K.hilbert_symbol(a + 5, 2, P)
            -1
            sage: K.<a> = NumberField(x^3 - 4*x + 2)
            sage: K.hilbert_symbol(2, -2, K.primes_above(2)[0])
            1

        Check that the bug reported at :issue:`16043` has been fixed::

            sage: K.<a> = NumberField(x^2 + 5)
            sage: p = K.primes_above(2)[0]; p
            Fractional ideal (2, a + 1)
            sage: K.hilbert_symbol(2*a, -1, p)
            1
            sage: K.hilbert_symbol(2*a, 2, p)
            -1
            sage: K.hilbert_symbol(2*a, -2, p)
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

        from sage.categories.map import Map
        from sage.categories.rings import Rings
        if isinstance(P, Map) and P.category_for().is_subcategory(Rings()):
            # P is a morphism of Rings
            if P.domain() is not self:
                raise ValueError("Domain of P (=%s) should be self (=%s) in self.hilbert_symbol" % (P, self))
            codom = P.codomain()
            from sage.rings.qqbar import AA, QQbar
            if isinstance(codom, (sage.rings.abc.ComplexField, sage.rings.abc.ComplexDoubleField, sage.rings.abc.ComplexIntervalField)) or \
               codom is QQbar:
                if P(self.gen()).imag() == 0:
                    raise ValueError("Possibly real place (=%s) given as complex embedding in hilbert_symbol. Is it real or complex?" % P)
                return 1
            if isinstance(codom, (sage.rings.abc.RealField, sage.rings.abc.RealDoubleField)) or codom is AA:
                if P(a) > 0 or P(b) > 0:
                    return 1
                return -1
        if not isinstance(P, NumberFieldIdeal):
            P = self.ideal(P)
        if P.number_field() is not self:
            raise ValueError("P (=%s) should be an ideal of self (=%s) in hilbert_symbol, not of %s" % (P, self, P.number_field()))
        if not P.is_prime():
            raise ValueError("Non-prime ideal P (=%s) in hilbert_symbol" % P)
        return pari(self).nfhilbert(a, b, P.pari_prime())

    def hilbert_symbol_negative_at_S(self, S, b, check=True):
        """
        Return `a` such that the Hilbert conductor of `a` and `b` is `S`.

        INPUT:

        - ``S`` -- list of places (or prime ideals) of even cardinality
        - ``b`` -- a nonzero rational number which is a non-square locally
          at every place in `S`
        - ``check`` -- boolean (default: ``True``); perform additional checks
          on the input and confirm the output

        OUTPUT:

        - an element `a` that has negative Hilbert symbol `(a,b)_p` for
          every (finite and infinite) place `p` in `S`.

        ALGORITHM:

        The implementation is following algorithm 3.4.1 in [Kir2016]_.
        We note that class and unit groups are computed using the generalized
        Riemann hypothesis. If it is false, this may result in an infinite loop.
        Nevertheless, if the algorithm terminates the output is correct.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 20072)
            sage: S = [K.primes_above(3)[0], K.primes_above(23)[0]]
            sage: b = K.hilbert_symbol_negative_at_S(S, a + 1)
            sage: [K.hilbert_symbol(b, a + 1, p) for p in S]
            [-1, -1]
            sage: K.<d> = CyclotomicField(11)
            sage: S = [K.primes_above(2)[0], K.primes_above(11)[0]]
            sage: b = d + 5
            sage: a = K.hilbert_symbol_negative_at_S(S, b)
            sage: [K.hilbert_symbol(a,b,p) for p in S]
            [-1, -1]
            sage: k.<c> = K.maximal_totally_real_subfield()[0]
            sage: S = [k.primes_above(3)[0], k.primes_above(5)[0]]
            sage: S += k.real_places()[:2]
            sage: b = 5 + c + c^9
            sage: a = k.hilbert_symbol_negative_at_S(S, b)
            sage: [k.hilbert_symbol(a, b, p) for p in S]
            [-1, -1, -1, -1]

        Note that the closely related Hilbert conductor
        takes only the finite places into account::

            sage: k.hilbert_conductor(a, b)
            Fractional ideal (15)


        TESTS::

            sage: K.<a> = NumberField(x^2 + 20072)
            sage: S = [K.primes_above(3)[0], K.primes_above(23)[0]]
            sage: s = K.hilbert_symbol_negative_at_S(S, a + 1)

        AUTHORS:

        - Simon Brandhorst, Anna Haensch (01-05-2018)
        """
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        from sage.modules.free_module import VectorSpace
        from sage.matrix.constructor import matrix
        from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup

        # input checks
        if not isinstance(S, list):
            raise TypeError("first argument must be a list")
        if b not in self:
            raise TypeError("second argument must be an element of this field")
        b = self(b)
        if b == 0:
            raise ValueError("second argument must be nonzero")
        if len(S) % 2:
            raise ValueError("the list should be of even cardinality")
        if check:
            for p in S:
                if p.parent() is self.ideal_monoid():
                    if not p.is_prime():
                        raise ValueError("not a prime ideal")
                    if self.quadratic_defect(b, p) == infinity.Infinity:
                        raise ValueError(f"{b} is a square in the completion "
                                         f"with respect to {p}")
                else:
                    if p not in self.real_places():
                        raise ValueError("entries of the list must be "
                                         "prime ideals or real places")
                    if p(b) > 0:
                        raise ValueError(f"{b} is a square in the completion "
                                         f"with respect to {p}")

        # L is the list of primes that we need to consider, b must have
        # nonzero valuation for each prime in L, this is the set S'
        # in Kirschmer's algorithm
        L = []
        for P in self.prime_factors(b) + self.prime_factors(self(2)):
            if not (P in S or P in L):
                L.append(P)

        # This adds some infinite places to L
        L.extend(sigma for sigma in self.real_places()
                 if sigma(b) < 0 and sigma not in S)
        Cl = self.class_group(proof=False)
        U = self.unit_group(proof=False).gens()
        SL = S + L
        # the finite places in SL
        P = [p for p in SL if p.parent() is self.ideal_monoid()]

        # v is the vector that we are searching for.
        # It represents the case when the Hilbert
        # symbol is negative for all primes in S and positive
        # at all primes in S'
        # For technical reasons, a Hilbert symbol of -1 is
        # respresented as 1 and a Hilbert symbol of 1
        # is represented as 0
        V = VectorSpace(GF(2), len(SL))
        v = V([1]*len(S) + [0]*len(L))

        # The algorithm terminates when the vector v is in the
        # subspace of V generated by the image of the phi map
        # on the set of generators

        def phi(x):
            return V([(1 - self.hilbert_symbol(x, b, p)) // 2 for p in SL])

        M = matrix([phi(g) for g in U])

        # we have to work around the inconvenience that multiplicative
        # abelian groups in sage do not yet have homomorphisms...
        Cl_additive = AdditiveAbelianGroup(Cl.gens_orders())
        n = len(Cl_additive.gens())
        A = AdditiveAbelianGroup([0]*(len(P)+1))
        PCl = []
        for p in P:
            pr = Cl(p).exponents()
            pr = Cl_additive.sum([Cl_additive.gens()[i]*pr[i] for i in range(n)])
            PCl.append(pr)

        # search through all the primes
        found = False
        for q in sage.sets.primes.Primes():
            for Q in self.primes_above(q):
                if Q in P:
                    continue
                pr = Cl(Q).exponents()
                pr = Cl_additive.sum([Cl_additive.gens()[i]*pr[i] for i in range(n)])
                # compute the kernel
                K = A.hom(PCl + [pr]).kernel().gens()
                K = [A(k) for k in K]
                # generate it by a list of ideals
                Pq = P + [Q]
                K = [prod([Pq[i]**k[i] for i in range(len(Pq))]) for k in K]
                # the ideals are principal
                # find a single generator for each
                K = [k.gens_reduced(proof=False)[0] for k in K]
                # we can now apply phi
                W = M.stack(matrix([phi(g) for g in K]))
                if v in W.row_space():
                    found = True
                    break
            if found:
                break
        J = list(U) + K
        l = W.solve_left(v)
        a = prod([J[i]**l[i] for i in range(len(J))])
        # let us double check the result
        if check:
            assert phi(a) == v, "oops"
        return a

    def hilbert_conductor(self, a, b):
        """
        This is the product of all (finite) primes where the Hilbert symbol is `-1`.
        What is the same, this is the (reduced) discriminant of the quaternion
        algebra `(a,b)` over a number field.

        INPUT:

        - ``a``, ``b`` -- elements of the number field ``self``

        OUTPUT: squarefree ideal of the ring of integers of ``self``

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F.<a> = NumberField(x^2 - x - 1)
            sage: F.hilbert_conductor(2*a, F(-1))
            Fractional ideal (2)
            sage: K.<b> = NumberField(x^3 - 4*x + 2)
            sage: K.hilbert_conductor(K(2), K(-2))
            Fractional ideal (1)
            sage: K.hilbert_conductor(K(2*b), K(-2))
            Fractional ideal (b^2 + b - 2)

        AUTHOR:

        - Aly Deines
        """
        a, b = self(a), self(b)
        d = self.ideal(1)
        for p in set(self.ideal(2).prime_factors()).union(self.ideal(a).prime_factors()).union(self.ideal(b).prime_factors()):
            if self.hilbert_symbol(a, b, p) == -1:
                d *= p
        return d

    def elements_of_bounded_height(self, **kwds):
        r"""
        Return an iterator over the elements of ``self`` with relative
        multiplicative height at most ``bound``.

        This algorithm computes 2 lists: `L` containing elements `x` in `K` such that
        `H_k(x) \leq B`, and a list `L'` containing elements `x` in `K` that, due to
        floating point issues,
        may be slightly larger than the bound. This can be controlled
        by lowering the tolerance.

        In the current implementation, both lists `(L,L')` are merged and returned in
        form of iterator.

        ALGORITHM:

        This is an implementation of the revised algorithm (Algorithm 4) in
        [DK2013]_. Algorithm 5 is used for imaginary quadratic fields.

        INPUT: keyword arguments:

        - ``bound`` -- a real number

        - ``tolerance`` -- (default: 0.01) a rational number in `(0,1]`

        - ``precision`` -- (default: 53) a positive integer

        OUTPUT: an iterator of number field elements

        EXAMPLES:

        There are no elements in a number field with multiplicative height less
        than 1::

            sage: x = polygen(QQ, 'x')
            sage: K.<g> = NumberField(x^5 - x + 19)
            sage: list(K.elements_of_bounded_height(bound=0.9))
            []

        The only elements in a number field of height 1 are 0 and the roots of
        unity::

            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: list(K.elements_of_bounded_height(bound=1))
            [0, a + 1, a, -1, -a - 1, -a, 1]

        ::

            sage: K.<a> = CyclotomicField(20)
            sage: len(list(K.elements_of_bounded_height(bound=1)))
            21

        The elements in the output iterator all have relative multiplicative
        height at most the input bound::

            sage: K.<a> = NumberField(x^6 + 2)
            sage: L = K.elements_of_bounded_height(bound=5)
            sage: for t in L:
            ....:     exp(6*t.global_height())
            1.00000000000000
            1.00000000000000
            1.00000000000000
            2.00000000000000
            2.00000000000000
            2.00000000000000
            2.00000000000000
            4.00000000000000
            4.00000000000000
            4.00000000000000
            4.00000000000000

        ::

            sage: K.<a> = NumberField(x^2 - 71)
            sage: L = K.elements_of_bounded_height(bound=20)
            sage: all(exp(2*t.global_height()) <= 20 for t in L) # long time (5 s)
            True

        ::

            sage: K.<a> = NumberField(x^2 + 17)
            sage: L = K.elements_of_bounded_height(bound=120)
            sage: len(list(L))
            9047

        ::

            sage: K.<a> = NumberField(x^4 - 5)
            sage: L = K.elements_of_bounded_height(bound=50)
            sage: len(list(L)) # long time (2 s)
            2163

        ::

            sage: K.<a> = CyclotomicField(13)
            sage: L = K.elements_of_bounded_height(bound=2)
            sage: len(list(L)) # long time (3 s)
            27

        ::

            sage: K.<a> = NumberField(x^6 + 2)
            sage: L = K.elements_of_bounded_height(bound=60, precision=100)
            sage: len(list(L)) # long time (5 s)
            1899

        ::

            sage: K.<a> = NumberField(x^4 - x^3 - 3*x^2 + x + 1)
            sage: L = K.elements_of_bounded_height(bound=10, tolerance=0.1)
            sage: len(list(L))
            99

        AUTHORS:

        - John Doyle (2013)

        - David Krumm (2013)

        - Raman Raghukul (2018)
        """
        from sage.rings.number_field.bdd_height import bdd_height, bdd_height_iq
        r1, r2 = self.signature()
        r = r1 + r2 - 1
        B = kwds.pop('bound')
        if self.degree() == 2 and r == 0:
            return bdd_height_iq(self, B)
        else:
            tol = kwds.pop('tolerance', 1e-2)
            prec = kwds.pop('precision', 53)
            return bdd_height(self, B, tolerance=tol, precision=prec)

    def _factor_univariate_polynomial(self, poly, **kwargs):
        """
        Factorisation of univariate polynomials over absolute number fields.

        This is called by the ``factor`` method of univariate polynomials.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x**2 + 1)
            sage: x = polygen(K, 'x')
            sage: factor(x*x + 4)  # indirect doctest
            (x - 2*i) * (x + 2*i)

        TESTS::

            sage: with seed(0):
            ....:     P.<y> = NumberField(QQ['a'].random_element(30,10^100), 'a')[]
            sage: (3*y^4/4).factor()
            (3/4) * y^4
        """
        if self.degree() == 1:
            factors = poly.change_ring(QQ).factor()
            return Factorization([(p.change_ring(self), e)
                                  for p, e in factors], self(factors.unit()))
        elif poly.is_term():
            return Factorization([(poly.parent().gen(), poly.degree())],
                                 poly.leading_coefficient())

        # Convert the polynomial we want to factor to PARI
        f = poly._pari_with_name()
        try:
            # Try to compute the PARI nf structure with important=False.
            # This will raise RuntimeError if the computation is too
            # difficult.
            Rpari = self.pari_nf(important=False)
        except RuntimeError:
            # Cannot easily compute the nf structure, use the defining
            # polynomial instead.
            Rpari = self.pari_polynomial("y")
        G = list(Rpari.nffactor(f))
        # PARI's nffactor() ignores the unit, _factor_pari_helper()
        # adds back the unit of the factorization.
        return poly._factor_pari_helper(G)


class NumberField_cyclotomic(NumberField_absolute, sage.rings.abc.NumberField_cyclotomic):
    """
    Create a cyclotomic extension of the rational field.

    The command ``CyclotomicField(n)`` creates the `n`-th cyclotomic field,
    obtained by adjoining an `n`-th root of unity to the rational field.

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

        sage: cf12 = CyclotomicField(12)
        sage: z12 = cf12.0
        sage: cf6 = CyclotomicField(6)
        sage: z6 = cf6.0
        sage: FF = Frac(cf12['x'])
        sage: x = FF.0
        sage: z6*x^3/(z6 + x)
        zeta12^2*x^3/(x + zeta12^2)

    ::

        sage: cf6 = CyclotomicField(6); z6 = cf6.gen(0)
        sage: cf3 = CyclotomicField(3); z3 = cf3.gen(0)
        sage: cf3(z6)
        zeta3 + 1
        sage: cf6(z3)
        zeta6 - 1
        sage: type(cf6(z3))
        <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
        sage: cf1 = CyclotomicField(1); z1 = cf1.0
        sage: cf3(z1)
        1
        sage: type(cf3(z1))
        <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
    """
    def __init__(self, n, names, embedding=None, assume_disc_small=False, maximize_at_primes=None):
        """
        A cyclotomic field, i.e., a field obtained by adjoining an `n`-th
        root of unity to the rational numbers.

        EXAMPLES::

            sage: k = CyclotomicField(3)
            sage: type(k)
            <class 'sage.rings.number_field.number_field.NumberField_cyclotomic_with_category'>

        TESTS:

        The ``gcd`` and ``xgcd`` methods do not agree on this field, see
        :issue:`23274`::

            sage: TestSuite(k).run()
            Failure in _test_gcd_vs_xgcd:
            ...
            AssertionError:... The methods gcd and xgcd disagree on Cyclotomic Field of order 3 and degree 2:
              gcd(0,2) = 1
             xgcd(0,2) = (2, 0, 1)
            ------------------------------------------------------------
            The following tests failed: _test_gcd_vs_xgcd

        ::

            sage: type(CyclotomicField(4).zero())
            <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_gaussian'>
            sage: type(CyclotomicField(6).one())
            <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: type(CyclotomicField(6).an_element())
            <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: type(CyclotomicField(15).zero())
            <class 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>
        """
        f = QQ['x'].cyclotomic_polynomial(n)
        if names[0].startswith('zeta'):
            latex_name = "\\zeta_{%s}" % n
        else:
            latex_name = latex_variable_name(names[0])
        self.__n = n = Integer(n)
        NumberField_absolute.__init__(self, f,
                                      name=names,
                                      latex_name=latex_name,
                                      check=False,
                                      embedding=embedding,
                                      assume_disc_small=assume_disc_small,
                                      maximize_at_primes=maximize_at_primes)
        if n % 2:
            self.__zeta_order = 2 * n
        else:
            self.__zeta_order = n
        #  quadratic number fields require this:
        if f.degree() == 2:
            # define a boolean flag as for NumberField_quadratic to know, which
            # square root we choose (True means no embedding or positive
            # imaginary value).
            # Note that the test is done with NumberFieldElement and not with
            # NumberFieldElement_quadratic which requires somehow this flag.
            # As a consequence, a result of _an_element_() with the wrong class
            # is cached during the call to has_coerce_map_from. We reset the
            # cache afterwards.
            self._standard_embedding = not CDF.has_coerce_map_from(self) or CDF(self.gen()).imag() > 0
            self._cache_an_element = None

            if n == 4:
                self._element_class = number_field_element_quadratic.NumberFieldElement_gaussian
                self._D = ZZ(-1)
                self._NumberField_generic__gen = self._element_class(self, (QQ(0), QQ.one()))
            else:
                # n is 3 or 6
                self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic
                self._D = ZZ(-3)
                one_half = QQ((1, 2))
                if n == 3:
                    self._NumberField_generic__gen = self._element_class(self, (one_half-1, one_half))
                else:
                    self._NumberField_generic__gen = self._element_class(self, (one_half, one_half))

            # NumberField_absolute.__init__(...) set _zero_element and
            # _one_element to NumberFieldElement_absolute values, which is
            # wrong (and dangerous; such elements can actually be used to
            # crash Sage: see #5316).  Overwrite them with correct values.
            self._zero_element = self._element_class(self, (QQ(0), QQ(0)))
            self._one_element = self._element_class(self, (QQ.one(), QQ(0)))

        zeta = self.gen()
        zeta._set_multiplicative_order(n)
        self._init_embedding_approx()

    def construction(self):
        """
        Return data defining a functorial construction of ``self``.

        EXAMPLES::

            sage: F, R = CyclotomicField(5).construction()
            sage: R
            Rational Field
            sage: F.polys
            [x^4 + x^3 + x^2 + x + 1]
            sage: F.names
            ['zeta5']
            sage: F.embeddings
            [0.309016994374948? + 0.951056516295154?*I]
            sage: F.structures
            [None]
        """
        F, R = NumberField_generic.construction(self)
        F.cyclotomic = self.__n
        return F, R

    def _pushout_(self, other):
        r"""
        TESTS:

        Pushout is implemented for cyclotomic fields::

            sage: K.<a> = CyclotomicField(3)
            sage: L.<b> = CyclotomicField(4)
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(K,L,operator.add)
            Coercion on left operand via
                Generic morphism:
                  From: Cyclotomic Field of order 3 and degree 2
                  To:   Cyclotomic Field of order 12 and degree 4
                  Defn: a -> zeta12^2 - 1
            Coercion on right operand via
                Generic morphism:
                  From: Cyclotomic Field of order 4 and degree 2
                  To:   Cyclotomic Field of order 12 and degree 4
                  Defn: b -> zeta12^3
            Arithmetic performed after coercions.
            Result lives in Cyclotomic Field of order 12 and degree 4
            Cyclotomic Field of order 12 and degree 4

        As a consequence, operations work nicely::

            sage: a + b
            zeta12^3 + zeta12^2 - 1
            sage: a * b
            -zeta12
        """
        if isinstance(other, NumberField_cyclotomic):
            return CyclotomicField((self.__n).lcm(other.__n))

    def _magma_init_(self, magma):
        """
        Function returning a string to create this cyclotomic field in
        Magma.

        .. NOTE::

           The Magma generator name is also initialized to be the same
           as for the Sage field.

        EXAMPLES::

            sage: K = CyclotomicField(7,'z')
            sage: K._magma_init_(magma)                                # optional - magma
            'SageCreateWithNames(CyclotomicField(7),["z"])'
            sage: K = CyclotomicField(7,'zeta')
            sage: K._magma_init_(magma)                                # optional - magma
            'SageCreateWithNames(CyclotomicField(7),["zeta"])'
        """
        s = 'CyclotomicField(%s)' % self.__n
        return magma._with_names(s, self.variable_names())

    def _gap_init_(self):
        """
        Return a string that provides a representation of ``self`` in GAP.

        TESTS::

            sage: # needs sage.libs.gap
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
        GAP interface -- see :issue:`5618`. ::

            sage: # needs sage.groups
            sage: H = AlternatingGroup(4)
            sage: g = H((1,4,3))
            sage: K = H.subgroup([g])
            sage: z = CyclotomicField(3).an_element(); z
            zeta3
            sage: c = K.character([1,z,z**2]); c
            Character of Subgroup generated by [(1,4,3)] of
             (Alternating group of order 4!/2 as a permutation group)
            sage: c(g^2); z^2
            zeta3
            -zeta3 - 1
        """
        return 'CyclotomicField(%s)' % self.__n

    def _libgap_(self):
        """
        Return a LibGAP representation of ``self``.

        TESTS::

            sage: # needs sage.libs.gap
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
        n = self.__n
        return f"Cyclotomic Field of order {n} and degree {self.degree()}"

    def _n(self):
        """
        Return the `n` used to create this cyclotomic field.

        EXAMPLES::

            sage: CyclotomicField(3).zeta_order()
            6
            sage: CyclotomicField(3)._n()
            3
        """
        return self.__n

    def _latex_(self):
        r"""
        Return the latex representation of this cyclotomic field.

        EXAMPLES::

            sage: Z = CyclotomicField(4)
            sage: Z.gen()
            zeta4
            sage: latex(Z)  # indirect doctest
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

        We check that the bug reported on :issue:`8938` is fixed::

            sage: C5.<z> = CyclotomicField(5)
            sage: P.<s, t> = C5[]
            sage: f = (z^2 + z)*s
            sage: f
            (z^2 + z)*s
            sage: latex(f)
            \left(z^{2} + z\right) s
        """
        v = self.latex_variable_names()[0]
        if v.startswith('\\zeta_'):
            return "%s(%s)" % (latex(QQ), v)
        return NumberField_generic._latex_(self)

    def _coerce_map_from_(self, K):
        r"""
        Return a coercion map from `K` to ``self``, or ``None``.

        The cyclotomic field `\QQ(\zeta_n)` coerces into the
        cyclotomic field `\QQ(\zeta_m)` if and only if `n' \mid m`,
        where `n'` is the odd part of `n` if `4 \nmid n` and `n' = n`
        otherwise.

        The morphism is consistent with the chosen embedding into `\CC`.

        If `K` is not a cyclotomic field, the normal coercion rules
        for number fields are used.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(12)
            sage: L.<b> = CyclotomicField(132)
            sage: L.coerce_map_from(K)  # indirect doctest
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

        Check that :issue:`12632` is fixed::

            sage: K1 = CyclotomicField(1); K2 = CyclotomicField(2)
            sage: K1.coerce_map_from(K2)
            Generic morphism:
              From: Cyclotomic Field of order 2 and degree 1
              To:   Cyclotomic Field of order 1 and degree 1
              Defn: zeta2 -> -1

        Check that custom embeddings are respected (:issue:`13765`)::

            sage: # needs sage.symbolic
            sage: z105 = CDF(exp(2*pi*I/105))
            sage: Ka.<a> = CyclotomicField(105, embedding=z105^11)
            sage: Kb.<b> = CyclotomicField(35, embedding=z105^6)
            sage: Ka.coerce_map_from(Kb)
            Generic morphism:
              From: Cyclotomic Field of order 35 and degree 24
              To:   Cyclotomic Field of order 105 and degree 48
              Defn: b -> -a^44 - a^42 + a^39 + a^37 + a^35 - a^29 - a^27 - a^25 + a^24
                          - a^23 + a^22 - a^21 + a^20 + a^18 + a^16 - a^12 - a^10
                          - a^8 - a^6 + a^5 + a^3 + a
            sage: CC(b)
            0.936234870639737 + 0.351374824081343*I
            sage: CC(-a^44 - a^42 + a^39 + a^37 + a^35 - a^29 - a^27 - a^25 + a^24
            ....:     - a^23 + a^22 - a^21 + a^20 + a^18 + a^16 - a^12 - a^10
            ....:     - a^8 - a^6 + a^5 + a^3 + a)
            0.936234870639731 + 0.351374824081341*I

            sage: z15 = CDF(exp(2*pi*I/15))                                             # needs sage.symbolic
            sage: K6 = CyclotomicField(6, embedding=-z15^5)                             # needs sage.symbolic
            sage: CyclotomicField(15).coerce_map_from(K6)                               # needs sage.symbolic
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 15 and degree 8
              Defn: zeta6 -> -zeta15^5

            sage: CyclotomicField(15, embedding=z15^4).coerce_map_from(K6)              # needs sage.symbolic
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 15 and degree 8
              Defn: zeta6 -> -zeta15^5

        Check transitivity of coercion embeddings (:issue:`20513`)::

            sage: K60.<zeta60> = CyclotomicField(60)
            sage: K30.<zeta30> = CyclotomicField(30, embedding=zeta60**14)
            sage: K15.<zeta15> = CyclotomicField(15, embedding=zeta30**26)
            sage: K5.<zeta5> = CyclotomicField(5, embedding=zeta15**12)
            sage: K60.has_coerce_map_from(K5)
            True
            sage: K60(zeta5)
            -zeta60^14 - zeta60^12 + zeta60^6 + zeta60^4 - 1
            sage: _ == zeta60**(14*26*12)
            True
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
        Return an integer `e` such that `self.gen()^e == x`, or ``None``
        if no such integer exists. This is primarily used to construct
        embedding-respecting coercions.

        If `x` is complex, the result is either an integer `e` such
        that the absolute value of ``self.gen()^e - x`` is small or
        ``None`` if no such `e` is found.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(5)
            sage: K._log_gen(CDF(a))
            1
            sage: K._log_gen(CDF(a^4))
            4

            sage: # needs sage.symbolic
            sage: zeta105 = CC(exp(2*pi*i/105))
            sage: K.<a> = CyclotomicField(105, embedding=zeta105^13)
            sage: zeta105^13, CC(a)
            (0.712376096951345 + 0.701797902883992*I,
             0.712376096951345 + 0.701797902883991*I)
            sage: K._log_gen(zeta105^26)
            2
            sage: K._log_gen(zeta105)
            97
            sage: zeta105, CC(a^97)
            (0.998210129767735 + 0.0598041539450342*I,
             0.998210129767736 + 0.0598041539450313*I)
            sage: K._log_gen(zeta105^3)
            81
            sage: zeta105^3, CC(a)^81
            (0.983929588598630 + 0.178556894798637*I,
             0.983929588598631 + 0.178556894798635*I)

            sage: K.<a> = CyclotomicField(5, embedding=None)
            sage: K._log_gen(CDF(.5, -.8)) is None
            True

            sage: # needs sage.rings.padics
            sage: zeta5 = cyclotomic_polynomial(5).change_ring(Qp(11)).roots()[0][0]
            sage: zeta5 ^ 5
            1 + O(11^20)
            sage: K.<a> = CyclotomicField(5, embedding=zeta5^2)
            sage: K._log_gen(zeta5)
            3

            sage: K60.<zeta60> = CyclotomicField(60)
            sage: K30.<zeta30> = CyclotomicField(30, embedding=zeta60**2)
            sage: K15.<zeta15> = CyclotomicField(15, embedding=zeta30**2)
            sage: K5.<zeta5> = CyclotomicField(5, embedding=zeta15**12)
            sage: K60._log_gen(zeta30)
            2
            sage: K60._log_gen(zeta15)
            4
            sage: K60._log_gen(zeta5)
            48
            sage: K5._log_gen(zeta15**3)
            4
        """
        X = x.parent()
        gen = self.gen()

        if self.has_coerce_map_from(X):
            Y = self
            x = self(x)
        elif X.has_coerce_map_from(self):
            Y = X
            gen = X(self.gen())
        else:
            return

        n = self._n()
        if CDF.has_coerce_map_from(Y):
            x = CDF(x)
            gen = CDF(gen)
            # Let zeta = e^(2*pi*i/n)
            two_pi = 2*RDF.pi()
            a = (n * x.arg() / two_pi).round()        # x = zeta^a
            b = (n * gen.arg() / two_pi).round()      # gen = zeta^b
            e = mod(a/b, n).lift()          # e is the expected result
            if abs(gen**e-x) < 1/n:         # a sanity check
                return e
        else:
            # NOTE: this can be *very* slow!
            gen_pow_e = 1
            for e in range(n):
                if gen_pow_e == x:
                    return e
                gen_pow_e *= gen

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of this cyclotomic field from `x`.

        EXAMPLES:

        The following example illustrates coercion from the
        cyclotomic field `\QQ(\zeta_{42})` to the cyclotomic field `\QQ(\zeta_6)`, in
        a case where such coercion is defined::

            sage: k42 = CyclotomicField(42)
            sage: k6 = CyclotomicField(6)
            sage: a = k42.gen(0)
            sage: b = a^7; b
            zeta42^7
            sage: k6(b)  # indirect doctest
            zeta6
            sage: b^2
            zeta42^7 - 1
            sage: k6(b^2)
            zeta6 - 1

        Conversion of elements of the :class:`~sage.rings.universal_cyclotomic_field.UniversalCyclotomicField`::

            sage: CF = CyclotomicField(5)
            sage: UCF.<E> = UniversalCyclotomicField()                                  # needs sage.libs.gap
            sage: CF(E(5))                                                              # needs sage.libs.gap
            zeta5

            sage: CF = CyclotomicField(10)
            sage: CF(E(5))                                                              # needs sage.libs.gap
            zeta10^2

        Coercion of GAP cyclotomic elements is also supported::

            sage: CyclotomicField(18)(gap('E(3)'))  # indirect doctest                  # needs sage.libs.gap
            zeta18^3 - 1

        Converting from rings of integers::

            sage: K.<z> = CyclotomicField(7)
            sage: O = K.maximal_order()
            sage: K(O.1)
            z
            sage: K(O.1^2 + O.1 - 2)
            z^2 + z - 2
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            if isinstance(x.parent(), NumberField_cyclotomic):
                return self._coerce_from_other_cyclotomic_field(x)
            else:
                return NumberField_absolute._element_constructor_(self, x)
        elif isinstance(x, pari_gen):
            return NumberField_absolute._element_constructor_(self, x, check=check)
        elif isinstance(x, (LibGapElement, GapElement)):
            return self._coerce_from_gap(x)
        elif isinstance(x, str):
            return self._convert_from_str(x)
        elif isinstance(x, UniversalCyclotomicFieldElement):
            return x.to_cyclotomic_field(self)
        else:
            return self._convert_non_number_field_element(x)

    def _coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
        """
        Coerce an element `x` of a cyclotomic field into ``self``, if at all
        possible.

        INPUT:

        - ``x`` -- number field element

        - ``only_canonical`` -- boolean (default: ``False``); attempt
          to work, even in some cases when `x` is not in a subfield of the
          cyclotomics (as long as `x` is a root of unity)

        EXAMPLES::

            sage: K = CyclotomicField(24); L = CyclotomicField(48)
            sage: L._coerce_from_other_cyclotomic_field(K.0 + 1)
            zeta48^2 + 1
            sage: K(L.0**2)
            zeta24
        """
        K = x.parent()
        if K is self:
            return x
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
                for r in range(y.multiplicative_order()):
                    if z == x:
                        return self.zeta(m)**(r+1)
                    z *= y
            raise TypeError("cannot coerce %s into %s" % (x, self))
        return self._element_class(self, x)

    def _coerce_from_gap(self, x):
        """
        Attempt to coerce a GAP number field element into this cyclotomic
        field.

        EXAMPLES::

            sage: k5.<z> = CyclotomicField(5)
            sage: w = libgap.eval('E(5)^7 + 3'); w                                      # needs sage.libs.gap
            -3*E(5)-2*E(5)^2-3*E(5)^3-3*E(5)^4
            sage: z^7 + 3                                                               # needs sage.libs.gap
            z^2 + 3
            sage: k5(w)  # indirect doctest                                             # needs sage.libs.gap
            z^2 + 3

        It may be that GAP uses a name for the generator of the cyclotomic field.
        We can deal with this case, if this name coincides with the name in Sage::

            sage: F = CyclotomicField(8)
            sage: z = F.gen()
            sage: a = libgap(z + 1/z); a                                                # needs sage.libs.gap
            E(8)-E(8)^3
            sage: F(a)                                                                  # needs sage.libs.gap
            -zeta8^3 + zeta8

        Matrices over cyclotomic fields are correctly dealt with it as well::

            sage: b = libgap.eval('[[E(4), 1], [0, 1+E(8)-E(8)^3]]')                    # needs sage.libs.gap
            sage: matrix(F, b)                                                          # needs sage.libs.gap
            [             zeta8^2                    1]
            [                   0 -zeta8^3 + zeta8 + 1]

        It also works with the old pexpect interface to GAP::

            sage: # needs sage.libs.gap
            sage: a = gap(z + 1/z)
            sage: b = gap(Matrix(F,[[z^2,1],[0,a+1]])); b
            [ [ E(4), 1 ], [ 0, 1+E(8)-E(8)^3 ] ]
            sage: b[1,2]
            1
            sage: F(b[1,2])
            1
            sage: matrix(F, b)
            [             zeta8^2                    1]
            [                   0 -zeta8^3 + zeta8 + 1]
        """
        if x.IsRat():
            return self(QQ(x))
        coeffs = x.CoeffsCyc(self.__n)
        zeta = self.gen()
        return sum(QQ(c) * zeta**i for i, c in enumerate(coeffs))

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from the cyclotomic field ``self`` to
        the number field codomain.

        The ``cat`` option is currently ignored.

        EXAMPLES:

        This function is implicitly called by the Hom method or
        function.

        ::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 3); K
            Number Field in a with defining polynomial x^2 + 3
            sage: CyclotomicField(3).Hom(K)  # indirect doctest
            Set of field embeddings
             from Cyclotomic Field of order 3 and degree 2
               to Number Field in a with defining polynomial x^2 + 3
            sage: End(CyclotomicField(21))
            Automorphism group of Cyclotomic Field of order 21 and degree 12
        """
        if is_NumberFieldHomsetCodomain(codomain):
            from sage.rings.number_field.homset import CyclotomicFieldHomset
            return CyclotomicFieldHomset(self, codomain)
        else:
            raise TypeError

    def is_galois(self):
        """
        Return ``True`` since all cyclotomic fields are automatically Galois.

        EXAMPLES::

            sage: CyclotomicField(29).is_galois()
            True
        """
        return True

    def is_abelian(self):
        """
        Return ``True`` since all cyclotomic fields are automatically abelian.

        EXAMPLES::

            sage: CyclotomicField(29).is_abelian()
            True
        """
        return True

    def is_isomorphic(self, other):
        """
        Return ``True`` if the cyclotomic field ``self`` is isomorphic as a number
        field to ``other``.

        EXAMPLES::

            sage: CyclotomicField(11).is_isomorphic(CyclotomicField(22))
            True
            sage: CyclotomicField(11).is_isomorphic(CyclotomicField(23))
            False
            sage: x = polygen(QQ, 'x')
            sage: CyclotomicField(3).is_isomorphic(NumberField(x^2 + x + 1, 'a'))
            True
            sage: CyclotomicField(18).is_isomorphic(CyclotomicField(9))
            True
            sage: CyclotomicField(10).is_isomorphic(NumberField(x^4 - x^3 + x^2 - x + 1, 'b'))
            True

        Check :issue:`14300`::

            sage: K = CyclotomicField(4)
            sage: N = K.extension(x^2 - 5, 'z')
            sage: K.is_isomorphic(N)
            False
            sage: K.is_isomorphic(CyclotomicField(8))
            False
        """
        if isinstance(other, NumberField_cyclotomic):
            return self.zeta_order() == other.zeta_order()
        return NumberField_generic.is_isomorphic(self, other)

    def complex_embedding(self, prec=53):
        r"""
        Return the embedding of this cyclotomic field into the approximate
        complex field with precision ``prec`` obtained by sending the generator
        `\zeta` of ``self`` to exp(2\*pi\*i/n), where `n` is
        the multiplicative order of `\zeta`.

        EXAMPLES::

            sage: C = CyclotomicField(4)
            sage: C.complex_embedding()
            Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Field with 53 bits of precision
              Defn: zeta4 |--> 6.12323399573677e-17 + 1.00000000000000*I

        Note in the example above that the way zeta is computed (using sine
        and cosine in MPFR) means that only the ``prec`` bits of the number
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
        CC = sage.rings.complex_mpfr.ComplexField(prec)
        return self.hom([CC.zeta(self._n())], check=False)

    @cached_method
    def embeddings(self, K):
        r"""
        Compute all field embeddings of this field into the field `K`.

        INPUT:

        - ``K`` -- a field

        EXAMPLES::

            sage: CyclotomicField(5).embeddings(ComplexField(53))[1]
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Field with 53 bits of precision
              Defn: zeta5 |--> -0.809016994374947 + 0.587785252292473*I
            sage: CyclotomicField(5).embeddings(Qp(11, 4, print_mode='digits'))[1]      # needs sage.rings.padics
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   11-adic Field with capped relative precision 4
              Defn: zeta5 |--> ...1525
        """
        n = self._n()
        if K.characteristic() == 0:
            try:
                z = K.zeta(n)
            except ValueError:
                # No nth root of unity
                v = []
            except AttributeError:
                # zeta not defined
                return super().embeddings(K)
            else:
                X = (m for m in range(n) if gcd(m, n) == 1)
                v = [self.hom([z**i], check=False) for i in X]
        else:
            v = []
        return Sequence(v, cr=True, immutable=True,
                        check=False, universe=self.Hom(K))

    def complex_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the approximate
        complex field with precision ``prec``.

        If you want 53-bit double precision, which is faster but less
        reliable, then do ``self.embeddings(CDF)``.

        EXAMPLES::

            sage: CyclotomicField(5).complex_embeddings()
            [Ring morphism:
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
               Defn: zeta5 |--> 0.309016994374947 - 0.951056516295154*I]
        """
        CC = sage.rings.complex_mpfr.ComplexField(prec)
        return self.embeddings(CC)

    def real_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the approximate
        real field with precision ``prec``.

        Mostly, of course, there are no such embeddings.

        EXAMPLES::

            sage: len(CyclotomicField(4).real_embeddings())
            0
            sage: CyclotomicField(2).real_embeddings()
            [Ring morphism:
               From: Cyclotomic Field of order 2 and degree 1
               To:   Real Field with 53 bits of precision
               Defn: -1 |--> -1.00000000000000]
        """
        K = sage.rings.real_mpfr.RealField(prec)
        return self.embeddings(K)

    def signature(self):
        r"""
        Return `(r_1, r_2)`, where `r_1` and `r_2` are the number of real embeddings
        and pairs of complex embeddings of this cyclotomic field,
        respectively.

        Trivial since, apart from `\QQ`, cyclotomic fields are totally
        complex.

        EXAMPLES::

            sage: CyclotomicField(5).signature()
            (0, 2)
            sage: CyclotomicField(2).signature()
            (1, 0)
        """
        m = ZZ(self.degree())
        if m == 1:
            return (ZZ.one(), ZZ(0))
        else:
            return (ZZ(0), m // 2)

    def different(self):
        """
        Return the different ideal of the cyclotomic field ``self``.

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
        Return the discriminant of the ring of integers of the cyclotomic
        field ``self``, or if ``v`` is specified, the determinant of the trace
        pairing on the elements of the list ``v``.

        Uses the formula for the discriminant of a prime power cyclotomic
        field and Hilbert Theorem 88 on the discriminant of composita.

        INPUT:

        - ``v`` -- (optional) list of elements of this number field

        OUTPUT: integer if ``v`` is omitted, and Rational otherwise

        EXAMPLES::

            sage: CyclotomicField(20).discriminant()
            4000000
            sage: CyclotomicField(18).discriminant()
            -19683
        """
        if v is None:
            try:
                return self.__disc
            except AttributeError:
                n = self._n()
                deg = self.degree()
                d = ZZ.one()  # so that CyclotomicField(1).disc() has the right type
                factors = n.factor()
                for (p, r) in factors:
                    e = (r * p - r - 1) * deg // (p - 1)
                    d *= p**e
                sign = 1
                if len(factors) == 1 and (n == 4 or factors[0][0].mod(4) == 3):
                    sign = -1
                elif len(factors) == 2 and factors[0] == (2, 1) and factors[1][0].mod(4) == 3:
                    sign = -1
                self.__disc = sign * d
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
            p = next_prime(p)
            if p % n == 1:
                return p

    def _pari_integral_basis(self, v=None, important=True):
        """
        Internal function returning an integral basis of this number field in
        PARI format.

        This field is cyclotomic, so this is a trivial computation,
        since the power basis on the generator is an integral basis.
        Thus the ``v`` and ``important`` parameters are ignored.

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
            for n in range(self.degree()):
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
        r"""
        Return a dictionary that maps powers of `\zeta` to their order. This
        makes computing the orders of the elements of finite order in this
        field faster.

        EXAMPLES::

            sage: v = CyclotomicField(6)._multiplicative_order_table()
            sage: w = sorted(v.items()); w
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
                t[x.polynomial()] = n // gcd(m, n)  # multiplicative_order of (zeta_n)**m
                x *= zeta
                m += 1
            self.__multiplicative_order_table = t
            return t

    def zeta(self, n=None, all=False):
        """
        Return an element of multiplicative order `n` in this
        cyclotomic field.

        If there is no such element, raise a :exc:`ValueError`.

        INPUT:

        - ``n`` -- integer (default: ``None``, returns element of
          maximal order)

        - ``all`` -- boolean (default: ``False``); whether to return
          a list of all primitive `n`-th roots of unity

        OUTPUT: root of unity or list

        EXAMPLES::

            sage: k = CyclotomicField(4)
            sage: k.zeta()
            zeta4
            sage: k.zeta(2)
            -1
            sage: k.zeta().multiplicative_order()
            4

        ::

            sage: k = CyclotomicField(21)
            sage: k.zeta().multiplicative_order()
            42
            sage: k.zeta(21).multiplicative_order()
            21
            sage: k.zeta(7).multiplicative_order()
            7
            sage: k.zeta(6).multiplicative_order()
            6
            sage: k.zeta(84)
            Traceback (most recent call last):
            ...
            ValueError: 84 does not divide order of generator (42)

        ::

            sage: K.<a> = CyclotomicField(7)
            sage: K.zeta(all=True)
            [-a^4, -a^5, a^5 + a^4 + a^3 + a^2 + a + 1, -a, -a^2, -a^3]
            sage: K.zeta(14, all=True)
            [-a^4, -a^5, a^5 + a^4 + a^3 + a^2 + a + 1, -a, -a^2, -a^3]
            sage: K.zeta(2, all=True)
            [-1]
            sage: K.<a> = CyclotomicField(10)
            sage: K.zeta(20, all=True)
            Traceback (most recent call last):
            ...
            ValueError: 20 does not divide order of generator (10)

        ::

            sage: K.<a> = CyclotomicField(5)
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 does not divide order of generator (10)
            sage: v = K.zeta(5, all=True); v
            [a, a^2, a^3, -a^3 - a^2 - a - 1]
            sage: [b^5 for b in v]
            [1, 1, 1, 1]
        """
        if n is None:
            n = self.zeta_order()
        else:
            n = Integer(n)

        z = self.gen()
        m = self._n()
        if n % 2 == 0 and m % 2 == 1:
            # In the n-th cyclotomic field, n odd, there are
            # actually 2*n-th roots of unity, so we include them.
            z = -z**((m+1)//2)  # -z
            m = 2 * m
        if m % n != 0:
            raise ValueError("%s does not divide order of generator (%s)" %
                             (n, self.zeta_order()))
        a = z**(m // n)
        if not all:
            return a

        v = [a]
        b = a * a
        for i in range(2, n):
            if n.gcd(i).is_one():
                v.append(b)
            b = b * a
        return v

    def number_of_roots_of_unity(self):
        """
        Return number of roots of unity in this cyclotomic field.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(21)
            sage: K.number_of_roots_of_unity()
            42
        """
        n = self._n()
        if n % 2:
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
            sage: [z**K.number_of_roots_of_unity() for z in zs]
            [1, 1, 1, 1, 1, 1]
        """
        n = self._n()
        v = self.gen().powers(n)
        if n % 2:
            v += [-x for x in v]
        return v


class NumberField_quadratic(NumberField_absolute, sage.rings.abc.NumberField_quadratic):
    r"""
    Create a quadratic extension of the rational field.

    The command ``QuadraticField(a)`` creates the field `\QQ(\sqrt{a})`.

    EXAMPLES::

        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
        sage: QuadraticField(-4, 'b')
        Number Field in b with defining polynomial x^2 + 4 with b = 2*I
    """
    def __init__(self, polynomial, name=None, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None, structure=None):
        """
        Create a quadratic number field.

        EXAMPLES::

            sage: k.<a> = QuadraticField(5, check=False); k
            Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?

        Don't do this::

            sage: k.<a> = QuadraticField(4, check=False); k
            Number Field in a with defining polynomial x^2 - 4 with a = 2

        TESTS::

            sage: k.<a> = QuadraticField(7)
            sage: type(k.zero())
            <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic_sqrt'>
            sage: type(k.one())
            <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic_sqrt'>

            sage: TestSuite(k).run()

        Check that :issue:`23008` is fixed::

            sage: z = polygen(ZZ, 'z')
            sage: K.<phi> = NumberField(z^2 - z - 1, embedding=QQbar(golden_ratio))     # needs sage.symbolic
            sage: floor(phi)                                                            # needs sage.symbolic
            1
        """
        NumberField_absolute.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes, structure=structure)
        self._standard_embedding = True

        # set the generator and element class
        c, b, a = (QQ(t) for t in self.defining_polynomial().list())
        Dpoly = b*b - 4*a*c
        D = (Dpoly.numer() * Dpoly.denom()).squarefree_part(bound=10000)
        self._D = D
        parts = -b/(2*a), (Dpoly/D).sqrt()/(2*a)

        if a.is_one() and b.is_zero() and c.is_one():
            self._element_class = number_field_element_quadratic.NumberFieldElement_gaussian
        else:
            if number_field_element_quadratic.is_sqrt_disc(parts[0], parts[1]):
                self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic_sqrt
            else:
                self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic

        self._NumberField_generic__gen = self._element_class(self, parts)

        # we must set the flag _standard_embedding *before* any element creation
        # Note that in the following code, no element is built.
        if self.coerce_embedding() is not None and CDF.has_coerce_map_from(self):
            rootD = CDF(number_field_element_quadratic.NumberFieldElement_quadratic(self, (QQ(0), QQ(1))))
            if D > 0:
                self._standard_embedding = rootD.real() > 0
            else:
                self._standard_embedding = rootD.imag() > 0

        # we reset _NumberField_generic__gen has the flag standard_embedding
        # might be modified
        self._NumberField_generic__gen = self._element_class(self, parts)

        # NumberField_absolute.__init__(...) set _zero_element and
        # _one_element to NumberFieldElement_absolute values, which is
        # wrong (and dangerous; such elements can actually be used to
        # crash Sage: see #5316).  Overwrite them with correct values.
        self._zero_element = self._element_class(self, (QQ(0), QQ(0)))
        self._one_element = self._element_class(self, (QQ(1), QQ(0)))

    def _coerce_map_from_(self, K):
        """
        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: f = K.coerce_map_from(QQ); f  # indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3 with a = 1.732050807568878?*I
            sage: f(3/5)
            3/5
            sage: parent(f(3/5)) is K
            True

            sage: g = K.coerce_map_from(ZZ); g  # indirect doctest
            Natural morphism:
              From: Integer Ring
              To:   Number Field in a with defining polynomial x^2 + 3 with a = 1.732050807568878?*I
            sage: g(1)
            1
            sage: parent(g(1)) is K
            True
        """
        if K is ZZ:
            return number_field_element_quadratic.Z_to_quadratic_field_element(self)
        if K is int:
            return self._coerce_map_via([ZZ], int)  # faster than direct
        if K is QQ:
            return number_field_element_quadratic.Q_to_quadratic_field_element(self)
        return NumberField_absolute._coerce_map_from_(self, K)

    def _latex_(self):
        r"""
        Return the latex representation of this quadratic field.

        EXAMPLES::

            sage: Z = QuadraticField(7)
            sage: latex(Z)  # indirect doctest
            \Bold{Q}(\sqrt{7})

            sage: Z = QuadraticField(7, latex_name='x')
            sage: latex(Z)  # indirect doctest
            \Bold{Q}[x]/(x^{2} - 7)
        """
        v = self.latex_variable_names()[0]
        if v.startswith('\\sqrt'):
            return "%s(%s)" % (latex(QQ), v)
        else:
            return NumberField_generic._latex_(self)

    def _polymake_init_(self):
        r"""
        Return the polymake representation of this quadratic field.

        This is merely a string, and does not represent a specific quadratic field.
        In polymake, only the elements know which field they belong to.

        EXAMPLES::

            sage: Z = QuadraticField(7)
            sage: polymake(Z)    # optional - jupymake  # indirect doctest
            QuadraticExtension
        """
        return '"QuadraticExtension"'

    def discriminant(self, v=None):
        """
        Return the discriminant of the ring of integers of the number
        field, or if ``v`` is specified, the determinant of the trace pairing
        on the elements of the list ``v``.

        INPUT:

        - ``v`` -- (optional) list of element of this number field

        OUTPUT: integer if ``v`` is omitted, and Rational otherwise

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.discriminant()
            -4
            sage: K.<a> = NumberField(x^2 + 5)
            sage: K.discriminant()
            -20
            sage: K.<a> = NumberField(x^2 - 5)
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

    def is_galois(self):
        """
        Return ``True`` since all quadratic fields are automatically Galois.

        EXAMPLES::

            sage: QuadraticField(1234,'d').is_galois()
            True
        """
        return True

    def class_number(self, proof=None):
        r"""
        Return the size of the class group of ``self``.

        INPUT:

        - ``proof`` -- boolean (default: ``True``, unless you called
          :meth:`proof.number_field` and set it otherwise).  If
          ``proof`` is ``False`` (*not* the default!), and the
          discriminant of the field is negative, then the following
          warning from the PARI manual applies:

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

            sage: [d for d in prime_range(2,300)
            ....:  if not is_square(d) and QuadraticField(-d,'a').class_number() == 1]
            [2, 3, 7, 11, 19, 43, 67, 163]

        It is an open problem to *prove* that there are infinity many
        positive square-free `d` such that
        `\QQ(\sqrt{d})` has class number `1`:

        ::

            sage: len([d for d in range(2,200)
            ....:      if not is_square(d) and QuadraticField(d,'a').class_number() == 1])
            121

        TESTS::

            sage: type(QuadraticField(-23,'a').class_number())
            <class 'sage.rings.integer.Integer'>
            sage: x = polygen(QQ, 'x')
            sage: type(NumberField(x^3 + 23, 'a').class_number())
            <class 'sage.rings.integer.Integer'>
            sage: type(NumberField(x^3 + 23, 'a').extension(x^2 + 5, 'b').class_number())
            <class 'sage.rings.integer.Integer'>
            sage: type(CyclotomicField(10).class_number())
            <class 'sage.rings.integer.Integer'>
        """
        proof = proof_flag(proof)
        try:
            return self.__class_number
        except AttributeError:
            self.__class_number = self.discriminant().class_number(proof)
            return self.__class_number

    def hilbert_class_field_defining_polynomial(self, name='x'):
        r"""
        Return a polynomial over `\QQ` whose roots generate the
        Hilbert class field of this quadratic field as an extension of
        this quadratic field.

        .. NOTE::

            Computed using PARI via Schertz's method. This
            implementation is quite fast.

        EXAMPLES::

            sage: K.<b> = QuadraticField(-23)
            sage: K.hilbert_class_field_defining_polynomial()
            x^3 - x^2 + 1

        Note that this polynomial is not the actual Hilbert class
        polynomial: see ``hilbert_class_polynomial``::

            sage: K.hilbert_class_polynomial()                                          # needs sage.schemes
            x^3 + 3491750*x^2 - 5151296875*x + 12771880859375

        ::

            sage: K.<a> = QuadraticField(-431)
            sage: K.class_number()
            21
            sage: K.hilbert_class_field_defining_polynomial(name='z')
            z^21 + 6*z^20 + 9*z^19 - 4*z^18 + 33*z^17 + 140*z^16 + 220*z^15 + 243*z^14
             + 297*z^13 + 461*z^12 + 658*z^11 + 743*z^10 + 722*z^9 + 681*z^8 + 619*z^7
             + 522*z^6 + 405*z^5 + 261*z^4 + 119*z^3 + 35*z^2 + 7*z + 1
        """
        f = pari(self.discriminant()).quadhilbert()
        return QQ[name](f)

    def hilbert_class_field(self, names):
        r"""
        Return the Hilbert class field of this quadratic field as a
        relative extension of this field.

        .. NOTE::

            For the polynomial that defines this field as a relative
            extension, see the method :meth:`hilbert_class_field_defining_polynomial`,
            which is vastly faster than this method, since it doesn't
            construct a relative extension.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2 + 23)
            sage: L = K.hilbert_class_field('b'); L
            Number Field in b with defining polynomial x^3 - x^2 + 1 over its base field
            sage: L.absolute_field('c')
            Number Field in c with defining polynomial
             x^6 - 2*x^5 + 70*x^4 - 90*x^3 + 1631*x^2 - 1196*x + 12743
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
            sage: K.hilbert_class_polynomial()                                          # needs sage.schemes
            x

            sage: K.<a> = QuadraticField(-31)
            sage: K.hilbert_class_polynomial(name='z')                                  # needs sage.schemes
            z^3 + 39491307*z^2 - 58682638134*z + 1566028350940383
        """
        D = self.discriminant()

        if D > 0:
            raise NotImplementedError("Hilbert class polynomial is not implemented for real quadratic fields.")

        from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial as HCP
        return QQ[name](HCP(D))

    def number_of_roots_of_unity(self):
        """
        Return the number of roots of unity in this quadratic field.

        This is always 2 except when `d` is `-3` or `-4`.

        EXAMPLES::

            sage: QF = QuadraticField
            sage: [QF(d).number_of_roots_of_unity() for d in range(-7, -2)]
            [2, 2, 2, 4, 6]
        """
        d = self.discriminant()
        if d == -4:
            return 4
        if d == -3:
            return 6
        return 2

    def order_of_conductor(self, f):
        r"""
        Return the unique order with the given conductor in this quadratic field.

        .. SEEALSO ::

            :meth:`sage.rings.number_field.order.Order.conductor`

        EXAMPLES::

            sage: K.<t> = QuadraticField(-123)
            sage: K.order_of_conductor(1) is K.maximal_order()
            True
            sage: K.order_of_conductor(2).gens()
            (1, t)
            sage: K.order_of_conductor(44).gens()
            (1, 22*t)
            sage: K.order_of_conductor(9001).conductor()
            9001
        """
        f = ZZ(f)
        if f <= 0:
            raise ValueError('conductor must be a positive integer')
        return self.order([f * g for g in self.maximal_order().gens()])


def is_fundamental_discriminant(D):
    r"""
    Return ``True`` if the integer `D` is a fundamental
    discriminant, i.e., if `D \cong 0,1\pmod{4}`, and
    `D\neq 0, 1` and either (1) `D` is square free or
    (2) we have `D\cong 0\pmod{4}` with
    `D/4 \cong 2,3\pmod{4}` and `D/4` square free. These
    are exactly the discriminants of quadratic fields.

    EXAMPLES::

        sage: [D for D in range(-15,15) if is_fundamental_discriminant(D)]
        ...
        DeprecationWarning: is_fundamental_discriminant(D) is deprecated;
        please use D.is_fundamental_discriminant()
        ...
        [-15, -11, -8, -7, -4, -3, 5, 8, 12, 13]
        sage: [D for D in range(-15,15)
        ....:  if not is_square(D) and QuadraticField(D,'a').disc() == D]
        [-15, -11, -8, -7, -4, -3, 5, 8, 12, 13]
    """
    deprecation(35147, "is_fundamental_discriminant(D) is deprecated; please use D.is_fundamental_discriminant()")
    return Integer(D).is_fundamental_discriminant()


###################
# For pickling
###################


def NumberField_absolute_v1(poly, name, latex_name, canonical_embedding=None):
    """
    Used for unpickling old pickles.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import NumberField_absolute_v1
        sage: R.<x> = QQ[]
        sage: NumberField_absolute_v1(x^2 + 1, 'i', 'i')
        Number Field in i with defining polynomial x^2 + 1
    """
    return NumberField(polynomial=poly, name=name, latex_name=latex_name, check=False, embedding=canonical_embedding)


NumberField_generic_v1 = NumberField_absolute_v1  # for historical reasons only (so old objects unpickle)


def NumberField_cyclotomic_v1(zeta_order, name, canonical_embedding=None):
    """
    Used for unpickling old pickles.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import NumberField_cyclotomic_v1
        sage: NumberField_cyclotomic_v1(5,'a')
        Cyclotomic Field of order 5 and degree 4
        sage: NumberField_cyclotomic_v1(5,'a').variable_name()
        'a'
    """
    return CyclotomicField(n=zeta_order, names=name, embedding=canonical_embedding)


def NumberField_quadratic_v1(poly, name, canonical_embedding=None):
    """
    Used for unpickling old pickles.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import NumberField_quadratic_v1
        sage: R.<x> = QQ[]
        sage: NumberField_quadratic_v1(x^2 - 2, 'd')
        Number Field in d with defining polynomial x^2 - 2
    """
    return NumberField(polynomial=poly, name=name, check=False, embedding=canonical_embedding)


def put_natural_embedding_first(v):
    """
    Helper function for embeddings() functions for number fields.

    INPUT:

    - ``v`` -- list of embeddings of a number field

    OUTPUT: none; the list is altered in-place, so that, if possible, the first
    embedding has been switched with one of the others, so that if there is an
    embedding which preserves the generator names then it appears first.

    EXAMPLES::

        sage: K.<a> = CyclotomicField(7)
        sage: embs = K.embeddings(K)
        sage: [e(a) for e in embs]  # random - there is no natural sort order
        [a, a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1]
        sage: id = [e for e in embs if e(a) == a][0]; id
        Ring endomorphism of Cyclotomic Field of order 7 and degree 6
          Defn: a |--> a
        sage: permuted_embs = list(embs); permuted_embs.remove(id); permuted_embs.append(id)
        sage: [e(a) for e in permuted_embs]  # random - but natural map is not first
        [a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1, a]
        sage: permuted_embs[0] != a
        True
        sage: from sage.rings.number_field.number_field import put_natural_embedding_first
        sage: put_natural_embedding_first(permuted_embs)
        sage: [e(a) for e in permuted_embs]  # random - but natural map is first
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
    `\CC`, return an equivalent embedding with higher precision.

    INPUT:

    - ``e`` -- an embedding of a number field into either `\RR` or `\CC`
      (with some precision)

    - ``prec`` -- (default: ``None``) the desired precision; if ``None``,
      current precision is doubled; if ``Infinity``, the equivalent
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

        sage: x = polygen(QQ, 'x')
        sage: K.<a> = NumberField(x^3 - 2)
        sage: K.signature()
        (1, 1)
        sage: e = K.embeddings(RR)[0]; e
        Ring morphism:
        From: Number Field in a with defining polynomial x^3 - 2
        To:   Real Field with 53 bits of precision
        Defn: a |--> 1.25992104989487
        sage: e = refine_embedding(e, Infinity); e
        Ring morphism:
        From: Number Field in a with defining polynomial x^3 - 2
        To:   Algebraic Real Field
        Defn: a |--> 1.259921049894873?

    Now we can obtain arbitrary precision values with no trouble::

        sage: RealField(150)(e(a))
        1.2599210498948731647672106072782283505702515
        sage: _^3
        2.0000000000000000000000000000000000000000000
        sage: RealField(200)(e(a^2 - 3*a + 7))
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
        -0.62996052494743658238360530363911417528512573235075399004099
         - 1.0911236359717214035600726141898088813258733387403009407036*I
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
          Defn: zeta7 |--> 0.623489801858733530525004884004239810632274730896402105365549439096853652456487284575942507
                            + 0.781831482468029808708444526674057750232334518708687528980634958045091731633936441700868007*I
        sage: refine_embedding(x, infinity)
        Ring morphism:
          From: Cyclotomic Field of order 7 and degree 6
          To:   Algebraic Field
          Defn: zeta7 |--> 0.6234898018587335? + 0.7818314824680299?*I

    When the old embedding is into the real lazy field,
    then only real embeddings should be considered.
    See :issue:`17495`::

        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3 + x - 1, embedding=0.68)
        sage: from sage.rings.number_field.number_field import refine_embedding
        sage: refine_embedding(K.specified_complex_embedding(), 100)
        Ring morphism:
          From: Number Field in a with defining polynomial x^3 + x - 1 with a = 0.6823278038280193?
          To:   Real Field with 100 bits of precision
          Defn: a |--> 0.68232780382801932736948373971
        sage: refine_embedding(K.specified_complex_embedding(), Infinity)
        Ring morphism:
          From: Number Field in a with defining polynomial x^3 + x - 1 with a = 0.6823278038280193?
          To:   Algebraic Real Field
          Defn: a |--> 0.6823278038280193?
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
    if isinstance(RC, (sage.rings.abc.RealField, sage.rings.abc.RealDoubleField)) or RC == RLF:
        if prec == Infinity:
            elist = K.embeddings(sage.rings.qqbar.AA)
        else:
            elist = K.real_embeddings(prec)
    else:
        if prec == Infinity:
            elist = K.embeddings(sage.rings.qqbar.QQbar)
        else:
            elist = K.complex_embeddings(prec)

    # Now we determine which is an extension of the old one; this
    # relies on the fact that coercing a high-precision root into a
    # field with lower precision will equal the lower-precision root!
    diffs = [(RC(ee(K.gen()))-old_root).abs() for ee in elist]
    return elist[min(zip(diffs, count()))[1]]


def is_real_place(v):
    r"""
    Return ``True`` if `v` is real, ``False`` if `v` is complex.

    INPUT:

    - ``v`` -- an infinite place of ``self``

    OUTPUT:

    A boolean indicating whether a place is real (``True``) or complex (``False``).

    EXAMPLES::

        sage: x = polygen(QQ, 'x')
        sage: K.<xi> = NumberField(x^3 - 3)
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: is_real_place(phi_real)
        True

        sage: is_real_place(phi_complex)
        False

    It is an error to put in a finite place

    ::

        sage: is_real_place(v_fin)
        Traceback (most recent call last):
        ...
        AttributeError: 'NumberFieldFractionalIdeal' object has no attribute 'im_gens'...
    """
    RR = sage.rings.real_mpfr.RealField(53)
    try:
        RR(v.im_gens()[0])
        return True
    except TypeError:
        return False


def _splitting_classes_gens_(K, m, d):
    r"""
    Given a number field `K` of conductor `m` and degree `d`,
    this returns a set of multiplicative generators of the
    subgroup of `(\ZZ/m\ZZ)^{\times}`
    containing exactly the classes that contain the primes splitting
    completely in `K`.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import _splitting_classes_gens_
        sage: K = CyclotomicField(101)
        sage: L = K.subfields(20)[0][0]
        sage: L.conductor()                                                             # needs sage.groups
        101
        sage: _splitting_classes_gens_(L,101,20)        # optional - gap_package_polycyclic, needs sage.libs.gap
        [95]

        sage: K = CyclotomicField(44)
        sage: L = K.subfields(4)[0][0]
        sage: _splitting_classes_gens_(L,44,4)  # optional - gap_package_polycyclic, needs sage.libs.gap
        [37]

        sage: K = CyclotomicField(44)
        sage: L = K.subfields(5)[0][0]
        sage: K.degree()
        20
        sage: L
        Number Field in zeta44_0 with defining polynomial x^5 - 2*x^4 - 16*x^3 + 24*x^2 + 48*x - 32
         with zeta44_0 = 3.837971894457990?
        sage: L.conductor()                                                             # needs sage.groups
        11
        sage: _splitting_classes_gens_(L,11,5)  # optional - gap_package_polycyclic, needs sage.libs.gap
        [10]
    """
    from sage.groups.abelian_gps.abelian_group import AbelianGroup

    Zm = IntegerModRing(m)
    unit_gens = Zm.unit_gens()
    Zmstar = AbelianGroup(len(unit_gens), [x.multiplicative_order() for x in unit_gens])

    def map_Zmstar_to_Zm(h):
        li = h.list()
        return prod(unit_gens[i]**li[i] for i in range(len(unit_gens)))

    Hgens = []
    H = Zmstar.subgroup([])
    p = 0
    Horder = euler_phi(m) / d
    for g in Zmstar:
        if H.order() == Horder:
            break
        if g not in H:
            u = map_Zmstar_to_Zm(g)
            p = u.lift()
            while not p.is_prime():
                p += m
            f = K.fractional_ideal(p).prime_factors()[0].residue_class_degree()
            h = g**f
            if h not in H:
                Hgens += [h]
                H = Zmstar.subgroup(Hgens)

    return [map_Zmstar_to_Zm(h) for h in Hgens]
