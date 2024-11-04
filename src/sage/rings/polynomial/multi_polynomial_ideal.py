# sage.doctest: needs sage.libs.singular
r"""
Ideals in multivariate polynomial rings

Sage has a powerful system to compute with multivariate polynomial
rings. Most algorithms dealing with these ideals are centered on the
computation of *Groebner bases*. Sage mainly uses Singular to
implement this functionality. Singular is widely regarded as the best
open-source system for Groebner basis calculation in multivariate
polynomial rings over fields.

EXAMPLES:

We compute a Groebner basis for some given ideal. The type returned by
the ``groebner_basis`` method is ``PolynomialSequence``, i.e., it is not a
:class:`MPolynomialIdeal`::

    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: B = I.groebner_basis()
    sage: type(B)
    <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>

Groebner bases can be used to solve the ideal membership problem::

    sage: f,g,h = B
    sage: (2*x*f + g).reduce(B)
    0

    sage: (2*x*f + g) in I
    True

    sage: (2*x*f + 2*z*h + y^3).reduce(B)
    y^3

    sage: (2*x*f + 2*z*h + y^3) in I
    False

We compute a Groebner basis for Cyclic 6, which is a standard
benchmark and test ideal. ::

    sage: R.<x,y,z,t,u,v> = QQ['x,y,z,t,u,v']
    sage: I = sage.rings.ideal.Cyclic(R,6)
    sage: B = I.groebner_basis()
    sage: len(B)
    45

We compute in a quotient of a polynomial ring over `\ZZ/17\ZZ`::

    sage: R.<x,y> = ZZ[]
    sage: S.<a,b> = R.quotient((x^2 + y^2, 17))
    sage: S
    Quotient of Multivariate Polynomial Ring in x, y over Integer Ring
    by the ideal (x^2 + y^2, 17)

    sage: a^2 + b^2 == 0
    True
    sage: a^3 - b^2
    -a*b^2 - b^2

Note that the result of a computation is not necessarily reduced::

    sage: (a+b)^17
    a*b^16 + b^17
    sage: S(17) == 0
    True

Or we can work with `\ZZ/17\ZZ` directly::

    sage: R.<x,y> = Zmod(17)[]
    sage: S.<a,b> = R.quotient((x^2 + y^2,))
    sage: S
    Quotient of Multivariate Polynomial Ring in x, y over Ring of
    integers modulo 17 by the ideal (x^2 + y^2)

    sage: a^2 + b^2 == 0
    True
    sage: a^3 - b^2 == -a*b^2 - b^2 == 16*a*b^2 + 16*b^2
    True
    sage: (a+b)^17
    a*b^16 + b^17
    sage: S(17) == 0
    True


Working with a polynomial ring over `\ZZ`::

    sage: R.<x,y,z,w> = ZZ[]
    sage: I = ideal(x^2 + y^2 - z^2 - w^2, x-y)
    sage: J = I^2
    sage: J.groebner_basis()
    [4*y^4 - 4*y^2*z^2 + z^4 - 4*y^2*w^2 + 2*z^2*w^2 + w^4,
     2*x*y^2 - 2*y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2,
     x^2 - 2*x*y + y^2]

    sage: y^2 - 2*x*y + x^2 in J
    True
    sage: 0 in J
    True

We do a Groebner basis computation over a number field::

    sage: K.<zeta> = CyclotomicField(3)
    sage: R.<x,y,z> = K[]; R
    Multivariate Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2

    sage: i = ideal(x - zeta*y + 1, x^3 - zeta*y^3); i
    Ideal (x + (-zeta)*y + 1, x^3 + (-zeta)*y^3) of Multivariate
    Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2

    sage: i.groebner_basis()
    [y^3 + (2*zeta + 1)*y^2 + (zeta - 1)*y + (-1/3*zeta - 2/3), x + (-zeta)*y + 1]

    sage: S = R.quotient(i); S
    Quotient of Multivariate Polynomial Ring in x, y, z over
    Cyclotomic Field of order 3 and degree 2 by the ideal (x +
    (-zeta)*y + 1, x^3 + (-zeta)*y^3)

    sage: S.0  - zeta*S.1
    -1
    sage: S.0^3 - zeta*S.1^3
    0

Two examples from the Mathematica documentation (done in Sage):

    We compute a Groebner basis::

        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: ideal(x^2 - 2*y^2, x*y - 3).groebner_basis()
        [x - 2/3*y^3, y^4 - 9/2]

    We show that three polynomials have no common root::

        sage: R.<x,y> = QQ[]
        sage: ideal(x+y, x^2 - 1, y^2 - 2*x).groebner_basis()
        [1]

The next example shows how we can use Groebner bases over `\ZZ` to
find the primes modulo which a system of equations has a solution,
when the system has no solutions over the rationals.

    We first form a certain ideal `I` in `\ZZ[x, y, z]`, and note that
    the Groebner basis of `I` over `\QQ` contains 1, so there are no
    solutions over `\QQ` or an algebraic closure of it (this is not
    surprising as there are 4 equations in 3 unknowns). ::

        sage: P.<x,y,z> = PolynomialRing(ZZ,order='lex')
        sage: I = ideal(-y^2 - 3*y + z^2 + 3, -2*y*z + z^2 + 2*z + 1,
        ....:           x*z + y*z + z^2, -3*x*y + 2*y*z + 6*z^2)
        sage: I.change_ring(P.change_ring(QQ)).groebner_basis()
        [1]

    However, when we compute the Groebner basis of `I` (defined over
    `\ZZ`), we note that there is a certain integer in the ideal
    which is not 1. ::

        sage: I.groebner_basis()
        [x + y + 57119*z + 4, y^2 + 3*y + 17220, y*z + ...,
         2*y + 158864, z^2 + 17223, 2*z + 41856, 164878]

    Now for each prime `p` dividing this integer 164878, the Groebner
    basis of `I` modulo `p` will be non-trivial and will thus give a
    solution of the original system modulo `p`. ::


        sage: factor(164878)
        2 * 7 * 11777

        sage: I.change_ring(P.change_ring(GF(2))).groebner_basis()                      # needs sage.rings.finite_rings
        [x + y + z, y^2 + y, y*z + y, z^2 + 1]
        sage: I.change_ring(P.change_ring(GF(7))).groebner_basis()                      # needs sage.rings.finite_rings
        [x - 1, y + 3, z - 2]
        sage: I.change_ring(P.change_ring(GF(11777))).groebner_basis()                  # needs sage.rings.finite_rings
        [x + 5633, y - 3007, z - 2626]

    The Groebner basis modulo any product of the prime factors is also non-trivial::

        sage: I.change_ring(P.change_ring(IntegerModRing(2 * 7))).groebner_basis()
        [x + ..., y^2 + 3*y, y*z + 7*y + 6, 2*y + 6, z^2 + 3, 2*z + 10]

    Modulo any other prime the Groebner basis is trivial so there are
    no other solutions. For example::

        sage: I.change_ring(P.change_ring(GF(3))).groebner_basis()                      # needs sage.rings.finite_rings
        [1]

TESTS::

    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: I == loads(dumps(I))
    True

.. NOTE::

    Sage distinguishes between lists or sequences of polynomials and
    ideals. Thus an ideal is not identified with a particular set of
    generators. For sequences of multivariate polynomials see
    :class:`sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic`.

AUTHORS:

- William Stein: initial version

- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of some Singular
  features

- Martin Albrecht (2007,2008): refactoring, many Singular related functions,
  added plot()

- Martin Albrecht (2009): added Groebner basis over rings functionality from
  Singular 3.1

- John Perry (2012): bug fixing equality & containment of ideals
"""

# ****************************************************************************
#
#                               Sage
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2008,2009 Martin Albrecht <malb@informatik.uni-bremen.de>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from warnings import warn

import sage.rings.abc
import sage.rings.polynomial.toy_buchberger as toy_buchberger
import sage.rings.polynomial.toy_d_basis as toy_d_basis
import sage.rings.polynomial.toy_variety as toy_variety
from sage.misc.cachefunc import cached_method
from sage.misc.method_decorator import MethodDecorator
from sage.misc.misc_c import prod
from sage.misc.verbose import get_verbose, verbose
from sage.rings.ideal import Ideal_generic
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.noncommutative_ideals import Ideal_nc
from sage.rings.qqbar_decorators import handle_AA_and_QQbar
from sage.structure.richcmp import (op_EQ, op_GE, op_GT, op_LE, op_LT, op_NE,
                                    rich_to_bool, richcmp_method)
from sage.structure.sequence import Sequence

try:
    from sage.interfaces.expect import StdOutContext
    from sage.interfaces.singular import singular as singular_default
    from sage.interfaces.singular import singular_gb_standard_options
    from sage.libs.singular.standard_options import \
        libsingular_gb_standard_options
except ImportError:
    singular_gb_standard_options = libsingular_gb_standard_options = MethodDecorator

try:
    from sage.interfaces.magma import magma as magma_default
    from sage.interfaces.magma import magma_gb_standard_options
except ImportError:
    magma_default = None
    magma_gb_standard_options = MethodDecorator


class RequireField(MethodDecorator):
    """
    Decorator which throws an exception if a computation over a
    coefficient ring which is not a field is attempted.

    .. NOTE::

        This decorator is used automatically internally so the user
        does not need to use it manually.
    """
    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(ZZ)
            sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
            sage: from sage.rings.polynomial.multi_polynomial_ideal import RequireField
            sage: class Foo(I.__class__):
            ....:     @RequireField
            ....:     def bar(I):
            ....:         return I.groebner_basis()
            ....:
            sage: J = Foo(I.ring(), I.gens())
            sage: J.bar()
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring must be a field for function 'bar'.
        """
        R = self._instance.ring()
        if not R.base_ring().is_field():
            raise ValueError("Coefficient ring must be a field for function '%s'." % (self.f.__name__))
        return self.f(self._instance, *args, **kwds)


require_field = RequireField


def is_MPolynomialIdeal(x) -> bool:
    """
    Return ``True`` if the provided argument ``x`` is an ideal in a
    multivariate polynomial ring.

    INPUT:

    - ``x`` -- an arbitrary object

    EXAMPLES::

        sage: from sage.rings.polynomial.multi_polynomial_ideal import is_MPolynomialIdeal
        sage: P.<x,y,z> = PolynomialRing(QQ)
        sage: I = [x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y + 2*y*z - y]

    Sage distinguishes between a list of generators for an ideal and
    the ideal itself. This distinction is inconsistent with Singular
    but matches Magma's behavior. ::

        sage: is_MPolynomialIdeal(I)
        doctest:warning...
        DeprecationWarning: The function is_MPolynomialIdeal is deprecated;
        use 'isinstance(..., MPolynomialIdeal)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        False

    ::

        sage: I = Ideal(I)
        sage: is_MPolynomialIdeal(I)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38266,
                "The function is_MPolynomialIdeal is deprecated; "
                "use 'isinstance(..., MPolynomialIdeal)' instead.")
    return isinstance(x, MPolynomialIdeal)


class MPolynomialIdeal_magma_repr:
    def _magma_init_(self, magma):
        """
        Return a Magma ideal matching this ideal if the base ring
        coerces to Magma and Magma is available.

        INPUT:

        - ``magma`` -- Magma instance

        EXAMPLES::

            sage: # optional - magma
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127), 10)
            sage: I = sage.rings.ideal.Cyclic(R,4)  # indirect doctest
            sage: magma(I)
            Ideal of Polynomial ring of rank 10 over GF(127)
            Order: Graded Reverse Lexicographical
            Variables: a, b, c, d, e, f, g, h, i, j
            Basis:
            [
            a + b + c + d,
            a*b + b*c + a*d + c*d,
            a*b*c + a*b*d + a*c*d + b*c*d,
            a*b*c*d + 126
            ]
        """
        P = magma(self.ring())
        G = magma(self.gens())
        return 'ideal<%s|%s>' % (P.name(), G._ref())

    @magma_gb_standard_options
    def _groebner_basis_magma(self, deg_bound=None, prot=False, magma=magma_default):
        """
        Compute a Groebner Basis for this ideal using Magma if
        available.

        INPUT:

        - ``deg_bound`` -- only compute to degree ``deg_bound``, that
          is, ignore all S-polynomials of higher degree (default: ``None``)

        - ``prot`` -- if ``True`` Magma's protocol is printed to stdout

        - ``magma`` -- Magma instance or ``None`` (default instance) (default: ``None``)

        EXAMPLES::

            sage: # optional - magma
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127), 10)
            sage: I = sage.rings.ideal.Cyclic(R, 6)
            sage: gb = I.groebner_basis('magma:GroebnerBasis')
            sage: len(gb)
            45

         We may also pass a degree bound to Magma::

            sage: # optional - magma
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127), 10)
            sage: I = sage.rings.ideal.Cyclic(R, 6)
            sage: gb = I.groebner_basis('magma:GroebnerBasis', deg_bound=4)
            sage: len(gb)
            5
        """
        R = self.ring()
        if not deg_bound:
            mself = magma(self)
        else:
            mself = magma(list(self.gens()))  # PolynomialSequence converts to a Magma Ideal too, so we force a list

        if get_verbose() >= 2:
            prot = True

        from sage.interfaces.magma import MagmaGBLogPrettyPrinter

        if prot:
            log_parser = MagmaGBLogPrettyPrinter(verbosity=get_verbose() + 1, style="sage" if prot == "sage" else "magma")
        else:
            log_parser = None

        ctx = StdOutContext(magma, silent=not prot, stdout=log_parser)
        if prot:
            magma.SetVerbose('Groebner', 1)
        with ctx:
            if deg_bound:
                mgb = mself.GroebnerBasis(deg_bound)
            else:
                mgb = mself.GroebnerBasis()

        if prot == "sage":
            print("")
            print("Highest degree reached during computation: %2d." % log_parser.max_deg)

        # TODO: rewrite this to be much more sophisticated in multi-level nested cases.
        mgb = [str(mgb[i+1]) for i in range(len(mgb))]
        if R.base_ring().degree() > 1:
            a = str(R.base_ring().gen())
            mgb = [e.replace("$.1",a) for e in mgb]

        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence

        B = PolynomialSequence([R(e) for e in mgb], R, immutable=True)
        return B


class MPolynomialIdeal_singular_base_repr:
    @require_field
    def syzygy_module(self):
        r"""
        Compute the first syzygy (i.e., the module of relations of the
        given generators) of the ideal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 2*x^2 + y
            sage: g = y
            sage: h = 2*f + g
            sage: I = Ideal([f,g,h])
            sage: M = I.syzygy_module(); M
            [       -2        -1         1]
            [       -y 2*x^2 + y         0]
            sage: G = vector(I.gens())
            sage: M*G
            (0, 0)

        ALGORITHM: Uses Singular's syz command
        """
        from sage.libs.singular.function_factory import ff
        syz = ff.syz
        from sage.matrix.constructor import matrix

        #return self._singular_().syz().transpose().sage_matrix(self.ring())
        S = syz(self)
        return matrix(self.ring(), S)

    @libsingular_gb_standard_options
    def _groebner_basis_libsingular(self, algorithm='groebner', *args, **kwds):
        """
        Return the reduced Groebner basis of this ideal.

        If the Groebner basis for this ideal has been calculated
        before, the cached Groebner basis is returned regardless
        of the requested algorithm.

        INPUT:

        - ``algorithm`` -- see below for available algorithms
        - ``redsb`` -- boolean (default: ``True``); return a reduced Groebner basis
        - ``red_tail`` -- boolean (default: ``True``); perform tail reduction

        ALGORITHMS:

        - ``'groebner'`` -- Singular's heuristic script (default)
        - ``'std'`` -- Buchberger's algorithm
        - ``'slimgb'`` -- the *SlimGB* algorithm
        - ``'stdhilb'`` -- Hilbert Basis driven Groebner basis
        - ``'stdfglm'`` -- Buchberger and FGLM

        EXAMPLES:

        We compute a Groebner basis of 'cyclic 4' relative to
        lexicographic ordering. ::

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4); I
            Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d
            + a*c*d + b*c*d, a*b*c*d - 1) of Multivariate Polynomial
            Ring in a, b, c, d over Rational Field

        ::

            sage: I._groebner_basis_libsingular()
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d,
            b*d^4 - b + d^5 - d, b*c - b*d + c^2*d^4 + c*d - 2*d^2,
            b^2 + 2*b*d + d^2, a + b + c + d]

        TESTS:

        We check that :issue:`17676` is fixed::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 3, order='lex')
            sage: I = Ideal(13*x*y*z+6*x*y+78*x*z+36*x-11*y^2*z-66*y*z,
            ....:           168*x*y*z+84*x*y+1008*x*z+504*x+12*y^3
            ....:            -154*y^2*z+72*y^2-924*y*z,
            ....:           -168*x^2*y*z^2-84*x^2*y*z-1008*x^2*z^2-504*x^2*z
            ....:            +x*y^3*z+6*x*y^3+154*x*y^2*z^2+6*x*y^2*z
            ....:            +36*x*y^2+924*x*y*z^2-11*y^4*z-66*y^3*z);
            sage: gI = R * (I.groebner_basis())
            sage: len(gI.gens())
            3
        """
        from sage.libs.singular.function import singular_function
        from sage.libs.singular.function_factory import ff
        from sage.libs.singular.option import opt
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import (
            slimgb_libsingular, std_libsingular)
        groebner = ff.groebner

        if get_verbose() >= 2:
            opt['prot'] = True
        for name, value in kwds.items():
            if value is not None:
                opt[name] = value

        if algorithm == "std":
            S = std_libsingular(self)
        elif algorithm == "slimgb":
            S = slimgb_libsingular(self)
        elif algorithm == "groebner":
            S = groebner(self)
        else:
            try:
                fnc = singular_function(algorithm)
                S = fnc(self)
            except NameError:
                raise NameError("Algorithm '%s' unknown" % algorithm)
        return S

    @libsingular_gb_standard_options
    def _groebner_cover(self):
        r"""
        Compute the Gröbner cover of the ideal.

        The Gröbner cover is a partition of the space of parameters,
        such that the Gröbner basis is constant for each of the parts.

        OUTPUT:

        A list of parts. Each element of this list contains:

        - The leading monomials of the Gröbner basis in the part
        - The Gröbner basis in the part
        - A list of components of the part. Each component is
          two lists of equations. The parameters in the
          part are those that satisfy the first list of
          equations, but do not satisfy the second one.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function, lib
            sage: F = PolynomialRing(QQ,'a').fraction_field()
            sage: F.inject_variables()
            Defining a
            sage: R.<x,y,z> = F[]
            sage: I = R.ideal([-x+3*y+z-5,2*x+a*z+4,4*x-3*z-a-1])
            sage: I._groebner_cover()
            [[[z, y, x],
            [(2*a + 3)*z + (a + 9), (12*a + 18)*y + (-a^2 - 23*a - 36), (4*a + 6)*x + (-a^2 - a + 12)],
            [[[0], [[(2*a + 3)]]]]],
            [[1], [1], [[[(2*a + 3)], [[1]]]]]]
        """
        from sage.rings.fraction_field import FractionField_generic
        from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        F = self.base_ring()
        if (not isinstance(F, FractionField_generic) or
            not isinstance(F.ring(), (MPolynomialRing_base, PolynomialRing_general))):
            raise TypeError("the base ring must be a field with parameters")
        from sage.arith.functions import lcm
        from sage.libs.singular.function import lib, singular_function
        lib("grobcov.lib")
        grobcov = singular_function("grobcov")
        polynomials = [f * lcm([c.denominator() for c in f.coefficients()])
                       for f in self.gens()]
        return grobcov(self.ring().ideal(polynomials))


class MPolynomialIdeal_singular_repr(
        MPolynomialIdeal_singular_base_repr):
    """
    An ideal in a multivariate polynomial ring, which has an
    underlying Singular ring associated to it.
    """
    def _singular_(self, singular=None):
        """
        Return Singular ideal corresponding to this ideal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal([x^3 + y, y])
            sage: S = I._singular_()
            sage: S
            x^3+y,
            y
        """
        if singular is None:
            singular = singular_default
        try:
            self.ring()._singular_(singular).set_ring()
            I = self.__singular
            if I.parent() is not singular:
                raise ValueError
            I._check_valid()
            return I
        except (AttributeError, ValueError):
            self.ring()._singular_(singular).set_ring()
            try:
                # this may fail for quotient ring elements, but is
                # faster
                gens = [str(x) for x in self.gens()]
                if len(gens) == 0:
                    gens = ['0']
                self.__singular = singular.ideal(gens)
            except TypeError:
                gens = [str(x.lift()) for x in self.gens()]
                if len(gens) == 0:
                    gens = ['0']
                self.__singular = singular.ideal(gens)
        return self.__singular

    @cached_method
    def _groebner_strategy(self):
        """
        Return Singular's Groebner Strategy object for the Groebner
        basis of this ideal which implements some optimized functions.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: I = R.ideal([y^3 - x^2])
            sage: I._groebner_strategy()
            Groebner Strategy for ideal generated by 1 elements over
            Multivariate Polynomial Ring in x, y over Rational Field

        .. NOTE::

            This function is mainly used internally.
        """
        from sage.libs.singular.groebner_strategy import GroebnerStrategy

        return GroebnerStrategy(MPolynomialIdeal(self.ring(), self.groebner_basis()))

    def plot(self, singular=None):
        r"""
        If you somehow manage to install surf, perhaps you can use
        this function to implicitly plot the real zero locus of this
        ideal (if principal).

        INPUT:

        - ``self`` -- must be a principal ideal in 2 or 3 vars over `\QQ`

        EXAMPLES:

        Implicit plotting in 2-d::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: I = R.ideal([y^3 - x^2])
            sage: I.plot()        # cusp
            Graphics object consisting of 1 graphics primitive
            sage: I = R.ideal([y^2 - x^2 - 1])
            sage: I.plot()        # hyperbola
            Graphics object consisting of 1 graphics primitive
            sage: I = R.ideal([y^2 + x^2*(1/4) - 1])
            sage: I.plot()        # ellipse
            Graphics object consisting of 1 graphics primitive
            sage: I = R.ideal([y^2-(x^2-1)*(x-2)])
            sage: I.plot()        # elliptic curve
            Graphics object consisting of 1 graphics primitive

        Implicit plotting in 3-d::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = R.ideal([y^2 + x^2*(1/4) - z])
            sage: I.plot()          # a cone; optional - surf
            sage: I = R.ideal([y^2 + z^2*(1/4) - x])
            sage: I.plot()          # same code, from a different angle; optional - surf
            sage: I = R.ideal([x^2*y^2+x^2*z^2+y^2*z^2-16*x*y*z])
            sage: I.plot()          # Steiner surface; optional - surf

        AUTHORS:

        - David Joyner (2006-02-12)
        """
        if self.ring().characteristic() != 0:
            raise TypeError("base ring must have characteristic 0")
        if not self.is_principal():
            raise TypeError("self must be principal")
        if singular is None:
            singular = singular_default
        singular.lib('surf')
        I = singular(self)
        I.plot()

    @require_field
    @cached_method
    @libsingular_gb_standard_options
    def complete_primary_decomposition(self, algorithm='sy'):
        r"""
        Return a list of primary ideals such that their intersection
        is ``self``, together with the associated prime ideals.

        An ideal `Q` is called primary if it is a proper ideal of the
        ring `R`, and if whenever `ab \in Q` and `a \not\in Q`, then
        `b^n \in Q` for some `n \in \ZZ`.

        If `Q` is a primary ideal of the ring `R`, then the radical
        ideal `P` of `Q` (i.e., the ideal consisting of all `a \in R`
        with a^n \in Q` for some `n \in \ZZ`), is called the
        associated prime of `Q`.

        If `I` is a proper ideal of a Noetherian ring `R`, then there
        exists a finite collection of primary ideals `Q_i` such that
        the following hold:

        - the intersection of the `Q_i` is `I`;

        - none of the `Q_i` contains the intersection of the others;

        - the associated prime ideals `P_i` of the `Q_i` are pairwise
          distinct.

        INPUT:

        - ``algorithm`` -- string; one of

          - ``'sy'`` -- (default) use the Shimoyama-Yokoyama
            algorithm

          - ``'gtz'`` -- use the Gianni-Trager-Zacharias algorithm

        OUTPUT:

        - a list of pairs `(Q_i, P_i)`, where the `Q_i` form a primary
          decomposition of ``self`` and `P_i` is the associated prime
          of `Q_i`.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2) * R
            sage: pd = I.complete_primary_decomposition(); sorted(pd, key=str)
            [(Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field),
             (Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field)]

            sage: pdc = I.primary_decomposition_complete(algorithm = 'gtz'); sorted(pdc, key=str)
            [(Ideal (z^2 + 1, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^2 + 1, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field),
             (Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field)]

            sage: from functools import reduce
            sage: reduce(lambda Qi,Qj: Qi.intersection(Qj), [Qi for (Qi,radQi) in pd]) == I
            True

            sage: [Qi.radical() == radQi for (Qi,radQi) in pd]
            [True, True]

            sage: P.<x,y,z> = PolynomialRing(ZZ)
            sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
            sage: I.complete_primary_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring must be a field for function 'complete_primary_decomposition'.

        ALGORITHM:

        Uses Singular.

        .. NOTE::

            See [BW1993]_ for an introduction to primary decomposition.

        TESTS:

        Check that :issue:`15745` is fixed::

            sage: R.<x,y>= QQ[]
            sage: I = Ideal(R(1))
            sage: I.complete_primary_decomposition()
            []
            sage: I.is_prime()
            False
        """

        # Avoid a bug in Singular (see #15745).
        if self.is_one():
            return []

        from sage.libs.singular.function_factory import ff

        if algorithm == 'sy':
            primdecSY = ff.primdec__lib.primdecSY
            P = primdecSY(self)
        elif algorithm == 'gtz':
            primdecGTZ = ff.primdec__lib.primdecGTZ
            P = primdecGTZ(self)

        R = self.ring()
        V = [(R.ideal(X[0]), R.ideal(X[1])) for X in P]
        return Sequence(V)

    # Seems useful for Tab-Completion
    primary_decomposition_complete = complete_primary_decomposition

    @require_field
    def primary_decomposition(self, algorithm='sy'):
        r"""
        Return a list of primary ideals such that their intersection
        is ``self``.

        An ideal `Q` is called primary if it is a proper ideal of the
        ring `R`, and if whenever `ab \in Q` and `a \not\in Q`, then
        `b^n \in Q` for some `n \in \ZZ`.

        If `Q` is a primary ideal of the ring `R`, then the radical
        ideal `P` of `Q` (i.e., the ideal consisting of all `a \in R`
        with `a^n \in Q` for some `n \in \ZZ`), is called the
        associated prime of `Q`.

        If `I` is a proper ideal of a Noetherian ring `R`, then there
        exists a finite collection of primary ideals `Q_i` such that
        the following hold:

        - the intersection of the `Q_i` is `I`;

        - none of the `Q_i` contains the intersection of the others;

        - the associated prime ideals of the `Q_i` are pairwise
          distinct.

        INPUT:

        - ``algorithm`` -- string; one of

          - ``'sy'`` -- (default) use the Shimoyama-Yokoyama
            algorithm

          - ``'gtz'`` -- use the Gianni-Trager-Zacharias algorithm

        OUTPUT:

        - a list of primary ideals `Q_i` forming a primary
          decomposition of ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2) * R
            sage: pd = I.primary_decomposition(); sorted(pd, key=str)
            [Ideal (z^2 + 1, y + 1)
              of Multivariate Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^6 + 4*z^3 + 4, y - z^2)
              of Multivariate Polynomial Ring in x, y, z over Rational Field]

        ::

            sage: from functools import reduce
            sage: reduce(lambda Qi, Qj: Qi.intersection(Qj), pd) == I
            True

        ALGORITHM:

        Uses Singular.

        REFERENCES:

        - Thomas Becker and Volker Weispfenning. Groebner Bases - A
          Computational Approach To Commutative Algebra. Springer, New
          York 1993.
        """
        return [I for I, _ in self.complete_primary_decomposition(algorithm)]

    @require_field
    def associated_primes(self, algorithm='sy'):
        r"""
        Return a list of the associated primes of primary ideals of
        which the intersection is `I` = ``self``.

        An ideal `Q` is called primary if it is a proper ideal of
        the ring `R` and if whenever `ab \in Q` and
        `a \not\in Q` then `b^n \in Q` for some
        `n \in \ZZ`.

        If `Q` is a primary ideal of the ring `R`, then the
        radical ideal `P` of `Q`, i.e.
        `P = \{a \in R, a^n \in Q\}` for some
        `n \in \ZZ`, is called the
        *associated prime* of `Q`.

        If `I` is a proper ideal of the ring `R` then there
        exists a decomposition in primary ideals `Q_i` such that

        - their intersection is `I`

        - none of the `Q_i` contains the intersection of the
          rest, and

        - the associated prime ideals of `Q_i` are pairwise
          different.


        This method returns the associated primes of the `Q_i`.

        INPUT:

        - ``algorithm`` -- string; one of

          - ``'sy'`` -- (default) use the Shimoyama-Yokoyama algorithm

          - ``'gtz'`` -- use the Gianni-Trager-Zacharias algorithm

        OUTPUT: list of associated primes

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2) * R
            sage: pd = I.associated_primes(); sorted(pd, key=str)
            [Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field]

        ALGORITHM:

        Uses Singular.

        REFERENCES:

        - Thomas Becker and Volker Weispfenning. Groebner Bases - A
          Computational Approach To Commutative Algebra. Springer, New
          York 1993.
        """
        return [P for _,P in self.complete_primary_decomposition(algorithm)]

    def is_prime(self, **kwds):
        r"""
        Return ``True`` if this ideal is prime.

        INPUT:

        - keyword arguments are passed on to
          :meth:`complete_primary_decomposition`; in this way you can
          specify the algorithm to use.

        EXAMPLES::

            sage: R.<x, y> = PolynomialRing(QQ, 2)
            sage: I = (x^2 - y^2 - 1) * R
            sage: I.is_prime()
            True
            sage: (I^2).is_prime()
            False

            sage: J = (x^2 - y^2) * R
            sage: J.is_prime()
            False
            sage: (J^3).is_prime()
            False

            sage: (I * J).is_prime()
            False

        The following is :issue:`5982`.  Note that the quotient ring
        is not recognized as being a field at this time, so the
        fraction field is not the quotient ring itself::

            sage: Q = R.quotient(I); Q
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field
             by the ideal (x^2 - y^2 - 1)
            sage: Q.fraction_field()
            Fraction Field of
             Quotient of Multivariate Polynomial Ring in x, y over Rational Field
              by the ideal (x^2 - y^2 - 1)
        """
        if not self.ring().base_ring().is_field():
            raise NotImplementedError
        CPD = self.complete_primary_decomposition(**kwds)
        if len(CPD) != 1:
            return False
        _, P = CPD[0]
        return self == P

    @require_field
    @handle_AA_and_QQbar
    @singular_gb_standard_options
    @libsingular_gb_standard_options
    def triangular_decomposition(self, algorithm=None, singular=None):
        """
        Decompose zero-dimensional ideal ``self`` into triangular
        sets.

        This requires that the given basis is reduced w.r.t. to the
        lexicographical monomial ordering. If the basis of ``self`` does
        not have this property, the required Groebner basis is
        computed implicitly.

        INPUT:

        - ``algorithm`` -- string or ``None`` (default: ``None``)

        ALGORITHMS:

        - ``'singular:triangL'`` -- decomposition of ``self`` into triangular
          systems (Lazard)

        - ``'singular:triangLfak'`` -- decomposition of ``self`` into triangular systems
          plus factorization

        - ``'singular:triangM'`` -- decomposition of ``self`` into
          triangular systems (Moeller)

        OUTPUT: list `T` of lists `t` such that the variety of
        ``self`` is the union of the varieties of `t` in `L` and each
        `t` is in triangular form

        EXAMPLES::

            sage: P.<e,d,c,b,a> = PolynomialRing(QQ, 5, order='lex')
            sage: I = sage.rings.ideal.Cyclic(P)
            sage: GB = Ideal(I.groebner_basis('libsingular:stdfglm'))
            sage: GB.triangular_decomposition('singular:triangLfak')
            [Ideal (a - 1, b - 1, c - 1, d^2 + 3*d + 1, e + d + 3) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a - 1, b - 1, c^2 + 3*c + 1, d + c + 3, e - 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a - 1, b^2 + 3*b + 1, c + b + 3, d - 1, e - 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a - 1, b^4 + b^3 + b^2 + b + 1, -c + b^2, -d + b^3,
                    e + b^3 + b^2 + b + 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^2 + 3*a + 1, b - 1, c - 1, d - 1, e + a + 3) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^2 + 3*a + 1, b + a + 3, c - 1, d - 1, e - 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 - 4*a^3 + 6*a^2 + a + 1,
                    -11*b^2 + 6*b*a^3 - 26*b*a^2 + 41*b*a - 4*b - 8*a^3 + 31*a^2 - 40*a - 24,
                    11*c + 3*a^3 - 13*a^2 + 26*a - 2, 11*d + 3*a^3 - 13*a^2 + 26*a - 2,
                    -11*e - 11*b + 6*a^3 - 26*a^2 + 41*a - 4) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 + a^3 + a^2 + a + 1,
                    b - 1, c + a^3 + a^2 + a + 1, -d + a^3, -e + a^2) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 + a^3 + a^2 + a + 1,
                    b - a, c - a, d^2 + 3*d*a + a^2, e + d + 3*a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 + a^3 + a^2 + a + 1,
                    b - a, c^2 + 3*c*a + a^2, d + c + 3*a, e - a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 + a^3 + a^2 + a + 1,
                    b^2 + 3*b*a + a^2, c + b + 3*a, d - a, e - a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 + a^3 + a^2 + a + 1,
                    b^3 + b^2*a + b^2 + b*a^2 + b*a + b + a^3 + a^2 + a + 1,
                    c + b^2*a^3 + b^2*a^2 + b^2*a + b^2,
                    -d + b^2*a^2 + b^2*a + b^2 + b*a^2 + b*a + a^2,
                    -e + b^2*a^3 - b*a^2 - b*a - b - a^2 - a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
             Ideal (a^4 + a^3 + 6*a^2 - 4*a + 1,
                    -11*b^2 + 6*b*a^3 + 10*b*a^2 + 39*b*a + 2*b + 16*a^3 + 23*a^2 + 104*a - 24,
                    11*c + 3*a^3 + 5*a^2 + 25*a + 1, 11*d + 3*a^3 + 5*a^2 + 25*a + 1,
                    -11*e - 11*b + 6*a^3 + 10*a^2 + 39*a + 2) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field]

            sage: R.<x1,x2> = PolynomialRing(QQ, 2, order='lex')
            sage: f1 = 1/2*((x1^2 + 2*x1 - 4)*x2^2 + 2*(x1^2 + x1)*x2 + x1^2)
            sage: f2 = 1/2*((x1^2 + 2*x1 + 1)*x2^2 + 2*(x1^2 + x1)*x2 - 4*x1^2)
            sage: I = Ideal(f1,f2)
            sage: I.triangular_decomposition()
            [Ideal (x2, x1^2) of Multivariate Polynomial Ring in x1, x2 over Rational Field,
             Ideal (x2, x1^2) of Multivariate Polynomial Ring in x1, x2 over Rational Field,
             Ideal (x2, x1^2) of Multivariate Polynomial Ring in x1, x2 over Rational Field,
             Ideal (x2^4 + 4*x2^3 - 6*x2^2 - 20*x2 + 5, 8*x1 - x2^3 + x2^2 + 13*x2 - 5)
              of Multivariate Polynomial Ring in x1, x2 over Rational Field]

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: J = Ideal(x^2+y^2-2, y^2-1)
            sage: J.triangular_decomposition()
            [Ideal (y^2 - 1, x^2 - 1) of Multivariate Polynomial Ring in x, y over Rational Field]

        Check that this method works over QQbar (:issue:`25351`)::

            sage: R.<x,y> = QQbar[]                                                     # needs sage.rings.number_field
            sage: J = Ideal(x^2 + y^2 - 2, y^2 - 1)                                     # needs sage.rings.number_field
            sage: J.triangular_decomposition()                                          # needs sage.rings.number_field
            [Ideal (y^2 - 1, x^2 - 1) of Multivariate Polynomial Ring in x, y over Algebraic Field]
        """
        P = self.ring()

        is_groebner = self.basis_is_groebner()

        # make sure to work w.r.t. 'lex'
        if P.term_order() != 'lex':
            Q = P.change_ring(order='lex')
        else:
            Q = P

        # the Singular routines are quite picky about their input.
        if is_groebner:
            if Q == P:
                I = MPolynomialIdeal(P, self.interreduced_basis()[::-1])
            else:
                I = self
                I = MPolynomialIdeal(P, I.transformed_basis('fglm')[::-1])  # -> 'lex'
                I = I.change_ring(Q)  # transform to 'lex' GB
        else:
            if Q == P:
                I = MPolynomialIdeal(P, self.groebner_basis()[::-1])
            else:
                I = self.change_ring(Q)  # transform to 'lex' GB
                I = MPolynomialIdeal(Q, I.groebner_basis()[::-1])

        if I.dimension() != 0:
            raise TypeError("dimension must be zero")

        from sage.libs.singular.function import lib as singular_lib
        from sage.libs.singular.function import singular_function

        singular_lib('triang.lib')

        if algorithm is None:
            algorithm = "singular:triangL"

        if algorithm in ("singular:triangL","singular:triangLfak","singular:triangM"):
            f = singular_function(algorithm[9:])
            Tbar = f(I, attributes={I:{'isSB':1}})
        else:
            raise TypeError("algorithm '%s' unknown" % algorithm)

        T = Sequence([ MPolynomialIdeal(Q,t) for t in Tbar])
        return sorted(T, key=lambda x: x.gens())

    @require_field
    @handle_AA_and_QQbar
    def dimension(self, singular=None):
        """
        The dimension of the ring modulo this ideal.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(32003), order='degrevlex')              # needs sage.rings.finite_rings
            sage: I = ideal(x^2 - y, x^3)                                               # needs sage.rings.finite_rings
            sage: I.dimension()                                                         # needs sage.rings.finite_rings
            1

        If the ideal is the total ring, the dimension is `-1` by convention.

        ALGORITHM:

        For principal ideals, Theorem 3.5.1 of [Ger2008]_ is used. Otherwise
        Singular is used, unless the characteristic is too large. This requires
        computation of a Groebner basis, which can be very expensive.

        For polynomials over a finite field of order too large for Singular,
        this falls back on a toy implementation of Buchberger to compute the
        Groebner basis, then uses the algorithm described in Chapter 9, Section
        1 of Cox, Little, and O'Shea's "Ideals, Varieties, and Algorithms".

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = PolynomialRing(GF(2147483659^2), order='lex')
            sage: I = R.ideal([x*y, x*y + 1])
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            -1
            sage: I=ideal([x*(x*y+1), y*(x*y+1)])
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            1
            sage: I = R.ideal([x^3*y, x*y^2])
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            1
            sage: R.<x,y> = PolynomialRing(GF(2147483659^2), order='lex')
            sage: I = R.ideal(0)
            sage: I.dimension()
            2

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: P.<x,y,z> = QQbar[]                                                   # needs sage.rings.number_field
            sage: I = ideal(x^2-y, x^3-QQbar(-1))                                       # needs sage.rings.number_field
            sage: I.dimension()                                                         # needs sage.rings.number_field
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            if not self.base_ring().is_field():
                raise NotImplementedError("implemented only over fields")
            if self.ngens() == 1:
                g = self.gen(0)
                if g.is_unit():
                    self.__dimension = Integer(-1)
                elif g:
                    self.__dimension = Integer(self.ring().ngens() - 1)
                else:
                    self.__dimension = Integer(self.ring().ngens())
                return self.__dimension
            try:
                from sage.libs.singular.function_factory import ff
                dim = ff.dim
                v = MPolynomialIdeal(self.ring(),self.groebner_basis())
                self.__dimension = Integer(dim(v, attributes={v:{'isSB':1}}))
            except TypeError:
                try:
                    v = self._groebner_basis_singular_raw()
                    self.__dimension = Integer(v.dim())
                except TypeError:
                    verbose("Warning: falling back to very slow toy implementation.", level=0)
                    if not self.ring().term_order().is_global():
                        raise TypeError("local/unknown ordering is not supported by the toy implementation")
                    # See Chapter 9, Section 1 of Cox, Little, O'Shea's
                    # "Ideals, Varieties, and Algorithms"
                    from sage.sets.set import Set
                    gb = toy_buchberger.buchberger_improved(self)
                    if self.ring().one() in gb:
                        self.__dimension = Integer(-1)
                        return self.__dimension
                    ring_vars = self.ring().gens()
                    n = len(ring_vars)
                    lms = [each.lm() for each in gb]
                    # compute M_j, denoted by var_lms
                    var_lms = [Set([]) for _ in lms]
                    for j in range(len(ring_vars)):
                        for i in range(len(lms)):
                            if lms[i].degree(ring_vars[j]) > 0:
                                var_lms[i] += Set([j+1])
                    # compute intersections of M_j and J
                    # we assume that the iterator starts with the empty set,
                    # then iterates through all subsets of order 1,
                    # then through all subsets of order 2, etc...
                    # the way Sage currently operates
                    all_J = Set([each + 1 for each in range(n)]).subsets()
                    min_dimension = -1
                    all_J = iter(all_J)
                    while min_dimension == -1:
                        try:
                            J = next(all_J)
                        except StopIteration:
                            min_dimension = n
                            break
                        J_intersects_all = True
                        i = 0
                        while J_intersects_all and i < len(var_lms):
                            J_intersects_all = J.intersection(var_lms[i]) != Set([])
                            i += 1
                        if J_intersects_all:
                            min_dimension = len(J)
                    self.__dimension = Integer(n - min_dimension)
            return self.__dimension

    @require_field
    @handle_AA_and_QQbar
    def vector_space_dimension(self):
        """
        Return the vector space dimension of the ring modulo this ideal. If
        the ideal is not zero-dimensional, a :exc:`TypeError` is raised.

        ALGORITHM:

        Uses Singular.

        EXAMPLES::

            sage: R.<u,v> = PolynomialRing(QQ)
            sage: g = u^4 + v^4 + u^3 + v^3
            sage: I = ideal(g) + ideal(g.gradient())
            sage: I.dimension()
            0
            sage: I.vector_space_dimension()
            4

        When the ideal is not zero-dimensional, we return infinity::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: I = R.ideal(x)
            sage: I.dimension()
            1
            sage: I.vector_space_dimension()
            +Infinity

        Due to integer overflow, the result is correct only modulo ``2^32``, see :issue:`8586`::

            sage: P.<x,y,z> = PolynomialRing(GF(32003), 3)                              # needs sage.rings.finite_rings
            sage: sage.rings.ideal.FieldIdeal(P).vector_space_dimension()       # known bug, needs sage.rings.finite_rings
            32777216864027

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: P.<x,y,z> = QQbar[]
            sage: I = ideal(x^2-y,x^3-QQbar(-1),z-y)
            sage: I.dimension()
            0
            sage: I.vector_space_dimension()
            3
        """
        R = self.ring()
        gb = R.ideal(self.groebner_basis())

        from sage.libs.singular.function_factory import ff
        vdim = ff.vdim
        vd = Integer(vdim(gb, attributes={gb:{'isSB':1}}))

        if vd == -1:
            from sage.rings.infinity import Infinity
            return Infinity
        else:
            return vd

    @require_field
    def _groebner_basis_ginv(self, algorithm='TQ', criteria='CritPartially', division_interface='Janet'):
        r"""
        Compute a Groebner basis using GINV.

        INPUT:

        - ``algorithm`` -- ``'TQ'``, ``'TQBlockHigh'``, ``'TQBlockLow'``, or
          ``'TQDegree'``
        - ``criteria`` -- ``'Without'`` (without any criteria),
          ``'C1'``, ``'CritPartially'`` (partial involutive criteria),
          ``'C1C2C3'``, ``'C1C2C3C4'`` (full involutive criteria)

        - ``division_interface`` -- either ``'Janet'`` or ``'JanetLike'``

        EXAMPLES:

        Currently, only `\GF{p}` and `\QQ` are supported as base fields::

            sage: # optional - ginv
            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I.groebner_basis(algorithm='ginv')
            [z^3 - 79/210*z^2 + 1/30*y + 1/70*z, y^2 - 3/5*z^2 - 1/5*y + 1/5*z, y*z + 6/5*z^2 - 1/10*y - 2/5*z, x + 2*y + 2*z - 1]
            sage: P.<x,y,z> = PolynomialRing(GF(127), order='degrevlex')
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I.groebner_basis(algorithm='ginv')
            ...
            [z^3 + 22*z^2 - 55*y + 49*z, y^2 - 26*z^2 - 51*y + 51*z, y*z + 52*z^2 + 38*y + 25*z, x + 2*y + 2*z - 1]

        .. NOTE::

            Criterion C1 is Buchberger's co-prime criterion. Criteria
            C2, C3 and C4 in the aggregate are equivalent to the second
            (chain) Buchberger's criterion. Supported term orderings are
            'lex' and 'degrevlex', supported base rings are `\GF{p}` with
            `p < 2^16` and `\QQ`.
        """
        P = self.ring()
        T = P.term_order()
        K = P.base_ring()

        import ginv

        st = ginv.SystemType("Polynomial")

        term_order_map = {'degrevlex':"DegRevLex",'lex':"Lex"}
        try:
            im = ginv.MonomInterface(term_order_map[T.name()], st, list(P.variable_names()))
        except KeyError:
            raise NotImplementedError("Term order '%s' not supported by Sage's GINV interface or GINV" % T.term_order())

        from sage.rings.rational_field import QQ
        if K is QQ:
            ic = ginv.CoeffInterface("GmpQ", st)
        elif K.order() <= 2**16 and K.order().is_prime():
            ic = ginv.CoeffInterface("ModularShort", st, modularShort=K.order())
        else:
            raise NotImplementedError("GINV interface for base ring '%s' is not implemented." % K)

        ip = ginv.PolyInterface("PolyList", st, im, ic)
        iw = ginv.WrapInterface(criteria, ip)
        iD = ginv.DivisionInterface(division_interface, iw)

        system = [ginv.Poly(ip, str(f)) for f in self.gens()]
        G = ginv.basisBuild(algorithm, iD, system)
        G = Sequence([P(str(f)) for f in G.iterGB()])
        return G

    @singular_gb_standard_options
    def _groebner_basis_singular(self, algorithm='groebner', *args, **kwds):
        """
        Return the reduced Groebner basis of this ideal. If the
        Groebner basis for this ideal has been calculated before, the
        cached Groebner basis is returned regardless of the requested
        algorithm.

        INPUT:

        - ``algorithm`` -- see below for available algorithms

        ALGORITHMS:

        'groebner'
            use Singular's groebner heuristic to choose an algorithm (default)

        'std'
            Buchberger's algorithm

        'stdhilb'
            computes the standard basis of the homogeneous ideal in the
            base ring, via a Hilbert driven standard basis computation.

        'stdfglm'
            computes the standard basis of the ideal in the base ring via fglm
            (from the degrevlex ordering to the ordering of the base ring).

        'slimgb'
            the *SlimGB* algorithm

        EXAMPLES:

        We compute a Groebner basis of 'cyclic 4' relative to
        lexicographic ordering. ::

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4); I
            Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d
            + a*c*d + b*c*d, a*b*c*d - 1) of Multivariate Polynomial
            Ring in a, b, c, d over Rational Field

        ::

            sage: I._groebner_basis_singular()
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d,
             b*d^4 - b + d^5 - d, b*c - b*d + c^2*d^4 + c*d - 2*d^2,
             b^2 + 2*b*d + d^2, a + b + c + d]

        ALGORITHM:

        Uses Singular.

        .. NOTE::

            This method is called by the :meth:`.groebner_basis` method
            and the user usually doesn't need to bother with this one.
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence

        R = self.ring()
        S = self._groebner_basis_singular_raw(algorithm=algorithm, *args, **kwds)
        S = PolynomialSequence([R(S[i+1]) for i in range(len(S))], R, immutable=True)
        return S

    @cached_method
    def _groebner_basis_singular_raw(self, algorithm='groebner', singular=None, *args, **kwds):
        r"""
        Return a Groebner basis in Singular format.

        INPUT:

        - ``algorithm`` -- Singular function to call (default: ``groebner``)

        - ``singular`` -- Singular instance to use (default: ``singular_default``)

        - ``args`` -- ignored

        - ``kwds`` -- Singular options

        EXAMPLES::

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: I._groebner_basis_singular()  # indirect doctest
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d,
             b*d^4 - b + d^5 - d, b*c - b*d + c^2*d^4 + c*d - 2*d^2,
             b^2 + 2*b*d + d^2, a + b + c + d]
        """
        if singular is None:
            singular = singular_default
        #try:
        #    return self.__gb_singular
        #except AttributeError:
        #    pass
        # singular options are preserved by @singular_gb_standard_options,
        # so we don't need to do that here too
        from sage.libs.singular.option import _options_py_to_singular
        S = self._singular_()   # for degBound, we need to ensure
                                # that a ring is defined

        if get_verbose() >= 2:
            kwds['prot'] = True

        for o, v in kwds.items():
            o = _options_py_to_singular.get(o,o)
            if v:
                if o in ['degBound','multBound']:
                    singular.eval(o+'=%d' % v)
                else:
                    singular.option(o)
            else:
                if o in ['degBound','multBound']:
                    singular.eval(o+'=0')
                else:
                    singular.option("no"+o)

        obj = self._singular_()

        prot = kwds.get('prot',False)

        if prot == "sage":
            if algorithm == 'slimgb':
                warn("'slimgb' does not print sufficient information for prot='sage' to work reliably, the highest degree reached might be too low.")
            from sage.interfaces.singular import SingularGBLogPrettyPrinter
            log_parser = SingularGBLogPrettyPrinter(verbosity=get_verbose()+1)
        else:
            log_parser = None

        ctx = StdOutContext(singular, silent=not prot, stdout=log_parser)

        with ctx:
            if algorithm == "groebner":
                S = obj.groebner()
            elif algorithm == "std":
                S = obj.std()
            elif algorithm == "slimgb":
                S = obj.slimgb()
            elif algorithm == "stdhilb":
                S = obj.stdhilb()
            elif algorithm == "stdfglm":
                S = obj.stdfglm()
            else:
                raise TypeError("algorithm '%s' unknown" % algorithm)
        self.__gb_singular = S
        if prot == "sage":
            print("")
            print("Highest degree reached during computation: %2d." % log_parser.max_deg)
        return S

    @require_field
    @cached_method
    @handle_AA_and_QQbar
    def genus(self):
        r"""
        Return the geometric genus of the projective curve defined by this
        ideal.

        OUTPUT:

        If the ideal is homogeneous and defines a curve in a projective space,
        then the genus of the curve is returned. If the ideal is not
        homogeneous and defines a curve in an affine space, the genus of the
        projective closure of the curve is returned.

        EXAMPLES:

        Consider the hyperelliptic curve `y^2 = 4x^5 - 30x^3 + 45x - 22` over
        `\QQ`, it has genus 2::

            sage: P.<x> = QQ[]
            sage: f = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C = HyperellipticCurve(f); C
            Hyperelliptic Curve over Rational Field defined by y^2 = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C.genus()
            2

        ::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: f = y^2 - 4*x^5 - 30*x^3 + 45*x - 22
            sage: I = Ideal([f])
            sage: I.genus()
            2

        Geometric genus is only defined for geometrically irreducible curves.
        You may get a nonsensical answer if the condition is not met.  A curve
        reducible over a quadratic extension of `\QQ`::

            sage: R.<x, y, z> = QQ[]
            sage: C = Curve(x^2 - 2*y^2)
            sage: C.genus()
            -1

        TESTS:

        An ideal that does not define a curve but we get a result! ::

            sage: R.<x, y, z> = QQ[]
            sage: Ideal(x^4 + y^2*x + x).genus()
            0

        An ideal that defines a geometrically reducible affine curve::

            sage: T.<t1,t2,u1,u2> = QQ[]
            sage: TJ = Ideal([t1^2 + u1^2 - 1,t2^2 + u2^2 - 1, (t1-t2)^2 + (u1-u2)^2 -1])
            sage: TJ.genus()
            -1

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = QQbar[]
            sage: I = ideal(y^3*z + x^3*y + x*z^3)
            sage: I.genus()
            3
        """
        from sage.libs.singular.function_factory import ff
        genus = ff.normal__lib.genus
        return Integer(genus(self))

    @handle_AA_and_QQbar
    @libsingular_gb_standard_options
    def intersection(self, *others):
        """
        Return the intersection of the arguments with this ideal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2, order='lex')
            sage: I = x*R
            sage: J = y*R
            sage: I.intersection(J)
            Ideal (x*y) of Multivariate Polynomial Ring in x, y over Rational Field

        The following simple example illustrates that the product need
        not equal the intersection. ::

            sage: I = (x^2, y) * R
            sage: J = (y^2, x) * R
            sage: K = I.intersection(J); K
            Ideal (y^2, x*y, x^2) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: IJ = I*J; IJ
            Ideal (x^2*y^2, x^3, y^3, x*y)
             of Multivariate Polynomial Ring in x, y over Rational Field
            sage: IJ == K
            False

        Intersection of several ideals::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: I1 = x * R
            sage: I2 = y * R
            sage: I3 = (x, y) * R
            sage: I4 = (x^2 + x*y*z, y^2 - z^3*y, z^3 + y^5*x*z) * R
            sage: I1.intersection(I2, I3, I4).groebner_basis()
            [x^2*y + x*y*z^4, x*y^2 - x*y*z^3, x*y*z^20 - x*y*z^3]

        The ideals must share the same ring::

            sage: R2.<x,y> = PolynomialRing(QQ, 2, order='lex')
            sage: R3.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: I2 = x*R2
            sage: I3 = x*R3
            sage: I2.intersection(I3)
            Traceback (most recent call last):
            ...
            TypeError: Intersection is only available for ideals of the same ring.

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: R.<x,y> = QQbar[]
            sage: I = x*R
            sage: J = y*R
            sage: I.intersection(J)
            Ideal (x*y) of Multivariate Polynomial Ring in x, y over Algebraic Field
        """
        R = self.ring()

        for other in others:
            if not isinstance(other, MPolynomialIdeal_singular_repr) or other.ring() != R:
                raise TypeError("Intersection is only available for ideals of the same ring.")

        from sage.libs.singular.function_factory import ff
        intersect = ff.intersect

        K = intersect(self, *others)
        return R.ideal(K)

    @require_field
    @libsingular_gb_standard_options
    def minimal_associated_primes(self):
        """
        OUTPUT: list of prime ideals

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, 'xyz')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2) * R
            sage: sorted(I.minimal_associated_primes(), key=str)
            [Ideal (z^2 + 1, -z^2 + y)
              of Multivariate Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^3 + 2, -z^2 + y)
              of Multivariate Polynomial Ring in x, y, z over Rational Field]

        ALGORITHM:

        Uses Singular.
        """
        from sage.libs.singular.function_factory import ff
        minAssGTZ = ff.primdec__lib.minAssGTZ

        M = minAssGTZ(self)
        R = self.ring()
        return [R.ideal(J) for J in M]

    @require_field
    @libsingular_gb_standard_options
    def radical(self):
        r"""
        The radical of this ideal.

        EXAMPLES:

        This is an obviously not radical ideal::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: I = (x^2, y^3, (x*z)^4 + y^3 + 10*x^2) * R
            sage: I.radical()
            Ideal (y, x) of Multivariate Polynomial Ring in x, y, z over Rational Field

        That the radical is correct is clear from the Groebner basis. ::

            sage: I.groebner_basis()
            [y^3, x^2]

        This is the example from the Singular manual::

            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2) * R
            sage: I.radical()
            Ideal (z^2 - y, y^2*z + y*z + 2*y + 2)
             of Multivariate Polynomial Ring in x, y, z over Rational Field

        .. NOTE::

            From the Singular manual: A combination of the algorithms
            of Krick/Logar and Kemper is used. Works also in positive
            characteristic (Kempers algorithm).

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y,z> = PolynomialRing(GF(37), 3)
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2) * R
            sage: I.radical()
            Ideal (z^2 - y, y^2*z + y*z + 2*y + 2)
             of Multivariate Polynomial Ring in x, y, z over Finite Field of size 37
        """
        from sage.libs.singular.function_factory import ff
        radical = ff.primdec__lib.radical
        r = radical(self)

        S = self.ring()

        #I = self._singular_()
        #I.parent().lib('primdec.lib')
        #r = I.radical()

        return S.ideal(r)

    @require_field
    @libsingular_gb_standard_options
    def integral_closure(self, p=0, r=True, singular=None):
        """
        Let `I` = ``self``.

        Return the integral closure of `I, ..., I^p`, where `sI` is
        an ideal in the polynomial ring `R=k[x(1),...x(n)]`. If `p` is
        not given, or `p=0`, compute the closure of all powers up to
        the maximum degree in `t` occurring in the closure of `R[It]`
        (so this is the last power whose closure is not just the
        sum/product of the smaller). If `r` is given and ``r is
        True``, ``I.integral_closure()`` starts with a check whether `I`
        is already a radical ideal.

        INPUT:

        - ``p`` -- powers of `I` (default: 0)

        - ``r`` -- check whether ``self`` is a radical ideal first (default: ``True``)

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = ideal([x^2, x*y^4, y^5])
            sage: I.integral_closure()
            [x^2, x*y^4, y^5, x*y^3]

        ALGORITHM:

        Uses libSINGULAR.
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence

        R = self.ring()
        from sage.libs.singular.function_factory import ff
        normalI = ff.reesclos__lib.normalI
        ret = PolynomialSequence(normalI(self, p, int(r))[0], R, immutable=True)
        return ret

    @require_field
    @handle_AA_and_QQbar
    def syzygy_module(self):
        r"""
        Compute the first syzygy (i.e., the module of relations of the
        given generators) of the ideal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 2*x^2 + y
            sage: g = y
            sage: h = 2*f + g
            sage: I = Ideal([f,g,h])
            sage: M = I.syzygy_module(); M
            [       -2        -1         1]
            [       -y 2*x^2 + y         0]
            sage: G = vector(I.gens())
            sage: M*G
            (0, 0)

        ALGORITHM:

        Uses Singular's syz command.

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: R.<x,y> = QQbar[]
            sage: f = 2*x^2 + y
            sage: g = y
            sage: h = 2*f + g
            sage: I = Ideal([f,g,h])
            sage: M = I.syzygy_module(); M
            [       -2        -1         1]
            [       -y 2*x^2 + y         0]
            sage: G = vector(I.gens())
            sage: M*G
            (0, 0)
        """
        from sage.libs.singular.function_factory import ff
        syz = ff.syz
        from sage.matrix.constructor import matrix

        #return self._singular_().syz().transpose().sage_matrix(self.ring())
        S = syz(self)
        return matrix(self.ring(), S)

    @require_field
    def free_resolution(self, *args, **kwds):
        r"""
        Return a free resolution of ``self``.

        For input options, see
        :class:`~sage.homology.free_resolution.FreeResolution`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 2*x^2 + y
            sage: g = y
            sage: h = 2*f + g
            sage: I = R.ideal([f,g,h])
            sage: res = I.free_resolution()
            sage: res
            S^1 <-- S^2 <-- S^1 <-- 0
            sage: ascii_art(res.chain_complex())
                                        [-x^2]
                        [  y x^2]       [   y]
             0 <-- C_0 <---------- C_1 <------- C_2 <-- 0

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = q.parent()[]
            sage: I = R.ideal([x^2+y^2+z^2, x*y*z^4])
            sage: I.free_resolution()
            Traceback (most recent call last):
            ...
            NotImplementedError: the ring must be a polynomial ring using Singular
        """
        from sage.rings.polynomial.multi_polynomial_libsingular import \
            MPolynomialRing_libsingular
        if isinstance(self.ring(), MPolynomialRing_libsingular):
            from sage.homology.free_resolution import \
                FiniteFreeResolution_singular
            return FiniteFreeResolution_singular(self, *args, **kwds)
        raise NotImplementedError("the ring must be a polynomial ring using Singular")

    @require_field
    def graded_free_resolution(self, *args, **kwds):
        r"""
        Return a graded free resolution of ``self``.

        For input options, see
        :class:`~sage.homology.graded_resolution.GradedFreeResolution`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 2*x^2 + y^2
            sage: g = y^2
            sage: h = 2*f + g
            sage: I = R.ideal([f,g,h])
            sage: I.graded_free_resolution(shifts=[1])
            S(-1) <-- S(-3)⊕S(-3) <-- S(-5) <-- 0

            sage: f = 2*x^2 + y
            sage: g = y
            sage: h = 2*f + g
            sage: I = R.ideal([f,g,h])
            sage: I.graded_free_resolution(degrees=[1,2])
            S(0) <-- S(-2)⊕S(-2) <-- S(-4) <-- 0

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = q.parent()[]
            sage: I = R.ideal([x^2+y^2+z^2, x*y*z^4])
            sage: I.graded_free_resolution()
            Traceback (most recent call last):
            ...
            NotImplementedError: the ring must be a polynomial ring using Singular
        """
        from sage.rings.polynomial.multi_polynomial_libsingular import \
            MPolynomialRing_libsingular
        if isinstance(self.ring(), MPolynomialRing_libsingular):
            from sage.homology.graded_resolution import \
                GradedFiniteFreeResolution_singular
            return GradedFiniteFreeResolution_singular(self, *args, **kwds)
        raise NotImplementedError("the ring must be a polynomial ring using Singular")

    @handle_AA_and_QQbar
    @singular_gb_standard_options
    @libsingular_gb_standard_options
    def interreduced_basis(self):
        r"""
        If this ideal is spanned by `(f_1, ..., f_n)`,
        return `(g_1, ..., g_s)` such that:

        - `(f_1,...,f_n) = (g_1,...,g_s)`

        - `LT(g_i) \neq LT(g_j)` for all `i \neq j`

        - `LT(g_i)` does not divide `m` for all monomials `m` of
          `\{g_1,...,g_{i-1},g_{i+1},...,g_s\}`

        - `LC(g_i) = 1` for all `i` if the coefficient ring is a field.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([z*x + y^3, z + y^3, z + x*y])
            sage: I.interreduced_basis()
            [y^3 + z, x*y + z, x*z - z]

        Note that tail reduction for local orderings is not well-defined::

            sage: R.<x,y,z> = PolynomialRing(QQ,order='negdegrevlex')
            sage: I = Ideal([z*x + y^3, z + y^3, z + x*y])
            sage: I.interreduced_basis()
            [z + x*y, x*y - y^3, x^2*y - y^3]

        A fixed error with nonstandard base fields::

            sage: R.<t> = QQ['t']
            sage: K.<x,y> = R.fraction_field()['x,y']
            sage: I = t*x * K
            sage: I.interreduced_basis()
            [x]

        The interreduced basis of 0 is 0::

            sage: P.<x,y,z> = GF(2)[]
            sage: Ideal(P(0)).interreduced_basis()
            [0]

        ALGORITHM:

        Uses Singular's ``interred`` command or
        :func:`sage.rings.polynomial.toy_buchberger.inter_reduction`
        if conversion to Singular fails.

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: R.<x,y,z> = QQbar[]                                                   # needs sage.rings.number_field
            sage: I = Ideal([z*x + y^3, z + y^3, z + x*y])                              # needs sage.rings.number_field
            sage: I.interreduced_basis()                                                # needs sage.rings.number_field
            [y^3 + z, x*y + z, x*z - z]
        """
        return self.basis.reduced()

    @cached_method
    @handle_AA_and_QQbar
    @singular_gb_standard_options
    def basis_is_groebner(self, singular=None):
        r"""
        Return ``True`` if the generators of this ideal
        (``self.gens()``) form a Groebner basis.

        Let `I` be the set of generators of this ideal. The check is
        performed by trying to lift `Syz(LM(I))` to `Syz(I)` as `I`
        forms a Groebner basis if and only if for every element `S` in
        `Syz(LM(I))`:

        .. MATH::

            S * G = \sum_{i=0}^{m} h_ig_i ---->_G 0.

        ALGORITHM:

        Uses Singular.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127), 10)
            sage: I = sage.rings.ideal.Cyclic(R, 4)
            sage: I.basis_is_groebner()
            False
            sage: I2 = Ideal(I.groebner_basis())
            sage: I2.basis_is_groebner()
            True

        A more complicated example::

            sage: R.<U6,U5,U4,U3,U2, u6,u5,u4,u3,u2, h> = PolynomialRing(GF(7583))      # needs sage.rings.finite_rings
            sage: l = [u6 + u5 + u4 + u3 + u2 - 3791*h,                                 # needs sage.rings.finite_rings
            ....:      U6 + U5 + U4 + U3 + U2 - 3791*h,
            ....:      U2*u2 - h^2, U3*u3 - h^2, U4*u4 - h^2,
            ....:      U5*u4 + U5*u3 + U4*u3 + U5*u2 + U4*u2 + U3*u2 - 3791*U5*h
            ....:        - 3791*U4*h - 3791*U3*h - 3791*U2*h - 2842*h^2,
            ....:      U4*u5 + U3*u5 + U2*u5 + U3*u4 + U2*u4 + U2*u3 - 3791*u5*h
            ....:        - 3791*u4*h - 3791*u3*h - 3791*u2*h - 2842*h^2,
            ....:      U5*u5 - h^2, U4*U2*u3 + U5*U3*u2 + U4*U3*u2 + U3^2*u2 - 3791*U5*U3*h
            ....:        - 3791*U4*U3*h - 3791*U3^2*h - 3791*U5*U2*h- 3791*U4*U2*h + U3*U2*h
            ....:        - 3791*U2^2*h - 3791*U4*u3*h - 3791*U4*u2*h - 3791*U3*u2*h
            ....:        - 2843*U5*h^2 + 1897*U4*h^2 - 946*U3*h^2 - 947*U2*h^2 + 2370*h^3,
            ....:      U3*u5*u4 + U2*u5*u4 + U3*u4^2 + U2*u4^2 + U2*u4*u3 - 3791*u5*u4*h
            ....:        - 3791*u4^2*h - 3791*u4*u3*h - 3791*u4*u2*h + u5*h^2 - 2842*u4*h^2,
            ....:      U2*u5*u4*u3 + U2*u4^2*u3 + U2*u4*u3^2 - 3791*u5*u4*u3*h
            ....:        - 3791*u4^2*u3*h - 3791*u4*u3^2*h - 3791*u4*u3*u2*h + u5*u4*h^2
            ....:        + u4^2*h^2 + u5*u3*h^2 - 2842*u4*u3*h^2,
            ....:      U5^2*U4*u3 + U5*U4^2*u3 + U5^2*U4*u2 + U5*U4^2*u2 + U5^2*U3*u2
            ....:        + 2*U5*U4*U3*u2 + U5*U3^2*u2 - 3791*U5^2*U4*h - 3791*U5*U4^2*h
            ....:        - 3791*U5^2*U3*h + U5*U4*U3*h - 3791*U5*U3^2*h - 3791*U5^2*U2*h
            ....:        + U5*U4*U2*h+ U5*U3*U2*h - 3791*U5*U2^2*h - 3791*U5*U3*u2*h
            ....:        - 2842*U5^2*h^2 + 1897*U5*U4*h^2 - U4^2*h^2 - 947*U5*U3*h^2
            ....:        - U4*U3*h^2 - 948*U5*U2*h^2 - U4*U2*h^2 - 1422*U5*h^3 + 3791*U4*h^3,
            ....:      u5*u4*u3*u2*h + u4^2*u3*u2*h + u4*u3^2*u2*h + u4*u3*u2^2*h
            ....:        + 2*u5*u4*u3*h^2 + 2*u4^2*u3*h^2 + 2*u4*u3^2*h^2 + 2*u5*u4*u2*h^2
            ....:        + 2*u4^2*u2*h^2 + 2*u5*u3*u2*h^2 + 1899*u4*u3*u2*h^2,
            ....:      U5^2*U4*U3*u2 + U5*U4^2*U3*u2 + U5*U4*U3^2*u2 - 3791*U5^2*U4*U3*h
            ....:        - 3791*U5*U4^2*U3*h - 3791*U5*U4*U3^2*h - 3791*U5*U4*U3*U2*h
            ....:        + 3791*U5*U4*U3*u2*h + U5^2*U4*h^2 + U5*U4^2*h^2 + U5^2*U3*h^2
            ....:        - U4^2*U3*h^2 - U5*U3^2*h^2 - U4*U3^2*h^2 - U5*U4*U2*h^2 - U5*U3*U2*h^2
            ....:        - U4*U3*U2*h^2 + 3791*U5*U4*h^3 + 3791*U5*U3*h^3 + 3791*U4*U3*h^3,
            ....:      u4^2*u3*u2*h^2 + 1515*u5*u3^2*u2*h^2 + u4*u3^2*u2*h^2
            ....:        + 1515*u5*u4*u2^2*h^2 + 1515*u5*u3*u2^2*h^2 + u4*u3*u2^2*h^2
            ....:        + 1521*u5*u4*u3*h^3 - 3028*u4^2*u3*h^3 - 3028*u4*u3^2*h^3
            ....:        + 1521*u5*u4*u2*h^3 - 3028*u4^2*u2*h^3 + 1521*u5*u3*u2*h^3
            ....:        + 3420*u4*u3*u2*h^3,
            ....:      U5^2*U4*U3*U2*h + U5*U4^2*U3*U2*h + U5*U4*U3^2*U2*h + U5*U4*U3*U2^2*h
            ....:        + 2*U5^2*U4*U3*h^2 + 2*U5*U4^2*U3*h^2 + 2*U5*U4*U3^2*h^2
            ....:        + 2*U5^2*U4*U2*h^2 + 2*U5*U4^2*U2*h^2 + 2*U5^2*U3*U2*h^2
            ....:        - 2*U4^2*U3*U2*h^2 - 2*U5*U3^2*U2*h^2 - 2*U4*U3^2*U2*h^2
            ....:        - 2*U5*U4*U2^2*h^2 - 2*U5*U3*U2^2*h^2 - 2*U4*U3*U2^2*h^2
            ....:        - U5*U4*U3*h^3 - U5*U4*U2*h^3 - U5*U3*U2*h^3 - U4*U3*U2*h^3]

            sage: Ideal(l).basis_is_groebner()                                          # needs sage.rings.finite_rings
            False
            sage: gb = Ideal(l).groebner_basis()                                        # needs sage.rings.finite_rings
            sage: Ideal(gb).basis_is_groebner()                                         # needs sage.rings.finite_rings
            True

        .. NOTE::

            From the Singular Manual for the reduce function we use in
            this method: 'The result may have no meaning if the second
            argument (``self``) is not a standard basis'. I (malb)
            believe this refers to the mathematical fact that the
            results may have no meaning if ``self`` is no standard basis,
            i.e., Singular doesn't 'add' any additional 'nonsense' to
            the result. So we may actually use reduce to determine if
            ``self`` is a Groebner basis.

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: R.<a,b,c,d,e,f,g,h,i,j> = QQbar[]
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: I.basis_is_groebner()
            False
            sage: I2 = Ideal(I.groebner_basis())
            sage: I2.basis_is_groebner()
            True
        """
        from sage.libs.singular.function_factory import ff
        from sage.libs.singular.option import opt_verb_ctx
        from sage.matrix.constructor import matrix
        sing_reduce = ff.reduce
        syz = ff.syz

        R = self.ring()
        if not R.base_ring().is_field():
            raise ValueError("Coefficient ring must be a field for function 'basis_is_groebner'.")

        try:
            F = matrix(R, 1, self.ngens(), self.gens())
            LTF = matrix(R, 1, self.ngens(), [f.lt() for f in self.gens()])

            with opt_verb_ctx(notWarnSB=True):
                M = F * matrix(R, syz(LTF)).transpose()
                M.set_immutable()
                M = sing_reduce(M, self)

            if any(M):
                return False
            return True
        except TypeError:
            if singular is None:
                singular = singular_default
            R._singular_().set_ring()
            F = singular( tuple(self.gens()), "module" )
            LTF = singular( [f.lt() for f in self.gens()] , "module" )

            M = (F * LTF.syz()).reduce(self._singular_())

            for i in range(M.ncols()):
                if int(singular.eval("%s[1,%s+1]!=0" % (M.name(),i))):
                    return False
            self._singular_().attrib('isSB',1)
        return True

    @require_field
    @handle_AA_and_QQbar
    @singular_gb_standard_options
    @libsingular_gb_standard_options
    def transformed_basis(self, algorithm='gwalk', other_ring=None, singular=None):
        """
        Return a lex or ``other_ring`` Groebner Basis for this ideal.

        INPUT:

        - ``algorithm`` -- see below for options

        - ``other_ring`` -- only valid for ``algorithm='fglm'``; if
          provided, conversion will be performed to this
          ring. Otherwise a lex Groebner basis will be returned.


        ALGORITHMS:

        - ``'fglm'`` -- FGLM algorithm. The input ideal must be given with a reduced
          Groebner Basis of a zero-dimensional ideal

        - ``'gwalk'`` -- Groebner Walk algorithm (*default*)

        - ``'awalk1'`` -- 'first alternative' algorithm

        - ``'awalk2'`` -- 'second alternative' algorithm

        - ``'twalk'`` -- Tran algorithm

        - ``'fwalk'`` -- Fractal Walk algorithm

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = Ideal([y^3+x^2,x^2*y+x^2, x^3-x^2, z^4-x^2-y])
            sage: I = Ideal(I.groebner_basis())
            sage: S.<z,x,y> = PolynomialRing(QQ,3,order='lex')
            sage: J = Ideal(I.transformed_basis('fglm',S))
            sage: J
            Ideal (z^4 + y^3 - y, x^2 + y^3, x*y^3 - y^3, y^4 + y^3)
             of Multivariate Polynomial Ring in z, x, y over Rational Field

        ::

            sage: R.<z,y,x> = PolynomialRing(GF(32003), 3, order='lex')                 # needs sage.rings.finite_rings
            sage: I = Ideal([y^3 + x*y*z + y^2*z + x*z^3, 3 + x*y + x^2*y + y^2*z])     # needs sage.rings.finite_rings
            sage: I.transformed_basis('gwalk')                                          # needs sage.rings.finite_rings
            [z*y^2 + y*x^2 + y*x + 3,
             z*x + 8297*y^8*x^2 + 8297*y^8*x + 3556*y^7 - 8297*y^6*x^4 + 15409*y^6*x^3
               - 8297*y^6*x^2 - 8297*y^5*x^5 + 15409*y^5*x^4 - 8297*y^5*x^3 + 3556*y^5*x^2
               + 3556*y^5*x + 3556*y^4*x^3 + 3556*y^4*x^2 - 10668*y^4 - 10668*y^3*x
               - 8297*y^2*x^9 - 1185*y^2*x^8 + 14224*y^2*x^7 - 1185*y^2*x^6 - 8297*y^2*x^5
               - 14223*y*x^7 - 10666*y*x^6 - 10666*y*x^5 - 14223*y*x^4 + x^5 + 2*x^4 + x^3,
             y^9 - y^7*x^2 - y^7*x - y^6*x^3 - y^6*x^2 - 3*y^6 - 3*y^5*x - y^3*x^7
               - 3*y^3*x^6 - 3*y^3*x^5 - y^3*x^4 - 9*y^2*x^5 - 18*y^2*x^4 - 9*y^2*x^3
               - 27*y*x^3 - 27*y*x^2 - 27*x]

        ALGORITHM:

        Uses Singular.

        TESTS:

        Check that this method works over QQbar (:issue:`25351`).  We are not currently
        able to specify other_ring, due to the limitations of @handle_AA_and_QQbar::

            sage: # needs sage.rings.number_field
            sage: R.<x,y,z> = QQbar[]
            sage: I = Ideal([y^3 + x^2, x^2*y + x^2, x^3 - x^2, z^4 - x^2 - y])
            sage: I = Ideal(I.groebner_basis())
            sage: S.<z,x,y> = PolynomialRing(QQbar, 3, order='lex')
            sage: J = Ideal(I.transformed_basis('fglm', other_ring=S))  # known bug
            sage: J                                                     # known bug
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence
        R = self.ring()

        if self.basis_is_groebner():
            I = R.ideal(self.interreduced_basis())
        else:
            I = R.ideal(self.groebner_basis())

        if algorithm in ("gwalk","awalk1","awalk2","twalk","fwalk"):
            from sage.libs.singular.function import lib, singular_function
            lib("grwalk.lib")
            gb = singular_function(algorithm)(I)
            return PolynomialSequence(R, sorted(gb,reverse=True), immutable=True)

        elif algorithm == "fglm":
            # new ring
            if other_ring is None:
                nR = R.change_ring(order='lex')
            else:
                nR = other_ring
            if singular is None:
                singular = singular_default
            Rs = singular(R)
            Is = singular(I)
            Is.attrib('isSB',1)
            singular(nR).set_ring()
            nIs = singular.fglm(Rs,Is)

            return PolynomialSequence(nR, sorted([nR(f) for f in nIs],reverse=True), immutable=True)

        else:
            raise TypeError("cannot convert basis with given algorithm")

    @handle_AA_and_QQbar
    def elimination_ideal(self, variables, algorithm=None, *args, **kwds):
        r"""
        Return the elimination ideal of this ideal with respect to the
        variables given in ``variables``.

        INPUT:

        - ``variables`` -- list or tuple of variables in ``self.ring()``

        - ``algorithm`` -- determines the algorithm to use, see below
          for available algorithms

        ALGORITHMS:

        - ``'libsingular:eliminate'`` -- libSingular's ``eliminate`` command
          (default)

        - ``'giac:eliminate'`` -- Giac's ``eliminate`` command (if available)

        If only a system is given - e.g. 'giac' - the default algorithm is
        chosen for that system.

        EXAMPLES::

            sage: R.<x,y,t,s,z> = PolynomialRing(QQ,5)
            sage: I = R * [x - t, y - t^2, z - t^3, s - x + y^3]
            sage: J = I.elimination_ideal([t, s]); J
            Ideal (y^2 - x*z, x*y - z, x^2 - y)
             of Multivariate Polynomial Ring in x, y, t, s, z over Rational Field

        You can use Giac to compute the elimination ideal::

            sage: # needs sage.libs.giac
            sage: print("possible output from giac", flush=True); I.elimination_ideal([t, s], algorithm='giac') == J
            possible output...
            True

        The list of available Giac options is provided at
        :func:`sage.libs.giac.groebner_basis`.

        ALGORITHM:

        Uses Singular, or Giac (if available).

        .. NOTE::

            Requires computation of a Groebner basis, which can be a very
            expensive operation.

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field sage.libs.giac
            sage: R.<x,y,t,s,z> = QQbar[]
            sage: I = R * [x - t, y - t^2, z - t^3, s - x + y^3]
            sage: J = I.elimination_ideal([t, s]); J
            Ideal (y^2 - x*z, x*y - z, x^2 - y) of Multivariate
            Polynomial Ring in x, y, t, s, z over Algebraic Field
            sage: print("possible output from giac", flush=True); I.elimination_ideal([t, s], algorithm='giac') == J
            possible output...
            True

        Check that the passed variables are actually variables of the ring
        (:issue:`31414`)::

            sage: R.<x,y,z> = QQ[]
            sage: I = R.ideal(x-y, z)
            sage: I.elimination_ideal([x, R(0)])
            Traceback (most recent call last):
            ...
            ValueError: not a ring variable: 0
        """
        if not isinstance(variables, (list, tuple)):
            variables = (variables,)
        gens = self.ring().gens()
        for v in variables:
            if v not in gens:
                raise ValueError("not a ring variable: %s" % v)

        if (algorithm is None or algorithm.lower() == 'libsingular'
                or algorithm == 'libsingular:eliminate'):
            return self._elimination_ideal_libsingular(variables)

        elif algorithm.lower() == 'giac' or algorithm == 'giac:eliminate':
            from sage.libs.giac import groebner_basis as groebner_basis_libgiac
            return groebner_basis_libgiac(
                    self, elim_variables=variables, *args, **kwds).ideal()

        else:
            raise NameError("Algorithm '%s' unknown." % algorithm)

    @libsingular_gb_standard_options
    def _elimination_ideal_libsingular(self, variables):
        r"""
        Compute the elimination ideal using libsingular.

        EXAMPLES::

            sage: R.<x,y,t,s,z> = PolynomialRing(QQ,5)
            sage: I = R * [x-t,y-t^2,z-t^3,s-x+y^3]
            sage: I.elimination_ideal([t,s], "libsingular")  # indirect doctest
            Ideal (y^2 - x*z, x*y - z, x^2 - y) of Multivariate
            Polynomial Ring in x, y, t, s, z over Rational Field
        """
        from sage.libs.singular.function_factory import ff
        eliminate = ff.eliminate

        R = self.ring()
        Is = MPolynomialIdeal(R,self.groebner_basis())
        return MPolynomialIdeal(R, eliminate(Is, prod(variables)) )

    @handle_AA_and_QQbar
    @libsingular_gb_standard_options
    def quotient(self, J):
        r"""
        Given ideals `I` = ``self`` and `J` in the same polynomial
        ring `P`, return the ideal quotient of `I` by `J` consisting
        of the polynomials `a` of `P` such that `\{aJ \subset I\}`.

        This is also referred to as the colon ideal
        (`I`:`J`).

        INPUT:

        - ``J`` -- multivariate polynomial ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y,z> = PolynomialRing(GF(181), 3)
            sage: I = Ideal([x^2 + x*y*z, y^2 - z^3*y, z^3 + y^5*x*z])
            sage: J = Ideal([x])
            sage: Q = I.quotient(J)
            sage: y*z + x in I
            False
            sage: x in J
            True
            sage: x * (y*z + x) in I
            True

        TESTS:

        This example checks :issue:`16301`::

            sage: R.<x,y,z> = ZZ[]
            sage: I = Ideal(R(2), x*y, x*z + x)
            sage: eD = Ideal(x, z^2-1)
            sage: I.quotient(eD).gens()
            [2, x*z + x, x*y]

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: R.<x,y,z> = QQbar[]
            sage: I = ideal(x, z)
            sage: J = ideal(R(1))
            sage: I.quotient(J)
            Ideal (z, x) of Multivariate Polynomial Ring in x, y, z over Algebraic Field

        Check that :issue:`12803` is fixed::

            sage: R.<xe,xv> = ZZ[]
            sage: J = ideal(4*xv^3 + 3*xv^2, 3*xe*xv^2 + xe - 2*xv)
            sage: I  = ideal(-3, -3*xv - 1, -3)
            sage: I2 = ideal(-3, -3*xv - 1)
            sage: I == I2
            True
            sage: Q1 = J.quotient(I)
            sage: Q2 = J.quotient(I2)
            sage: Q1 == Q2
            True
        """
        R = self.ring()

        if not isinstance(J, MPolynomialIdeal):
            raise TypeError("J needs to be a multivariate polynomial ideal")

        if R is not J.ring() and not R == J.ring():
            raise TypeError("base rings do not match")

        from sage.libs.singular.function_factory import ff
        quotient = ff.quotient
        return R.ideal(quotient(self, J))

    @handle_AA_and_QQbar
    def saturation(self, other):
        r"""
        Return the saturation (and saturation exponent) of the ideal ``self``
        with respect to the ideal ``other``.

        INPUT:

        - ``other`` -- another ideal in the same ring

        OUTPUT: a pair (ideal, integer)

        EXAMPLES::

            sage: R.<x, y, z> = QQ[]
            sage: I = R.ideal(x^5*z^3, x*y*z, y*z^4)
            sage: J = R.ideal(z)
            sage: I.saturation(J)
            (Ideal (y, x^5) of Multivariate Polynomial Ring in x, y, z over Rational Field, 4)

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: R.<x, y, z> = QQbar[]
            sage: I = R.ideal(x^5*z^3, x*y*z, y*z^4)
            sage: J = R.ideal(z)
            sage: I.saturation(other = J)
            (Ideal (y, x^5) of Multivariate Polynomial Ring in x, y, z over Algebraic Field, 4)
        """
        from sage.libs.singular.function_factory import ff

        # function renamed in singular > 4.3.2p4, see issue #35980
        try:
            sat = ff.elim__lib.sat_with_exp
        except NameError:
            sat = ff.elim__lib.sat
        R = self.ring()
        ideal, expo = sat(self, other)
        return (R.ideal(ideal), ZZ(expo))

    @require_field
    def variety(self, ring=None, *, algorithm='triangular_decomposition', proof=True):
        r"""
        Return the variety of this ideal.

        Given a zero-dimensional ideal `I` (= ``self``) of a
        polynomial ring `P` whose order is lexicographic, return the
        variety of `I` as a list of dictionaries with ``(variable, value)``
        pairs. By default, the variety of the ideal over its
        coefficient field `K` is returned; ``ring`` can be specified
        to find the variety over a different ring.

        These dictionaries have cardinality equal to the number of
        variables in `P` and represent assignments of values to these
        variables such that all polynomials in `I` vanish.

        If ``ring`` is specified, then a triangular decomposition of
        ``self`` is found over the original coefficient field `K`;
        then the triangular systems are solved using root-finding over
        ``ring``. This is particularly useful when `K` is ``QQ`` (to
        allow fast symbolic computation of the triangular
        decomposition) and ``ring`` is ``RR``, ``AA``, ``CC``, or
        ``QQbar`` (to compute the whole real or complex variety of the
        ideal).

        Note that with ``ring=RR`` or ``CC``, computation is done
        numerically and potentially inaccurately; in particular, the
        number of points in the real variety may be miscomputed. With
        ``ring=AA`` or ``QQbar``, computation is done exactly
        (which may be much slower, of course).

        INPUT:

        - ``ring`` -- return roots in the ``ring`` instead of the base
          ring of this ideal (default: ``None``)
        - ``algorithm`` -- algorithm or implementation to use; see below for
          supported values
        - ``proof`` -- return a provably correct result (default: ``True``)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<w> = GF(27)  # this example is from the MAGMA handbook
            sage: P.<x, y> = PolynomialRing(K, 2, order='lex')
            sage: I = Ideal([x^8 + y + 2, y^6 + x*y^5 + x^2])
            sage: I = Ideal(I.groebner_basis()); I
            Ideal (x - y^47 - y^45 + y^44 - y^43 + y^41 - y^39 - y^38 - y^37 - y^36
                     + y^35 - y^34 - y^33 + y^32 - y^31 + y^30 + y^28 + y^27 + y^26
                     + y^25 - y^23 + y^22 + y^21 - y^19 - y^18 - y^16 + y^15 + y^13
                     + y^12 - y^10 + y^9 + y^8 + y^7 - y^6 + y^4 + y^3 + y^2 + y - 1,
                   y^48 + y^41 - y^40 + y^37 - y^36 - y^33 + y^32 - y^29 + y^28
                     - y^25 + y^24 + y^2 + y + 1)
            of Multivariate Polynomial Ring in x, y over Finite Field in w of size 3^3
            sage: V = I.variety();
            sage: sorted(V, key=str)
            [{y: w^2 + 2*w, x: 2*w + 2}, {y: w^2 + 2, x: 2*w}, {y: w^2 + w, x: 2*w + 1}]
            sage: [f.subs(v)                      # check that all polynomials vanish
            ....:  for f in I.gens() for v in V]
            [0, 0, 0, 0, 0, 0]
            sage: [I.subs(v).is_zero() for v in V]  # same test, but nicer syntax
            [True, True, True]

        However, we only account for solutions in the ground field and not
        in the algebraic closure::

            sage: I.vector_space_dimension()                                            # needs sage.rings.finite_rings
            48

        Here we compute the points of intersection of a hyperbola and a
        circle, in several fields::

            sage: K.<x, y> = PolynomialRing(QQ, 2, order='lex')
            sage: I = Ideal([x*y - 1, (x-2)^2 + (y-1)^2 - 1])
            sage: I = Ideal(I.groebner_basis()); I
            Ideal (x + y^3 - 2*y^2 + 4*y - 4, y^4 - 2*y^3 + 4*y^2 - 4*y + 1)
             of Multivariate Polynomial Ring in x, y over Rational Field

        These two curves have one rational intersection::

            sage: I.variety()
            [{y: 1, x: 1}]

        There are two real intersections::

            sage: sorted(I.variety(ring=RR), key=str)
            [{y: 0.361103080528647, x: 2.76929235423863},
             {y: 1.00000000000000, x: 1.00000000000000}]
            sage: I.variety(ring=AA)                                                    # needs sage.rings.number_field
            [{y: 1, x: 1},
             {y: 0.3611030805286474?, x: 2.769292354238632?}]

        and a total of four intersections::

            sage: sorted(I.variety(ring=CC), key=str)
            [{y: 0.31944845973567... + 1.6331702409152...*I,
              x: 0.11535382288068... - 0.58974280502220...*I},
             {y: 0.31944845973567... - 1.6331702409152...*I,
              x: 0.11535382288068... + 0.58974280502220...*I},
             {y: 0.36110308052864..., x: 2.7692923542386...},
             {y: 1.00000000000000, x: 1.00000000000000}]
            sage: sorted(I.variety(ring=QQbar), key=str)                                # needs sage.rings.number_field
            [{y: 0.3194484597356763? + 1.633170240915238?*I,
              x: 0.11535382288068429? - 0.5897428050222055?*I},
             {y: 0.3194484597356763? - 1.633170240915238?*I,
              x: 0.11535382288068429? + 0.5897428050222055?*I},
             {y: 0.3611030805286474?, x: 2.769292354238632?},
             {y: 1, x: 1}]

        We can also use the :ref:`optional package msolve <spkg_msolve>`
        to compute the variety.
        See :mod:`~sage.rings.polynomial.msolve` for more information. ::

            sage: I.variety(RBF, algorithm='msolve', proof=False)   # optional - msolve
            [{x: [2.76929235423863 +/- 2.08e-15], y: [0.361103080528647 +/- 4.53e-16]},
             {x: 1.000000000000000, y: 1.000000000000000}]

        Computation over floating point numbers may compute only a partial solution,
        or even none at all. Notice that x values are missing from the following variety::

            sage: R.<x,y> = CC[]
            sage: I = ideal([x^2+y^2-1,x*y-1])
            sage: sorted(I.variety(), key=str)
            verbose 0 (...: multi_polynomial_ideal.py, variety) Warning: computations in the complex field are inexact; variety may be computed partially or incorrectly.
            verbose 0 (...: multi_polynomial_ideal.py, variety) Warning: falling back to very slow toy implementation.
            [{y: -0.86602540378443... + 0.500000000000000*I},
             {y: -0.86602540378443... - 0.500000000000000*I},
             {y: 0.86602540378443... + 0.500000000000000*I},
             {y: 0.86602540378443... - 0.500000000000000*I}]

        This is due to precision error,
        which causes the computation of an intermediate Groebner basis to fail.

        If the ground field's characteristic is too large for
        Singular, we resort to a toy implementation::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = PolynomialRing(GF(2147483659^3), order='lex')
            sage: I = ideal([x^3 - 2*y^2, 3*x + y^4])
            sage: I.variety()
            verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            verbose 0 (...: multi_polynomial_ideal.py, variety) Warning: falling back to very slow toy implementation.
            [{y: 0, x: 0}]

        The dictionary expressing the variety will be indexed by generators
        of the polynomial ring after changing to the target field.
        But the mapping will also accept generators of the original ring,
        or even generator names as strings, when provided as keys::

            sage: # needs sage.rings.number_field
            sage: K.<x,y> = QQ[]
            sage: I = ideal([x^2 + 2*y - 5, x + y + 3])
            sage: v = I.variety(AA)[0]; v[x], v[y]
            (4.464101615137755?, -7.464101615137755?)
            sage: list(v)[0].parent()
            Multivariate Polynomial Ring in x, y over Algebraic Real Field
            sage: v[x]
            4.464101615137755?
            sage: v["y"]
            -7.464101615137755?

        ``msolve`` also works over finite fields::

            sage: R.<x, y> = PolynomialRing(GF(536870909), 2, order='lex')              # needs sage.rings.finite_rings
            sage: I = Ideal([x^2 - 1, y^2 - 1])                                         # needs sage.rings.finite_rings
            sage: sorted(I.variety(algorithm='msolve',          # optional - msolve, needs sage.rings.finite_rings
            ....:                  proof=False),
            ....:        key=str)
            [{x: 1, y: 1},
             {x: 1, y: 536870908},
             {x: 536870908, y: 1},
             {x: 536870908, y: 536870908}]

        but may fail in small characteristic, especially with ideals of high
        degree with respect to the characteristic::

            sage: R.<x, y> = PolynomialRing(GF(3), 2, order='lex')
            sage: I = Ideal([x^2 - 1, y^2 - 1])
            sage: I.variety(algorithm='msolve', proof=False)    # optional - msolve
            Traceback (most recent call last):
            ...
            NotImplementedError: characteristic 3 too small

        ALGORITHM:

        - With ``algorithm`` = ``'triangular_decomposition'`` (default),
          uses triangular decomposition, via Singular if possible, falling back
          on a toy implementation otherwise.

        - With ``algorithm`` = ``'msolve'``, uses the
          :ref:`optional package msolve <spkg_msolve>`.
          Note that msolve uses heuristics and therefore
          requires setting the ``proof`` flag to ``False``. See
          :mod:`~sage.rings.polynomial.msolve` for more information.
        """
        if algorithm == "triangular_decomposition":
            return self._variety_triangular_decomposition(ring)
        elif algorithm == "msolve":
            from . import msolve
            return msolve.variety(self, ring, proof=proof)
        else:
            raise ValueError(f"unknown algorithm {algorithm!r}")

    def _variety_triangular_decomposition(self, ring):
        r"""
        Compute the variety of this ideal by triangular decomposition.

        The triangular decomposition is computed using Singular when conversion
        of the ideal to Singular is supported, falling back to a toy Python
        implementation otherwise.

        TESTS::

            sage: K.<w> = GF(27)                                                        # needs sage.rings.finite_rings
            sage: P.<x, y> = PolynomialRing(K, 2, order='lex')                          # needs sage.rings.finite_rings
            sage: I = Ideal([ x^8 + y + 2, y^6 + x*y^5 + x^2 ])                         # needs sage.rings.finite_rings

        Testing the robustness of the Singular interface::

            sage: T = I.triangular_decomposition('singular:triangLfak')                 # needs sage.rings.finite_rings
            sage: sorted(I.variety(), key=str)                                          # needs sage.rings.finite_rings
            [{y: w^2 + 2*w, x: 2*w + 2}, {y: w^2 + 2, x: 2*w}, {y: w^2 + w, x: 2*w + 1}]

        Testing that a bug is indeed fixed ::

            sage: R = PolynomialRing(GF(2), 30, ['x%d'%(i+1) for i in range(30)], order='lex')
            sage: R.inject_variables()
            Defining...
            sage: I = Ideal([x1 + 1, x2, x3 + 1, x5*x10 + x10 + x18, x5*x11 + x11, \
                             x5*x18, x6, x7 + 1, x9, x10*x11 + x10 + x18, x10*x18 + x18, \
                             x11*x18, x12, x13, x14, x15, x16 + 1, x17 + x18 + 1, x19, x20, \
                             x21 + 1, x22, x23, x24, x25 + 1, x28 + 1, x29 + 1, x30, x8, \
                             x26, x1^2 + x1, x2^2 + x2, x3^2 + x3, x4^2 + x4, x5^2 + x5, \
                             x6^2 + x6, x7^2 + x7, x8^2 + x8, x9^2 + x9, x10^2 + x10, \
                             x11^2 + x11, x12^2 + x12, x13^2 + x13, x14^2 + x14, x15^2 + x15, \
                             x16^2 + x16, x17^2 + x17, x18^2 + x18, x19^2 + x19, x20^2 + x20, \
                             x21^2 + x21, x22^2 + x22, x23^2 + x23, x24^2 + x24, x25^2 + x25, \
                             x26^2 + x26, x27^2 + x27, x28^2 + x28, x29^2 + x29, x30^2 + x30])
            sage: I.basis_is_groebner()
            True
            sage: sorted("".join(str(V[g]) for g in R.gens()) for V in I.variety())  # long time (6s on sage.math, 2011)
            ['101000100000000110001000100110',
             '101000100000000110001000101110',
             '101000100100000101001000100110',
             '101000100100000101001000101110',
             '101010100000000110001000100110',
             '101010100000000110001000101110',
             '101010100010000110001000100110',
             '101010100010000110001000101110',
             '101010100110000110001000100110',
             '101010100110000110001000101110',
             '101100100000000110001000100110',
             '101100100000000110001000101110',
             '101100100100000101001000100110',
             '101100100100000101001000101110',
             '101110100000000110001000100110',
             '101110100000000110001000101110',
             '101110100010000110001000100110',
             '101110100010000110001000101110',
             '101110100110000110001000100110',
             '101110100110000110001000101110']

        Check that the issue at :issue:`7425` is fixed::

            sage: R.<x, y, z> = QQ[]
            sage: I = R.ideal([x^2-y^3*z, x+y*z])
            sage: I.dimension()
            1
            sage: I.variety()
            Traceback (most recent call last):
            ...
            ValueError: The dimension of the ideal is 1, but it should be 0

        Check that the issue at :issue:`7425` is fixed::

            sage: S.<t> = PolynomialRing(QQ)
            sage: F.<q> = QQ.extension(t^4 + 1)
            sage: R.<x,y> = PolynomialRing(F)
            sage: I = R.ideal(x, y^4 + 1)
            sage: I.variety()
            [...{y: -q^3, x: 0}...]

        Check that computing the variety of the ``[1]`` ideal is allowed (:issue:`13977`)::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal(1)
            sage: I.variety()
            []

        Check that the issue at :issue:`16485` is fixed::

            sage: R.<a,b,c> = PolynomialRing(QQ, order='lex')
            sage: I = R.ideal(c^2 - 2, b - c, a)
            sage: I.variety(QQbar)                                                      # needs sage.rings.number_field
            [...a: 0...]

        An early version of :issue:`25351` broke this method by adding the
        @handle_AA_and_QQbar decorator to it.  To check for this bug, verify
        that this circle and this hyperbola have two real intersections and
        two more complex ones::

            sage: # needs sage.rings.number_field
            sage: K.<x, y> = PolynomialRing(AA)
            sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
            sage: len(I.variety())
            2
            sage: K.<x, y> = PolynomialRing(QQbar)
            sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
            sage: len(I.variety())
            4
        """

        def _variety(T, V, v=None):
            """
            Return variety ``V`` for one triangular set of
            polynomials ``T``.
            """
            if v is None:
                v = {}
            found = False
            for f in T:
                if f.is_univariate() and not f.is_constant():
                    T.remove(f)
                    found = True
                    break

            if found is False:
                V.append(v)
                return V

            variable = f.variable(0)
            roots = f.univariate_polynomial().roots(ring=ring, multiplicities=False)

            for root in roots:
                vbar = v.copy()
                vbar[variable] = root
                Tbar = [ f.subs({variable:root}) for f in T ]
                _variety(Tbar,V,vbar)

            return V

        d = self.dimension()
        if d > 0:
            raise ValueError("The dimension of the ideal is %s, but it should be 0" % d)
        if d == -1:
            return []

        if isinstance(self.base_ring(), sage.rings.abc.ComplexField):
            verbose("Warning: computations in the complex field are inexact; variety may be computed partially or incorrectly.", level=0, caller_name='variety')
        P = self.ring()
        if ring is not None:
            P = P.change_ring(ring)
        try:
            TI = self.triangular_decomposition('singular:triangLfak')
            T = [list(each.gens()) for each in TI]
        except TypeError:  # conversion to Singular not supported
            if self.ring().term_order().is_global():
                verbose("Warning: falling back to very slow toy implementation.", level=0, caller_name='variety')
                T = toy_variety.triangular_factorization(self.groebner_basis())
            else:
                raise TypeError("Local/unknown orderings not supported by 'toy_buchberger' implementation.")

        from sage.misc.converting_dict import KeyConvertingDict
        V = []
        for t in T:
            V.extend(KeyConvertingDict(P, v)
                     for v in _variety([P(f) for f in t], []))
        return V

    @require_field
    @handle_AA_and_QQbar
    def hilbert_polynomial(self, algorithm='sage'):
        r"""
        Return the Hilbert polynomial of this ideal.

        INPUT:

        - ``algorithm`` -- (default: ``'sage'``) must be either ``'sage'``
          or ``'singular'``

        Let `I` (which is ``self``) be a homogeneous ideal and
        `R = \bigoplus_d R_d` (which is ``self.ring()``) be a graded
        commutative algebra over a field `K`. The *Hilbert polynomial*
        is the unique polynomial `HP(t)` with rational coefficients
        such that `HP(d) = \dim_K R_d` for all but finitely many
        positive integers `d`.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I.hilbert_polynomial()                                                # needs sage.libs.flint
            5*t - 5

        Of course, the Hilbert polynomial of a zero-dimensional ideal
        is zero::

            sage: J0 = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5,
            ....:             y^3 - 2*x*z^2 + x*y, x^4 + x*y - y*z^2])
            sage: J = P*[m.lm() for m in J0.groebner_basis()]
            sage: J.dimension()
            0
            sage: J.hilbert_polynomial()                                                # needs sage.libs.flint
            0

        It is possible to request a computation using the Singular library::

            sage: I.hilbert_polynomial(algorithm='singular') == I.hilbert_polynomial()  # needs sage.libs.flint
            True
            sage: J.hilbert_polynomial(algorithm='singular') == J.hilbert_polynomial()  # needs sage.libs.flint
            True

        Here is a bigger examples::

            sage: n = 4; m = 11; P = PolynomialRing(QQ, n * m, "x"); x = P.gens()
            sage: M = Matrix(n, x)
            sage: Minors = P.ideal(M.minors(2))
            sage: hp = Minors.hilbert_polynomial(); hp                                  # needs sage.libs.flint
            1/21772800*t^13 + 61/21772800*t^12 + 1661/21772800*t^11
             + 26681/21772800*t^10 + 93841/7257600*t^9 + 685421/7257600*t^8
             + 1524809/3110400*t^7 + 39780323/21772800*t^6 + 6638071/1360800*t^5
             + 12509761/1360800*t^4 + 2689031/226800*t^3 + 1494509/151200*t^2
             + 12001/2520*t + 1

        Because Singular uses 32-bit integers, the above example would fail
        with Singular. We don't test it here, as it has a side-effect on
        other tests that is not understood yet (see :issue:`26300`)::

            sage: Minors.hilbert_polynomial(algorithm='singular')       # not tested
            Traceback (most recent call last):
            ...
            RuntimeError: error in Singular function call 'hilbPoly':
            int overflow in hilb 1
            error occurred in or before poly.lib::hilbPoly line 58: `   intvec v=hilb(I,2);`
            expected intvec-expression. type 'help intvec;'

        Note that in this example, the Hilbert polynomial gives the
        coefficients of the Hilbert-Poincaré series in all degrees::

            sage: P = PowerSeriesRing(QQ, 't', default_prec=50)
            sage: hs = Minors.hilbert_series()                                          # needs sage.libs.flint
            sage: list(P(hs.numerator()) / P(hs.denominator())) == [hp(t=k)             # needs sage.libs.flint
            ....:                                                   for k in range(50)]
            True

        TESTS:

        Check that :issue:`27483` and :issue:`28110` are fixed::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3, x*y^2, y^4, x^2*y*z, y^3*z, x^2*z^2, x*y*z^2, x*z^3])
            sage: I.hilbert_polynomial(algorithm='singular')
            3
            sage: I.hilbert_polynomial()                                                # needs sage.libs.flint
            3

        Check that this method works over ``QQbar`` (:issue:`25351`)::

            sage: P.<x,y,z> = QQbar[]                                                   # needs sage.rings.number_field
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])                    # needs sage.rings.number_field
            sage: I.hilbert_polynomial()                                                # needs sage.rings.number_field
            5*t - 5

        Check for :issue:`33597`::

            sage: R.<X, Y, Z> = QQ[]
            sage: I = R.ideal([X^2*Y^3, X*Z])
            sage: I.hilbert_polynomial()                                                # needs sage.libs.flint
            t + 5
        """
        if not self.is_homogeneous():
            raise TypeError("ideal must be homogeneous")
        if algorithm == 'sage':
            from sage.misc.misc_c import prod
            hilbert_poincare = self.hilbert_series()
            denom = hilbert_poincare.denominator()
            if denom.degree() == 0:
                return denom.parent().zero()
            t = denom.parent().gen()
            s = denom.valuation(t - 1)
            numerator = hilbert_poincare.numerator()
            # we assume the denominator of the Hilbert series is of
            # the form (1 - t)^s, need to scale numerator
            scalar = ~(denom[0] * (s - 1).factorial())
            st = s - 1 + t
            out = scalar * sum(c * prod(st - n - nu for nu in range(s - 1))
                               for n, c in enumerate(numerator))
            return t.parent().zero() + out
        if algorithm == 'singular':
            from sage.libs.singular.function_factory import ff
            hilbPoly = ff.polylib__lib.hilbPoly

            hp = hilbPoly(self)
            t = ZZ['t'].gen()
            fp = ZZ(len(hp) - 1).factorial()
            return sum(ZZ(coeff) * t**i for i, coeff in enumerate(hp)) / fp
        raise ValueError("'algorithm' must be 'sage' or 'singular'")

    @require_field
    @handle_AA_and_QQbar
    def hilbert_series(self, grading=None, algorithm='sage'):
        r"""
        Return the Hilbert series of this ideal.

        INPUT:

        - ``grading`` -- (optional) a list or tuple of integers
        - ``algorithm`` -- (default: ``'sage'``) must be either ``'sage'``
          or ``'singular'``

        Let `I` (which is ``self``) be a homogeneous ideal and
        `R = \bigoplus_d R_d` (which is ``self.ring()``) be a
        graded commutative algebra over a field `K`. Then the
        *Hilbert function* is defined as `H(d) = \dim_K R_d` and
        the *Hilbert series* of `I` is defined as the formal power
        series `HS(t) = \sum_{d=0}^{\infty} H(d) t^d`.

        This power series can be expressed as
        `HS(t) = Q(t) / (1-t)^n` where `Q(t)` is a polynomial
        over `\ZZ` and `n` the number of variables in `R`.
        This method returns `Q(t) / (1-t)^n`, normalised so
        that the leading monomial of the numerator is positive.

        An optional ``grading`` can be given, in which case
        the graded (or weighted) Hilbert series is given.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I.hilbert_series()                                                    # needs sage.libs.flint
            (t^4 + t^3 + t^2 + t + 1)/(t^2 - 2*t + 1)
            sage: R.<a,b> = PolynomialRing(QQ)
            sage: J = R.ideal([a^2*b, a*b^2])
            sage: J.hilbert_series()                                                    # needs sage.libs.flint
            (t^3 - t^2 - t - 1)/(t - 1)
            sage: J.hilbert_series(grading=(10,3))                                      # needs sage.libs.flint
            (t^25 + t^24 + t^23 - t^15 - t^14 - t^13 - t^12 - t^11
             - t^10 - t^9 - t^8 - t^7 - t^6 - t^5 - t^4 - t^3 - t^2
             - t - 1)/(t^12 + t^11 + t^10 - t^2 - t - 1)

            sage: K = R.ideal([a^2*b^3, a*b^4 + a^3*b^2])
            sage: K.hilbert_series(grading=[1,2])                                       # needs sage.libs.flint
            (t^11 + t^8 - t^6 - t^5 - t^4 - t^3 - t^2 - t - 1)/(t^2 - 1)
            sage: K.hilbert_series(grading=[2,1])                                       # needs sage.libs.flint
            (2*t^7 - t^6 - t^4 - t^2 - 1)/(t - 1)

        TESTS::

            sage: # needs sage.libs.flint
            sage: I.hilbert_series() == I.hilbert_series(algorithm='singular')
            True
            sage: J.hilbert_series() == J.hilbert_series(algorithm='singular')
            True
            sage: J.hilbert_series(grading=(10,3)) == J.hilbert_series(grading=(10,3), algorithm='singular')
            True
            sage: K.hilbert_series(grading=(1,2)) == K.hilbert_series(grading=(1,2), algorithm='singular')
            True
            sage: K.hilbert_series(grading=(2,1)) == K.hilbert_series(grading=(2,1), algorithm='singular')
            True

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I.hilbert_series(grading=5)                                                                           # needs sage.libs.flint
            Traceback (most recent call last):
            ...
            TypeError: grading must be a list or a tuple of integers

        Check that this method works over QQbar (:issue:`25351`)::

            sage: P.<x,y,z> = QQbar[]                                                                                   # needs sage.rings.number_field
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])                                                    # needs sage.rings.number_field
            sage: I.hilbert_series()                                                                                    # needs sage.rings.number_field
            (t^4 + t^3 + t^2 + t + 1)/(t^2 - 2*t + 1)
        """
        if not self.is_homogeneous():
            raise TypeError("ideal must be homogeneous")

        if algorithm == 'sage':
            from sage.rings.polynomial.hilbert import hilbert_poincare_series

            if grading is not None:
                if not isinstance(grading, (list, tuple)) or any(a not in ZZ for a in grading):
                    raise TypeError("grading must be a list or a tuple of integers")
            gb = MPolynomialIdeal(self.ring(), [mon.lm() for mon in self.groebner_basis()])

            return hilbert_poincare_series(gb, grading)
        elif algorithm == 'singular':
            t = ZZ['t'].gen()
            n = self.ring().ngens()

            if grading is None:
                return self.hilbert_numerator(algorithm='singular') / (1 - t)**n

            # The check that ``grading`` is valid input is done by ``hilbert_numerator()``
            return (self.hilbert_numerator(algorithm='singular', grading=grading)
                    / prod((1 - t**a) for a in grading))
        else:
            raise ValueError("'algorithm' must be one of 'sage' or 'singular'")

    @require_field
    @handle_AA_and_QQbar
    def hilbert_numerator(self, grading=None, algorithm='sage'):
        r"""
        Return the Hilbert numerator of this ideal.

        INPUT:

        - ``grading`` -- (optional) a list or tuple of integers
        - ``algorithm`` -- (default: ``'sage'``) must be either ``'sage'``
          or ``'singular'``

        Let `I` (which is ``self``) be a homogeneous ideal and
        `R = \bigoplus_d R_d` (which is ``self.ring()``) be a
        graded commutative algebra over a field `K`. Then the
        *Hilbert function* is defined as `H(d) = \dim_K R_d` and
        the *Hilbert series* of `I` is defined as the formal power
        series `HS(t) = \sum_{d=0}^{\infty} H(d) t^d`.

        This power series can be expressed as
        `HS(t) = Q(t) / (1-t)^n` where `Q(t)` is a polynomial
        over `\ZZ` and `n` the number of variables in `R`. This
        method returns `Q(t)`, the numerator; hence the name,
        ``hilbert_numerator``. An optional ``grading`` can be given, in
        which case the graded (or weighted) Hilbert numerator is given.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I.hilbert_numerator()                                                 # needs sage.libs.flint
            -t^5 + 1
            sage: R.<a,b> = PolynomialRing(QQ)
            sage: J = R.ideal([a^2*b, a*b^2])
            sage: J.hilbert_numerator()                                                 # needs sage.libs.flint
            t^4 - 2*t^3 + 1
            sage: J.hilbert_numerator(grading=(10,3))                                   # needs sage.libs.flint
            t^26 - t^23 - t^16 + 1

        TESTS::

            sage: I.hilbert_numerator() == I.hilbert_numerator(algorithm='singular')    # needs sage.libs.flint
            True
            sage: J.hilbert_numerator() == J.hilbert_numerator(algorithm='singular')    # needs sage.libs.flint
            True
            sage: J.hilbert_numerator(grading=(10,3)) == J.hilbert_numerator(grading=(10,3), algorithm='singular')      # needs sage.libs.flint
            True

        Check that this method works over QQbar (:issue:`25351`)::

            sage: P.<x,y,z> = QQbar[]                                                   # needs sage.rings.number_field
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])                    # needs sage.rings.number_field
            sage: I.hilbert_numerator()                                                 # needs sage.rings.number_field
            -t^5 + 1

        This example returns a wrong answer in singular < 4.3.2p4 due to an integer overflow::

            sage: n=4; m=11; P = PolynomialRing(QQ, n*m, "x"); x = P.gens(); M = Matrix(n, x)
            sage: I = P.ideal(M.minors(2))
            sage: J = P * [m.lm() for m in I.groebner_basis()]
            sage: J.hilbert_numerator(algorithm='singular') # known bug
            Traceback (most recent call last):
            ....
            RuntimeError: error in Singular function call 'hilb':
            overflow at t^22

        Our two algorithms should always agree; not tested until
        :issue:`33178` is fixed::

            sage: # not tested
            sage: m = ZZ.random_element(2,8)
            sage: nvars = m^2
            sage: R = PolynomialRing(QQ, "x", nvars)
            sage: M = matrix(R, m, R.gens())
            sage: I = R.ideal(M.minors(2))
            sage: n1 = I.hilbert_numerator()
            sage: n2 = I.hilbert_numerator(algorithm='singular')
            sage: n1 == n2
            True
        """
        if not self.is_homogeneous():
            raise TypeError("Ideal must be homogeneous.")

        if algorithm == 'sage':
            from sage.rings.polynomial.hilbert import first_hilbert_series

            if grading is not None:
                if not isinstance(grading, (list, tuple)) or any(a not in ZZ for a in grading):
                    raise TypeError("grading must be a list or a tuple of integers")
            gb = MPolynomialIdeal(self.ring(), [mon.lm() for mon in self.groebner_basis()])

            return first_hilbert_series(gb, grading)
        elif algorithm == 'singular':
            from sage.libs.singular.function_factory import ff
            hilb = ff.hilb

            gb = self.groebner_basis()
            t = ZZ['t'].gen()
            gb = MPolynomialIdeal(self.ring(), gb)
            if grading is not None:
                if not isinstance(grading, (list, tuple)) or any(a not in ZZ for a in grading):
                    raise TypeError("grading must be a list or a tuple of integers")

                hs = hilb(gb, 1, tuple(grading), attributes={gb: {'isSB': 1}})
            else:
                hs = hilb(gb, 1, attributes={gb: {'isSB': 1}})
            return sum(ZZ(hs[i]) * t**i for i in range(len(hs)-1))
        else:
            raise ValueError("'algorithm' must be one of 'sage' or 'singular'")

    @require_field
    def _normal_basis_libsingular(self, degree, weights=None):
        r"""
        Return the normal basis for a given Groebner basis.

        This will use the Groebner basis as computed by
        ``MPolynomialIdeal._groebner_basis_libsingular()``.

        INPUT:

        - ``degree`` -- ``None`` or integer

        - ``weights`` -- tuple of positive integers (default: ``None``); if not
          ``None``, compute the degree with respect to these weights

        OUTPUT:

        If ``degree`` is an integer, only the monomials of the given degree in
        the normal basis.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = R.ideal(x^2-2*x*z+5, x*y^2+y*z+1, 3*y^2-8*x*z)
            sage: I.normal_basis()  # indirect doctest
            [z^2, y*z, x*z, z, x*y, y, x, 1]
            sage: J = R.ideal(x^2-2*x*z+5)
            sage: J.normal_basis(3)  # indirect doctest
            [z^3, y*z^2, x*z^2, y^2*z, x*y*z, y^3, x*y^2]
            sage: [J._normal_basis_libsingular(d, (2, 2, 3)) for d in (0..8)]
            [[1], [], [x, y], [z], [x*y, y^2], [x*z, y*z], [x*y^2, y^3, z^2],
             [x*y*z, y^2*z], [x*y^3, y^4, x*z^2, y*z^2]]

        TESTS:

        Check that zero is not included in trivial results::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = R.ideal(x^2-2*x*z+5, x*y^2+y*z+1, 3*y^2-8*x*z)
            sage: I._normal_basis_libsingular(5)
            []

        Check what happens with weights but no degree (:issue:`34789`)::

            sage: TO = TermOrder('wdegrevlex', (2,))
            sage: R = PolynomialRing(QQ, 1, ['k'], order=TO)
            sage: k = R.gen()
            sage: I = R.ideal([k**2])
            sage: I.normal_basis()
            [k, 1]
        """
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import \
            kbase_libsingular
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence
        gb = self._groebner_basis_libsingular()
        J = self.ring().ideal(gb)
        if weights is None or degree is None:
            res = kbase_libsingular(J, degree)
        else:
            from sage.libs.singular.function_factory import ff
            res = ff.weightKB(J, -1 if degree is None else degree,
                              tuple(weights), attributes={J: {'isSB': 1}})
        if len(res) == 1 and res[0].is_zero():
            res = []
        return PolynomialSequence(self.ring(), res, immutable=True)

    @require_field
    @handle_AA_and_QQbar
    def normal_basis(self, degree=None, algorithm='libsingular',
                     singular=None):
        """
        Return a vector space basis of the quotient ring of this ideal.

        INPUT:

        - ``degree`` -- integer (default: ``None``)

        - ``algorithm`` -- string (default: ``'libsingular'``); if not the
          default, this will use the ``kbase()`` or ``weightKB()`` command from
          Singular

        - ``singular`` -- the singular interpreter to use when ``algorithm`` is
          not ``'libsingular'`` (default: the default instance)

        OUTPUT:

        Monomials in the basis. If ``degree`` is given, only the monomials of
        the given degree are returned.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = R.ideal(x^2 + y^2 + z^2 - 4, x^2 + 2*y^2 - 5, x*z - 1)
            sage: I.normal_basis()
            [y*z^2, z^2, y*z, z, x*y, y, x, 1]
            sage: I.normal_basis(algorithm='singular')
            [y*z^2, z^2, y*z, z, x*y, y, x, 1]

        The result can be restricted to monomials of a chosen degree, which is
        particularly useful when the quotient ring is not finite-dimensional as
        a vector space.  ::

            sage: J = R.ideal(x^2 + y^2 + z^2 - 4, x^2 + 2*y^2 - 5)
            sage: J.dimension()
            1
            sage: [J.normal_basis(d) for d in (0..3)]
            [[1], [z, y, x], [z^2, y*z, x*z, x*y], [z^3, y*z^2, x*z^2, x*y*z]]
            sage: [J.normal_basis(d, algorithm='singular') for d in (0..3)]
            [[1], [z, y, x], [z^2, y*z, x*z, x*y], [z^3, y*z^2, x*z^2, x*y*z]]

        In case of a polynomial ring with a weighted term order, the degree of
        the monomials is taken with respect to the weights.  ::

            sage: T = TermOrder('wdegrevlex', (1, 2, 3))
            sage: R.<x,y,z> = PolynomialRing(QQ, order=T)
            sage: B = R.ideal(x*y^2 + x^5, z*y + x^3*y).normal_basis(9); B
            [x^2*y^2*z, x^3*z^2, x*y*z^2, z^3]
            sage: all(f.degree() == 9 for f in B)
            True

        TESTS:

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: R.<x,y,z> = QQbar[]
            sage: I = R.ideal(x^2+y^2+z^2-4, x^2+2*y^2-5, x*z-1)
            sage: I.normal_basis()
            [y*z^2, z^2, y*z, z, x*y, y, x, 1]
            sage: J = R.ideal(x^2+y^2+z^2-4, x^2+2*y^2-5)
            sage: [J.normal_basis(d) for d in (0..3)]
            [[1], [z, y, x], [z^2, y*z, x*z, x*y], [z^3, y*z^2, x*z^2, x*y*z]]

        Check the option ``algorithm="singular"`` with a weighted term order::

            sage: T = TermOrder('wdegrevlex', (1, 2, 3))
            sage: S.<x,y,z> = PolynomialRing(GF(2), order=T)                            # needs sage.rings.finite_rings
            sage: S.ideal(x^6 + y^3 + z^2).normal_basis(6, algorithm='singular')        # needs sage.rings.finite_rings
            [x^4*y, x^2*y^2, y^3, x^3*z, x*y*z, z^2]
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence

        weights = tuple(x.degree() for x in self.ring().gens())
        if all(w == 1 for w in weights):
            weights = None

        if algorithm == 'libsingular':
            return self._normal_basis_libsingular(degree, weights=weights)
        else:
            if singular is None:
                singular = singular_default
            gb = self.groebner_basis()
            R = self.ring()
            if degree is None:
                res = singular.kbase(R.ideal(gb))
            elif weights is None:
                res = singular.kbase(R.ideal(gb), int(degree))
            else:
                res = singular.weightKB(R.ideal(gb), int(degree),
                                        singular(weights, type='intvec'))
            return PolynomialSequence(R, [R(f) for f in res], immutable=True)


class MPolynomialIdeal_macaulay2_repr:
    """
    An ideal in a multivariate polynomial ring, which has an underlying
    Macaulay2 ring associated to it.

    EXAMPLES::

        sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
        sage: I = ideal(x*y-z^2, y^2-w^2)
        sage: I
        Ideal (x*y - z^2, y^2 - w^2) of Multivariate Polynomial Ring in x, y, z, w over Integer Ring
    """

    def _groebner_basis_macaulay2(self, strategy=None):
        r"""
        Return the Groebner basis for this ideal, computed using
        Macaulay2.

        ALGORITHM:

        Compute the Groebner basis using the specified strategy in Macaulay2.
        With no strategy option, the Macaulay2 ``gb`` command is used; other
        possible strategies are "f4" and "mgb", which correspond to the "F4"
        and "MGB" strategies in Macaulay2.

        A big advantage of Macaulay2 is that it can compute the Groebner basis
        of ideals in polynomial rings over the integers.

        INPUT:

        - ``strategy`` -- (default: ``'gb'``) argument specifying the strategy
          to be used by Macaulay2; possibilities: ``'f4'``, ``'gb'``, ``'mgb'``

        EXAMPLES::

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = ideal(x*y-z^2, y^2-w^2)
            sage: I.groebner_basis('macaulay2')                # indirect doctest; optional - macaulay2
            [z^4 - x^2*w^2, y*z^2 - x*w^2, x*y - z^2, y^2 - w^2]

        The Groebner basis can be used to compute in
        `\ZZ/n\ZZ[x,\ldots]`.

        ::

            sage: R.<x,y,z> = ZZ[]
            sage: I = ideal([y^2*z - x^3 - 19*x*z, y^2, 19^2])
            sage: I.groebner_basis('macaulay2')               # optional - macaulay2
            [x^3 + 19*x*z, y^2, 361]
            sage: I = ideal([y^2*z - x^3 - 19^2*x*z, y^2, 19^2])
            sage: I.groebner_basis('macaulay2')               # optional - macaulay2
            [x^3, y^2, 361]

        Over finite fields, Macaulay2 supports different algorithms to compute
        Gröbner bases::

            sage: # needs sage.rings.finite_rings
            sage: R = PolynomialRing(GF(101), 'x', 4)
            sage: I = sage.rings.ideal.Cyclic(R)
            sage: gb1 = I.groebner_basis('macaulay2:gb')        # optional - macaulay2
            sage: I = sage.rings.ideal.Cyclic(R)
            sage: gb2 = I.groebner_basis('macaulay2:mgb')       # optional - macaulay2
            sage: I = sage.rings.ideal.Cyclic(R)
            sage: gb3 = I.groebner_basis('macaulay2:f4')        # optional - macaulay2
            sage: gb1 == gb2 == gb3             # optional - macaulay2
            True

        TESTS::

            sage: R.<x,y,z> = ZZ[]
            sage: I = ideal([y^2*z - x^3 - 19*x*z, y^2, 19^2])
            sage: I.groebner_basis('macaulay2:gibberish')     # optional - macaulay2
            Traceback (most recent call last):
            ...
            ValueError: unsupported Macaulay2 strategy
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence

        I = self._macaulay2_()
        if strategy == "gb" or strategy is None:
            m2G = I.gb().generators()
        elif strategy == 'f4':
            m2G = I.groebnerBasis('Strategy=>"F4"')
        elif strategy == 'mgb':
            m2G = I.groebnerBasis('Strategy=>"MGB"')
        else:
            raise ValueError("unsupported Macaulay2 strategy")
        G = str(m2G.external_string()).replace('\n','')
        i = G.rfind('{{')
        j = G.rfind('}}')
        G = G[i+2:j].split(',')
        R = self.ring()
        B = [R(f) for f in G]
        B = PolynomialSequence(self.ring(), B, immutable=True)
        return B

    def _reduce_using_macaulay2(self, f):
        """
        EXAMPLES::

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = ideal(x*y-z^2, y^2-w^2)
            sage: I._reduce_using_macaulay2(x*y-z^2 + y^2)    # optional  - macaulay2
            w^2
        """
        I = self._macaulay2_()
        M2 = I.parent()
        k = M2('(%r) %% %s' % (f, I.name()))
        R = self.ring()
        return R(k)


class NCPolynomialIdeal(MPolynomialIdeal_singular_repr, Ideal_nc):
    def __init__(self, ring, gens, coerce=True, side='left'):
        r"""
        Create a non-commutative polynomial ideal.

        INPUT:

        - ``ring`` -- the g-algebra to which this ideal belongs
        - ``gens`` -- the generators of this ideal
        - ``coerce`` -- boolean (default: ``True``); generators are
          coerced into the ring before creating the ideal
        - ``side`` -- (optional) string; either ``'left'`` (default)
          or ``'twosided'``. Defines whether this ideal is left or two-sided.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: I = H.ideal([y^2, x^2, z^2 - H.one()],             # indirect doctest
            ....:             coerce=False)
            sage: I  # random
            Left Ideal (y^2, x^2, z^2 - 1) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(I.gens(), key=str)
            [x^2, y^2, z^2 - 1]
            sage: H.ideal([y^2, x^2, z^2 - H.one()], side='twosided')  # random
            Twosided Ideal (y^2, x^2, z^2 - 1) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(H.ideal([y^2, x^2, z^2 - H.one()], side='twosided').gens(),
            ....:        key=str)
            [x^2, y^2, z^2 - 1]
            sage: H.ideal([y^2, x^2, z^2 - H.one()], side='right')
            Traceback (most recent call last):
            ...
            ValueError: Only left and two-sided ideals are allowed.
        """
        if side == "right":
            raise ValueError("Only left and two-sided ideals are allowed.")
        Ideal_nc.__init__(self, ring, gens, coerce=coerce, side=side)

    def __call_singular(self, cmd, arg=None):
        """
        Internal function for calling a Singular function.

        INPUT:

        - ``cmd`` - string, representing a Singular function
        - ``arg`` -- (default: ``None``) arguments for which cmd is called

        OUTPUT: result of the Singular function call

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: id = H.ideal(x + y, y + z)
            sage: id.std()  # indirect doctest  # random
            Left Ideal (z, y, x) of Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(id.std().gens(), key=str)
            [x, y, z]
        """
        from sage.libs.singular.function import singular_function
        fun = singular_function(cmd)
        if arg is None:
            return fun(self, ring=self.ring())

        return fun(self, arg, ring=self.ring())

    @cached_method
    def std(self):
        r"""
        Compute a GB of the ideal. It is two-sided if and only if the ideal is two-sided.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: I = H.ideal([y^2, x^2, z^2 - H.one()], coerce=False)
            sage: I.std()  #random
            Left Ideal (z^2 - 1, y*z - y, x*z + x, y^2, 2*x*y - z - 1, x^2) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(I.std().gens(), key=str)
            [2*x*y - z - 1, x*z + x, x^2, y*z - y, y^2, z^2 - 1]


        If the ideal is a left ideal, then :meth:`std` returns a left
        Groebner basis. But if it is a two-sided ideal, then
        the output of :meth:`std` and :meth:`twostd` coincide::

            sage: # needs sage.combinat sage.modules
            sage: JL = H.ideal([x^3, y^3, z^3 - 4*z])
            sage: JL  #random
            Left Ideal (x^3, y^3, z^3 - 4*z) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(JL.gens(), key=str)
            [x^3, y^3, z^3 - 4*z]
            sage: JL.std()  # random
            Left Ideal (z^3 - 4*z, y*z^2 - 2*y*z,
                        x*z^2 + 2*x*z, 2*x*y*z - z^2 - 2*z, y^3, x^3) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(JL.std().gens(), key=str)
            [2*x*y*z - z^2 - 2*z, x*z^2 + 2*x*z, x^3, y*z^2 - 2*y*z, y^3, z^3 - 4*z]
            sage: JT = H.ideal([x^3, y^3, z^3 - 4*z], side='twosided')
            sage: JT  #random
            Twosided Ideal (x^3, y^3, z^3 - 4*z) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(JT.gens(), key=str)
            [x^3, y^3, z^3 - 4*z]
            sage: JT.std()  #random
            Twosided Ideal (z^3 - 4*z, y*z^2 - 2*y*z, x*z^2 + 2*x*z,
                            y^2*z - 2*y^2, 2*x*y*z - z^2 - 2*z, x^2*z + 2*x^2,
                            y^3, x*y^2 - y*z, x^2*y - x*z - 2*x, x^3) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: sorted(JT.std().gens(), key=str)
            [2*x*y*z - z^2 - 2*z, x*y^2 - y*z, x*z^2 + 2*x*z, x^2*y - x*z - 2*x,
             x^2*z + 2*x^2, x^3, y*z^2 - 2*y*z, y^2*z - 2*y^2, y^3, z^3 - 4*z]
            sage: JT.std() == JL.twostd()
            True

        ALGORITHM: Uses Singular's ``std`` command
        """
        if self.side() == 'twosided':
            return self.twostd()
        return self.ring().ideal( self.__call_singular('std'), side=self.side())
#        return self.__call_singular('std')

    def elimination_ideal(self, variables):
        r"""
        Return the elimination ideal of this ideal with respect to the
        variables given in ``variables``.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: I = H.ideal([y^2, x^2, z^2 - H.one()], coerce=False)
            sage: I.elimination_ideal([x, z])
            Left Ideal (y^2) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {...}
            sage: J = I.twostd(); J
            Twosided Ideal (z^2 - 1, y*z - y, x*z + x, y^2, 2*x*y - z - 1, x^2) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {...}
            sage: J.elimination_ideal([x, z])
            Twosided Ideal (y^2) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {...}


        ALGORITHM: Uses Singular's ``eliminate`` command
        """
        from sage.misc.misc_c import prod
        if self.side() == 'twosided':
            J = self.twostd()
        else:
            J = self
        return J.ring().ideal(J.__call_singular('eliminate', prod(variables)), side=self.side())
#        return self.__call_singular('std')

    @cached_method
    def twostd(self):
        r"""
        Compute a two-sided GB of the ideal (even if it is a left ideal).

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: I = H.ideal([y^2, x^2, z^2 - H.one()], coerce=False)
            sage: I.twostd()  #random
            Twosided Ideal (z^2 - 1, y*z - y, x*z + x, y^2, 2*x*y - z - 1, x^2) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field...
            sage: sorted(I.twostd().gens(), key=str)
            [2*x*y - z - 1, x*z + x, x^2, y*z - y, y^2, z^2 - 1]

        ALGORITHM: Uses Singular's ``twostd`` command
        """
        return self.ring().ideal( self.__call_singular('twostd'), side='twosided')
#        return self.__call_singular('twostd')

#    def syz(self):
#        return self.__call_singular('syz')

    @cached_method
    def _groebner_strategy(self):
        """
        Return Singular's Groebner Strategy object for the Groebner
        basis of this ideal which implements some optimized functions.

        EXAMPLES::

           sage: # needs sage.combinat sage.modules
           sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
           sage: H.<x,y,z> = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
           sage: I = H.ideal([y^2, x^2, z^2-H.one()], coerce=False)
           sage: I._groebner_strategy()  # random
           Groebner Strategy for ideal generated by 6 elements over
            Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
            nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}

        .. NOTE::

            This function is mainly used internally.
        """
        from sage.libs.singular.groebner_strategy import NCGroebnerStrategy
        return NCGroebnerStrategy(self.std())

    def reduce(self,p):
        """
        Reduce an element modulo a Groebner basis for this ideal.

        It returns 0 if and only if the element is in this ideal. In any
        case, this reduction is unique up to monomial orders.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H.<x,y,z> = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: I = H.ideal([y^2, x^2, z^2 - H.one()],
            ....:             coerce=False, side='twosided')
            sage: Q = H.quotient(I); Q  #random
            Quotient of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
             by the ideal (y^2, x^2, z^2 - 1)
            sage: Q.2^2 == Q.one()   # indirect doctest
            True

        Here, we see that the relation that we just found in the quotient
        is actually a consequence of the given relations::

            sage: H.2^2 - H.one() in I.std().gens()                                     # needs sage.combinat sage.modules
            True

        Here is the corresponding direct test::

            sage: I.reduce(z^2)                                                         # needs sage.combinat sage.modules
            1
        """
        return self._groebner_strategy().normal_form(p)

    def _contains_(self,p):
        """
        EXAMPLES:

        We define a left and a two-sided ideal::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H.<x,y,z> = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: JL = H.ideal([x^3, y^3, z^3 - 4*z])
            sage: JL.std()  #random
            Left Ideal (z^3 - 4*z, y*z^2 - 2*y*z, x*z^2 + 2*x*z, 2*x*y*z - z^2 - 2*z, y^3, x^3) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}
            sage: JT = H.ideal([x^3, y^3, z^3 - 4*z], side='twosided')
            sage: JT.std()  #random
            Twosided Ideal (z^3 - 4*z, y*z^2 - 2*y*z, x*z^2 + 2*x*z, y^2*z - 2*y^2, 2*x*y*z - z^2 - 2*z,
                            x^2*z + 2*x^2, y^3, x*y^2 - y*z, x^2*y - x*z - 2*x, x^3) of
             Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
             nc-relations: {z*x: x*z + 2*x, z*y: y*z - 2*y, y*x: x*y - z}

        Apparently, ``x*y^2-y*z`` should be in the two-sided, but not
        in the left ideal::

            sage: x*y^2-y*z in JL   #indirect doctest                                   # needs sage.combinat sage.modules
            False
            sage: x*y^2-y*z in JT                                                       # needs sage.combinat sage.modules
            True
        """
        return self.reduce(p).is_zero()

    @require_field
    def syzygy_module(self):
        r"""
        Compute the first syzygy (i.e., the module of relations of the
        given generators) of the ideal.

        NOTE:

        Only left syzygies can be computed. So, even if the ideal is
        two-sided, then the syzygies are only one-sided. In that case,
        a warning is printed.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: I = H.ideal([y^2, x^2, z^2-H.one()], coerce=False)
            sage: G = vector(I.gens()); G
            d...: UserWarning: You are constructing a free module
            over a noncommutative ring. Sage does not have a concept
            of left/right and both sided modules, so be careful.
            It's also not guaranteed that all multiplications are
            done from the right side.
            d...: UserWarning: You are constructing a free module
            over a noncommutative ring. Sage does not have a concept
            of left/right and both sided modules, so be careful.
            It's also not guaranteed that all multiplications are
            done from the right side.
            (y^2, x^2, z^2 - 1)
            sage: M = I.syzygy_module(); M
            [                                                                         -z^2 - 8*z - 15                                                                                        0                                                                                      y^2]
            [                                                                                       0                                                                          -z^2 + 8*z - 15                                                                                      x^2]
            [                                                              x^2*z^2 + 8*x^2*z + 15*x^2                                                              -y^2*z^2 + 8*y^2*z - 15*y^2                                                                   -4*x*y*z + 2*z^2 + 2*z]
            [                 x^2*y*z^2 + 9*x^2*y*z - 6*x*z^3 + 20*x^2*y - 72*x*z^2 - 282*x*z - 360*x                                                              -y^3*z^2 + 7*y^3*z - 12*y^3                                                                                  6*y*z^2]
            [                                                              x^3*z^2 + 7*x^3*z + 12*x^3                 -x*y^2*z^2 + 9*x*y^2*z - 4*y*z^3 - 20*x*y^2 + 52*y*z^2 - 224*y*z + 320*y                                                                                 -6*x*z^2]
            [  x^2*y^2*z + 4*x^2*y^2 - 8*x*y*z^2 - 48*x*y*z + 12*z^3 - 64*x*y + 108*z^2 + 312*z + 288                                                                           -y^4*z + 4*y^4                                                                                        0]
            [                                                  2*x^3*y*z + 8*x^3*y + 9*x^2*z + 27*x^2                                   -2*x*y^3*z + 8*x*y^3 - 12*y^2*z^2 + 99*y^2*z - 195*y^2                                                                -36*x*y*z + 24*z^2 + 18*z]
            [                                                                           x^4*z + 4*x^4    -x^2*y^2*z + 4*x^2*y^2 - 4*x*y*z^2 + 32*x*y*z - 6*z^3 - 64*x*y + 66*z^2 - 240*z + 288                                                                                        0]
            [x^3*y^2*z + 4*x^3*y^2 + 18*x^2*y*z - 36*x*z^3 + 66*x^2*y - 432*x*z^2 - 1656*x*z - 2052*x                                      -x*y^4*z + 4*x*y^4 - 8*y^3*z^2 + 62*y^3*z - 114*y^3                                                                        48*y*z^2 - 36*y*z]

            sage: M*G                                                                   # needs sage.combinat sage.modules
            (0, 0, 0, 0, 0, 0, 0, 0, 0)

        ALGORITHM: Uses Singular's ``syz`` command
        """
        if self.side() == 'twosided':
            warn("The result of this Syzygy computation is one-sided (left)!")
        from sage.libs.singular.function_factory import ff
        syz = ff.syz
        from sage.matrix.constructor import matrix

        #return self._singular_().syz().transpose().sage_matrix(self.ring())
        S = syz(self)
        return matrix(self.ring(), S)

    def res(self, length):
        r"""
        Compute the resolution up to a given length of the ideal.

        NOTE:

        Only left syzygies can be computed. So, even if the ideal is
        two-sided, then the resolution is only one-sided. In that case,
        a warning is printed.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: H = A.g_algebra({y*x: x*y-z, z*x: x*z+2*x, z*y: y*z-2*y})
            sage: H.inject_variables()
            Defining x, y, z
            sage: I = H.ideal([y^2, x^2, z^2-H.one()], coerce=False)
            sage: I.res(3)
            <Resolution>
        """
        if self.side() == 'twosided':
            warn("The resulting resolution is one-sided (left)!")
        return self.__call_singular('res', length)


@richcmp_method
class MPolynomialIdeal(MPolynomialIdeal_singular_repr,
                       MPolynomialIdeal_macaulay2_repr,
                       MPolynomialIdeal_magma_repr,
                       Ideal_generic):
    def __init__(self, ring, gens, coerce=True):
        r"""
        Create an ideal in a multivariate polynomial ring.

        INPUT:

        - ``ring`` -- the ring the ideal is defined in

        - ``gens`` -- list of generators for the ideal

        - ``coerce`` -- whether to coerce elements to the ring ``ring``

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(IntegerRing(), 2, order='lex')
            sage: R.ideal([x, y])
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: R.<x0,x1> = GF(3)[]
            sage: R.ideal([x0^2, x1^3])
            Ideal (x0^2, x1^3) of Multivariate Polynomial Ring in x0, x1 over Finite Field of size 3
        """
        Ideal_generic.__init__(self, ring, gens, coerce=coerce)
        self._gb_by_ordering = dict()

    def __hash__(self):
        r"""
        Stupid constant hash function.

        TESTS::

            sage: R.<x,y> = PolynomialRing(IntegerRing(), 2, order='lex')
            sage: hash(R.ideal([x, y]))
            0
        """
        return 0

    @cached_method
    def gens(self):
        """
        Return a set of generators / a basis of this ideal. This is usually the
        set of generators provided during object creation.

        EXAMPLES::

           sage: P.<x,y> = PolynomialRing(QQ,2)
           sage: I = Ideal([x, y + 1]); I
           Ideal (x, y + 1) of Multivariate Polynomial Ring in x, y over Rational Field
           sage: I.gens()
           [x, y + 1]
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence
        return PolynomialSequence(self.ring(), Ideal_generic.gens(self), immutable=True)

    @property
    def basis(self):
        """
        Shortcut to ``gens()``.

        EXAMPLES::

           sage: P.<x,y> = PolynomialRing(QQ,2)
           sage: I = Ideal([x, y + 1])
           sage: I.basis
           [x, y + 1]
        """
        return self.gens()

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polynomial ideal

        OUTPUT: boolean

        ALGORITHM:

        Comparison for ``==`` and ``!=`` compares two Groebner bases.

        Comparison for ``<=`` and ``>=`` tests the inclusion of ideals
        using the usual ideal membership test, namely all generators
        of one ideal must reduce to zero in the other ideal's Groebner
        basis.

        Comparison for ``<`` and ``>`` tests for inclusion and different
        Groebner bases.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]; I = R*[x^2 + y, 2*y]; J = R*[x^2 + y]
            sage: I > J
            True
            sage: J < I
            True
            sage: I == I
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = GF(32003)[]
            sage: I = R*[x^2 + x, y]
            sage: J = R*[x + 1, y]
            sage: J < I
            False
            sage: I < J
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = GF(32003)[]
            sage: I = R*[x^2 + x, y]
            sage: J = R*[x + 1, y]
            sage: J > I
            True
            sage: I > J
            False

        ::

            sage: R = PolynomialRing(QQ,'x,y,z')
            sage: I = R.ideal()
            sage: I == R.ideal()
            True

        ::

            sage: R = PolynomialRing(QQ, names=[])
            sage: R.ideal(0) == R.ideal(0)
            True

        ::

            sage: R.<x,y> = QQ[]
            sage: I = (x^3 + y, y) * R
            sage: J = (x^3 + y, y, y*x^3 + y^2) * R
            sage: I == J
            True

        ::

            sage: R = PolynomialRing(QQ, 'x,y,z', order='degrevlex')
            sage: S = PolynomialRing(QQ, 'x,y,z', order='invlex')
            sage: I = R.ideal([R.0,R.1])
            sage: J = S.ideal([S.0,S.1])
            sage: I == J
            True

        TESTS:

        We test to make sure that pickling works with the cached
        Groebner basis::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = GF(32003)[]
            sage: I = R*[x^2 + x, y]
            sage: J = R*[x + 1, y]
            sage: J >= I
            True
            sage: I >= J
            False

            sage: loads(dumps(I)).__getstate__()                                        # needs sage.rings.finite_rings
            (Monoid of ideals of Multivariate Polynomial Ring in x, y over Finite Field of size 32003,
             {'_Ideal_generic__gens': (x^2 + x, y),
              '_Ideal_generic__ring': Multivariate Polynomial Ring in x, y over Finite Field of size 32003,
              '_cache__groebner_basis': {},
              '_gb_by_ordering': {'degrevlex': [x^2 + x, y]},
              'gens': Pickle of the cached method "gens",
              'groebner_basis': Pickle of the cached method "groebner_basis"})

        This example checks :issue:`12802`::

            sage: R.<x,y> = ZZ[]
            sage: I = R * [ x^2 + y, 2*y ]
            sage: J = R * [ x^2 - y, 2*y ]
            sage: I == J
            True

        Another good test from the discussion in :issue:`12802`::

            sage: Rx = PolynomialRing(QQ, 2, "x")
            sage: Ix = Rx.ideal(Rx.0)
            sage: Ry = PolynomialRing(QQ, 2, "y")
            sage: Iy = Ry.ideal(Ry.0)
            sage: Ix == Iy
            False

        However, this should work if only the orderings are different::

            sage: R = PolynomialRing(QQ, 'x', 2, order='degrevlex')
            sage: S = PolynomialRing(QQ, 'x', 2, order='lex')
            sage: R == S
            False
            sage: I = R*[R.0^2 + R.1, R.1]
            sage: J = S*[S.0^2 + S.1, S.1]
            sage: I == J
            True
        """
        # first check the type
        if not isinstance(other, MPolynomialIdeal):
            return NotImplemented

        if self is other:
            return rich_to_bool(op, 0)

        # comparison for >= and > : swap the arguments
        if op == op_GE:
            return other.__richcmp__(self, op_LE)
        elif op == op_GT:
            return other.__richcmp__(self, op_LT)

        # the ideals may be defined w.r.t. to different term orders
        # but are still the same.
        R = self.ring()
        S = other.ring()

        # separate next two tests to avoid unnecessary creation of
        # Groebner basis
        if S is not R:
            if S.change_ring(order=R.term_order()) != R:  # rings are unique
                return NotImplemented
            else:
                # at this point, the rings are the same, but for the term order,
                # and we can fix that easily
                other_new = other.change_ring(R)
        else:
            other_new = other

        s_gens = self.gens()
        o_gens = other_new.gens()
        try:
            if (s_gens == o_gens) or (set(s_gens) == set(o_gens)):
                # the first clause works in the non-hashable case
                return rich_to_bool(op, 0)
        except TypeError:
            pass

        # comparison for <=
        # needs just the Groebner basis for other
        if op == op_LE:
            l = self.gens()
            try:
                # first check whether the GB is cached already
                if other_new.groebner_basis.is_in_cache():
                    r = other_new.groebner_basis()
                elif len(other_new._gb_by_ordering):
                    o, r = next(iter(other_new._gb_by_ordering.items()))
                    l = self.change_ring(R.change_ring(order=o)).gens()
                else: # use easy GB otherwise
                    newR = R.change_ring(order='degrevlex')
                    l = self.change_ring(newR).gens()
                    r = other_new.change_ring(newR).groebner_basis()
                    # remember this Groebner basis for future reference
                    other_new._gb_by_ordering['degrevlex'] = r
            except AttributeError: # e.g. quotient rings
                r = other_new.groebner_basis()
            return all(f.reduce(r) == 0 for f in l)

        # comparison for == and != and <
        # needs both groebner bases for the same term order
        # first check whether the GBs are cached already
        if op in [op_EQ, op_NE, op_LT]:
            try:
                if (other_new.groebner_basis.is_in_cache()
                    or self.groebner_basis().is_in_cache()):
                    l = self.groebner_basis()
                    r = other_new.groebner_basis()
                else: # use easy GB otherwise
                    newR = R.change_ring(order='degrevlex')
                    l = self.change_ring(newR).groebner_basis()
                    r = other_new.change_ring(newR).groebner_basis()
            except AttributeError: # e.g. quotient rings
                l = self.groebner_basis()
                r = other_new.groebner_basis()
            contained = all(f.reduce(r) == 0 for f in l)
            contains = all(f.reduce(l) == 0 for f in r)
            if op == op_EQ:
                return contained and contains
            elif op == op_NE:
                return not (contained and contains)
            else:  # remaining case <
                return contained and not contains

    def groebner_fan(self, is_groebner_basis=False, symmetry=None, verbose=False):
        r"""
        Return the Groebner fan of this ideal.

        The base ring must be `\QQ` or a finite field
        `\GF{p}` of with `p \leq 32749`.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: i = ideal(x^2 - y^2 + 1)
            sage: g = i.groebner_fan()
            sage: g.reduced_groebner_bases()
            [[x^2 - y^2 + 1], [-x^2 + y^2 - 1]]

        INPUT:

        - ``is_groebner_basis`` -- boolean (default: ``False``); if
          ``True``, then ``I.gens()`` must be a Groebner basis with respect to
          the standard degree lexicographic term order

        - ``symmetry`` -- (default: ``None``) if not ``None``, describes
          symmetries of the ideal

        - ``verbose`` -- (default: ``False``) if ``True``, printout
          useful info during computations
        """
        import sage.rings.polynomial.groebner_fan as groebner_fan
        return groebner_fan.GroebnerFan(self, is_groebner_basis=is_groebner_basis,
                                        symmetry=symmetry, verbose=verbose)

    @cached_method(do_pickle=True)
    @handle_AA_and_QQbar
    def groebner_basis(self, algorithm='', deg_bound=None, mult_bound=None, prot=False, *args, **kwds):
        r"""
        Return the reduced Groebner basis of this ideal.

        A Groebner basis `g_1,...,g_n` for an ideal `I` is a
        generating set such that `<LM(g_i)> = LM(I)`, i.e., the
        leading monomial ideal of `I` is spanned by the leading terms
        of `g_1,...,g_n`. Groebner bases are the key concept in
        computational ideal theory in multivariate polynomial rings
        which allows a variety of problems to be solved.

        Additionally, a *reduced* Groebner basis `G` is a unique
        representation for the ideal `<G>` with respect to the chosen
        monomial ordering.

        INPUT:

        - ``algorithm`` -- determines the algorithm to use, see below
          for available algorithms

        - ``deg_bound`` -- only compute to degree ``deg_bound``, that
          is, ignore all S-polynomials of higher degree. (default:
          ``None``)

        - ``mult_bound`` -- the computation is stopped if the ideal is
          zero-dimensional in a ring with local ordering and its
          multiplicity is lower than ``mult_bound``. Singular
          only. (default: ``None``)

        - ``prot`` -- if set to ``True`` the computation protocol of
          the underlying implementation is printed. If an algorithm
          from the ``singular:`` or ``magma:`` family is used,
          ``prot`` may also be ``sage`` in which case the output is
          parsed and printed in a common format where the amount of
          information printed can be controlled via calls to
          :func:`set_verbose`.

        - ``*args`` -- additional parameters passed to the respective
          implementations

        - ``**kwds`` -- additional keyword parameters passed to the
          respective implementations

        ALGORITHMS:

        ``''``
            autoselect (default)

        ``'singular:groebner'``
            Singular's ``groebner`` command

        ``'singular:std'``
            Singular's ``std`` command

        ``'singular:stdhilb'``
            Singular's ``stdhib`` command

        ``'singular:stdfglm'``
            Singular's ``stdfglm`` command

        ``'singular:slimgb'``
            Singular's ``slimgb`` command

        ``'libsingular:groebner'``
            libSingular's ``groebner`` command

        ``'libsingular:std'``
            libSingular's ``std`` command

        ``'libsingular:slimgb'``
            libSingular's ``slimgb`` command

        ``'libsingular:stdhilb'``
            libSingular's ``stdhib`` command

        ``'libsingular:stdfglm'``
            libSingular's ``stdfglm`` command

        ``'toy:buchberger'``
            Sage's toy/educational buchberger without Buchberger criteria

        ``'toy:buchberger2'``
            Sage's toy/educational buchberger with Buchberger criteria

        ``'toy:d_basis'``
            Sage's toy/educational algorithm for computation over PIDs

        ``'macaulay2:gb'``
            Macaulay2's ``gb`` command (if available)

        ``'macaulay2:f4'``
            Macaulay2's ``GroebnerBasis`` command with the strategy "F4" (if available)

        ``'macaulay2:mgb'``
            Macaulay2's ``GroebnerBasis`` command with the strategy "MGB" (if available)

        ``'msolve'``
            :ref:`optional package msolve <spkg_msolve>` (degrevlex order)

        ``'magma:GroebnerBasis'``
            Magma's ``Groebnerbasis`` command (if available)

        ``'ginv:TQ'``, ``'ginv:TQBlockHigh'``, ``'ginv:TQBlockLow'`` and ``'ginv:TQDegree'``
            One of GINV's implementations (if available)

        ``'giac:gbasis'``
            Giac's ``gbasis`` command (if available)

        If only a system is given - e.g. ``'magma'`` - the default algorithm is
        chosen for that system.

        .. NOTE::

            The Singular and libSingular versions of the respective
            algorithms are identical, but the former calls an external
            Singular process while the latter calls a C function,
            and thus the calling overhead is smaller. However, the
            libSingular interface does not support pretty printing of
            computation protocols.

        EXAMPLES:

        Consider Katsura-3 over `\QQ` with lexicographical term
        ordering. We compute the reduced Groebner basis using every
        available implementation and check their equality.

        ::

            sage: P.<a,b,c> = PolynomialRing(QQ,3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis()
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:groebner')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:std')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:stdhilb')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:stdfglm')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:slimgb')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        Although Giac does support lexicographical ordering, we use degree
        reverse lexicographical ordering here, in order to test against
        :issue:`21884`::

            sage: # needs sage.libs.giac
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: J = I.change_ring(P.change_ring(order='degrevlex'))
            sage: gb = J.groebner_basis('giac')  # random
            sage: gb
            [c^3 - 79/210*c^2 + 1/30*b + 1/70*c, b^2 - 3/5*c^2 - 1/5*b + 1/5*c, b*c + 6/5*c^2 - 1/10*b - 2/5*c, a + 2*b + 2*c - 1]
            sage: J.groebner_basis.set_cache(gb)
            sage: ideal(J.transformed_basis()).change_ring(P).interreduced_basis()  # testing issue #21884
            ...[a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        Giac's gbasis over `\QQ` can benefit from a probabilistic lifting and
        multi threaded operations::

            sage: # needs sage.libs.giac
            sage: A9 = PolynomialRing(QQ, 9, 'x')
            sage: I9 = sage.rings.ideal.Katsura(A9)
            sage: print("possible output from giac", flush=True); I9.groebner_basis("giac", proba_epsilon=1e-7)  # long time (3s)
            possible output...
            Polynomial Sequence with 143 Polynomials in 9 Variables

        The list of available Giac options is provided at :func:`sage.libs.giac.groebner_basis`.

        Note that ``toy:buchberger`` does not return the reduced Groebner
        basis, ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: gb = I.groebner_basis('toy:buchberger')
            sage: gb.is_groebner()
            True
            sage: gb == gb.reduced()
            False

        but that ``toy:buchberger2`` does. ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: gb = I.groebner_basis('toy:buchberger2'); gb
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]
            sage: gb == gb.reduced()
            True

        Here we use Macaulay2 with three different strategies over a finite
        field. ::

            sage: # optional - macaulay2
            sage: R.<a,b,c> = PolynomialRing(GF(101), 3)
            sage: I = sage.rings.ideal.Katsura(R,3)  # regenerate to prevent caching
            sage: I.groebner_basis('macaulay2:gb')
            [c^3 + 28*c^2 - 37*b + 13*c, b^2 - 41*c^2 + 20*b - 20*c,
             b*c - 19*c^2 + 10*b + 40*c, a + 2*b + 2*c - 1]
            sage: I = sage.rings.ideal.Katsura(R,3)  # regenerate to prevent caching
            sage: I.groebner_basis('macaulay2:f4')
            [c^3 + 28*c^2 - 37*b + 13*c, b^2 - 41*c^2 + 20*b - 20*c,
             b*c - 19*c^2 + 10*b + 40*c, a + 2*b + 2*c - 1]
            sage: I = sage.rings.ideal.Katsura(R,3)  # regenerate to prevent caching
            sage: I.groebner_basis('macaulay2:mgb')
            [c^3 + 28*c^2 - 37*b + 13*c, b^2 - 41*c^2 + 20*b - 20*c,
             b*c - 19*c^2 + 10*b + 40*c, a + 2*b + 2*c - 1]

        Over prime fields of small characteristic, we can also use the
        :ref:`optional package msolve <spkg_msolve>`::

            sage: R.<a,b,c> = PolynomialRing(GF(101), 3)
            sage: I = sage.rings.ideal.Katsura(R,3)  # regenerate to prevent caching
            sage: I.groebner_basis('msolve')                    # optional - msolve
            [a + 2*b + 2*c - 1, b*c - 19*c^2 + 10*b + 40*c,
             b^2 - 41*c^2 + 20*b - 20*c, c^3 + 28*c^2 - 37*b + 13*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('magma:GroebnerBasis')       # optional - magma
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1,
             b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        Singular and libSingular can compute Groebner basis with degree
        restrictions. ::

            sage: R.<x,y> = QQ[]
            sage: I = R*[x^3 + y^2, x^2*y + 1]
            sage: I.groebner_basis(algorithm='singular')
            [x^3 + y^2, x^2*y + 1, y^3 - x]
            sage: I.groebner_basis(algorithm='singular', deg_bound=2)
            [x^3 + y^2, x^2*y + 1]
            sage: I.groebner_basis()
            [x^3 + y^2, x^2*y + 1, y^3 - x]
            sage: I.groebner_basis(deg_bound=2)
            [x^3 + y^2, x^2*y + 1]

        A protocol is printed, if the verbosity level is at least 2,
        or if the argument ``prot`` is provided. Historically, the
        protocol did not appear during doctests, so, we skip the
        examples with protocol output.  ::

            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(2)
            sage: I = R*[x^3+y^2,x^2*y+1]
            sage: I.groebner_basis()  # not tested
            std in (QQ),(x,y),(dp(2),C)
            [...:2]3ss4s6
            (S:2)--
            product criterion:1 chain criterion:0
            [x^3 + y^2, x^2*y + 1, y^3 - x]
            sage: I.groebner_basis(prot=False)
            std in (QQ),(x,y),(dp(2),C)
            [...:2]3ss4s6
            (S:2)--
            product criterion:1 chain criterion:0
            [x^3 + y^2, x^2*y + 1, y^3 - x]
            sage: set_verbose(0)
            sage: I.groebner_basis(prot=True)  # not tested
            std in (QQ),(x,y),(dp(2),C)
            [...:2]3ss4s6
            (S:2)--
            product criterion:1 chain criterion:0
            [x^3 + y^2, x^2*y + 1, y^3 - x]

        The list of available options is provided at
        :class:`~sage.libs.singular.option.LibSingularOptions`.

        Note that Groebner bases over `\ZZ` can also be computed. ::

            sage: P.<a,b,c> = PolynomialRing(ZZ,3)
            sage: I = P * (a + 2*b + 2*c - 1, a^2 - a + 2*b^2 + 2*c^2, 2*a*b + 2*b*c - b)
            sage: I.groebner_basis()
            [b^3 + b*c^2 + 12*c^3 + b^2 + b*c - 4*c^2,
             2*b*c^2 - 6*c^3 - b^2 - b*c + 2*c^2,
             42*c^3 + b^2 + 2*b*c - 14*c^2 + b,
             2*b^2 + 6*b*c + 6*c^2 - b - 2*c,
             10*b*c + 12*c^2 - b - 4*c,
             a + 2*b + 2*c - 1]

        ::

            sage: I.groebner_basis('macaulay2')                 # optional - macaulay2
            [b^3 + b*c^2 + 12*c^3 + b^2 + b*c - 4*c^2,
             2*b*c^2 - 6*c^3 + b^2 + 5*b*c + 8*c^2 - b - 2*c,
             42*c^3 + b^2 + 2*b*c - 14*c^2 + b,
             2*b^2 - 4*b*c - 6*c^2 + 2*c, 10*b*c + 12*c^2 - b - 4*c,
             a + 2*b + 2*c - 1]

        Groebner bases over `\ZZ/n\ZZ` are also supported::

            sage: P.<a,b,c> = PolynomialRing(Zmod(1000), 3)
            sage: I = P * (a + 2*b + 2*c - 1, a^2 - a + 2*b^2 + 2*c^2, 2*a*b + 2*b*c - b)
            sage: I.groebner_basis()
            [b*c^2 + 732*b*c + 808*b,
             2*c^3 + 884*b*c + 666*c^2 + 320*b,
             b^2 + 438*b*c + 281*b,
             5*b*c + 156*c^2 + 112*b + 948*c,
             50*c^2 + 600*b + 650*c,
             a + 2*b + 2*c + 999,
             125*b]

        ::

            sage: R.<x,y,z> = PolynomialRing(Zmod(2233497349584))
            sage: I = R.ideal([z*(x-3*y), 3^2*x^2-y*z, z^2+y^2])
            sage: I.groebner_basis()
            [2*z^4, y*z^2 + 81*z^3, 248166372176*z^3, 9*x^2 - y*z, y^2 + z^2, x*z +
            2233497349581*y*z, 248166372176*y*z]

        Sage also supports local orderings::

            sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='negdegrevlex')
            sage: I = P * (  x*y*z + z^5, 2*x^2 + y^3 + z^7, 3*z^5 + y^5 )
            sage: I.groebner_basis()
            [x^2 + 1/2*y^3, x*y*z + z^5, y^5 + 3*z^5, y^4*z - 2*x*z^5, z^6]

        We can represent every element in the ideal as a combination
        of the generators using the :meth:`~sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict.lift` method::

            sage: P.<x,y,z> = PolynomialRing(QQ, 3)
            sage: I = P * ( x*y*z + z^5, 2*x^2 + y^3 + z^7, 3*z^5 + y^5 )
            sage: J = Ideal(I.groebner_basis())
            sage: f = sum(P.random_element(terms=2)*f for f in I.gens())
            sage: f                        # random
            1/2*y^2*z^7 - 1/4*y*z^8 + 2*x*z^5 + 95*z^6 + 1/2*y^5 - 1/4*y^4*z + x^2*y^2 + 3/2*x^2*y*z + 95*x*y*z^2
            sage: f.lift(I.gens())         # random
            [2*x + 95*z, 1/2*y^2 - 1/4*y*z, 0]
            sage: l = f.lift(J.gens()); l  # random
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/2*y^2 + 1/4*y*z, 1/2*y^2*z^2 - 1/4*y*z^3 + 2*x + 95*z]
            sage: sum(map(mul, zip(l,J.gens()))) == f
            True

        Groebner bases over fraction fields of polynomial rings are also supported::

            sage: P.<t> = QQ[]
            sage: F = Frac(P)
            sage: R.<X,Y,Z> = F[]
            sage: I = Ideal([f + P.random_element() for f in sage.rings.ideal.Katsura(R).gens()])
            sage: I.groebner_basis().ideal() == I
            True

        In cases where a characteristic cannot be determined, we use a toy implementation of Buchberger's algorithm
        (see :issue:`6581`)::

            sage: R.<a,b> = QQ[]; I = R.ideal(a^2+b^2-1)
            sage: Q = QuotientRing(R,I); K = Frac(Q)
            sage: R2.<x,y> = K[]; J = R2.ideal([(a^2+b^2)*x + y, x+y])
            sage: J.groebner_basis()
            verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
            [x + y]

        ALGORITHM:

        Uses Singular, one of the other systems listed above (if available),
        or a toy implementation.

        TESTS:

        Check :issue:`27445`::

            sage: P = PolynomialRing(QQ, 't', 0)
            sage: P.ideal([P(2)]).groebner_basis()
            [1]
            sage: P.ideal([]).groebner_basis()
            [0]
            sage: P.ideal([0]).groebner_basis()
            [0]
            sage: P.ideal([3, 4, 0, 5]).groebner_basis()
            [1]

            sage: P = PolynomialRing(ZZ, 't', 0)
            sage: P.ideal([P(2)]).groebner_basis()
            [2]
            sage: P.ideal([]).groebner_basis()
            [0]
            sage: P.ideal([0]).groebner_basis()
            [0]
            sage: P.ideal([2, 4, 6]).groebner_basis()
            [2]

            sage: P = PolynomialRing(Zmod(8), 't', 0)
            sage: P.ideal([P(2)]).groebner_basis()
            [2]
            sage: P.ideal([]).groebner_basis()
            [0]
            sage: P.ideal([0]).groebner_basis()
            [0]
            sage: P.ideal([P(3)]).groebner_basis()
            [1]

        Check that this method works over QQbar (:issue:`25351`)::

            sage: # needs sage.rings.number_field
            sage: P.<a,b,c> = PolynomialRing(QQbar, 3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis()
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:groebner')
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:std')
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:stdhilb')
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:stdfglm')
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('libsingular:slimgb')
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]

            sage: # needs sage.rings.number_field sage.libs.giac
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: J = I.change_ring(P.change_ring(order='degrevlex'))
            sage: gb = J.groebner_basis('giac')  # random
            sage: gb
            [c^3 + (-79/210)*c^2 + 1/30*b + 1/70*c, b^2 + (-3/5)*c^2 + (-1/5)*b + 1/5*c, b*c + 6/5*c^2 + (-1/10)*b + (-2/5)*c, a + 2*b + 2*c - 1]

            sage: # needs sage.rings.number_field
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('toy:buchberger2')
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('macaulay2:gb')      # optional - macaulay2
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
            sage: I = sage.rings.ideal.Katsura(P,3)  # regenerate to prevent caching
            sage: I.groebner_basis('magma:GroebnerBasis')       # optional - magma
            [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 + (-79/7)*c^2 + 3/7*c, c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]

        msolve currently supports the degrevlex order only::

            sage: R.<a,b,c> = PolynomialRing(GF(101), 3, order='lex')
            sage: I = sage.rings.ideal.Katsura(R,3)  # regenerate to prevent caching
            sage: I.groebner_basis('msolve')                    # optional - msolve
            Traceback (most recent call last):
            ...
            NotImplementedError: msolve only supports the degrevlex order (use transformed_basis())
        """
        from sage.rings.polynomial.multi_polynomial_sequence import \
            PolynomialSequence
        from sage.rings.polynomial.polynomial_ring_constructor import \
            PolynomialRing

        if algorithm.lower() == "magma":
            algorithm = "magma:GroebnerBasis"
        elif algorithm.lower() == "singular":
            algorithm = "singular:groebner"
        elif algorithm.lower() == "libsingular":
            algorithm = "libsingular:groebner"
        elif algorithm.lower() == "macaulay2":
            algorithm = "macaulay2:gb"
        elif algorithm.lower() == "toy":
            algorithm = "toy:buchberger2"
        elif algorithm.lower() == "giac":
            algorithm = "giac:gbasis"

        if not algorithm:
            try:
                gb = self._groebner_basis_libsingular("groebner", deg_bound=deg_bound, mult_bound=mult_bound, *args, **kwds)
            except (TypeError, NameError, ImportError):  # conversion to Singular not supported
                try:
                    gb = self._groebner_basis_singular("groebner", deg_bound=deg_bound, mult_bound=mult_bound, *args, **kwds)
                except (TypeError, NameError, NotImplementedError, ImportError):  # conversion to Singular not supported
                    R = self.ring()
                    B = R.base_ring()
                    if R.ngens() == 0:
                        if R.base_ring().is_field():
                            if any(g != 0 for g in self.gens()):
                                gb = [R.one()]
                            else:
                                gb = [R.zero()]
                        else:
                            # To deal with this degenerated case,
                            # we take a detour to a polynomial ring
                            # with one variable and then go back.
                            Rt = PolynomialRing(B, 't', 1)
                            It = Rt.ideal([Rt(g) for g in self.gens()])
                            gb = [R(g) for g in It.groebner_basis(
                                algorithm=algorithm,
                                deg_bound=deg_bound, mult_bound=mult_bound,
                                prot=prot, *args, **kwds)]
                    elif (R.term_order().is_global()
                          and isinstance(B, sage.rings.abc.IntegerModRing)
                          and not B.is_field()):
                        verbose("Warning: falling back to very slow toy implementation.", level=0)

                        ch = B.characteristic()
                        R_ZZ = R.change_ring(ZZ)
                        I = R_ZZ.ideal([R_ZZ(f) for f in self.gens()] + [R_ZZ(ch)])
                        gb_ZZ = toy_d_basis.d_basis(I, *args, **kwds)
                        gb = [r for r in (R(f) for f in gb_ZZ) if r]
                    elif R.term_order().is_global():
                        verbose("Warning: falling back to very slow toy implementation.", level=0)
                        gb = toy_buchberger.buchberger_improved(self, *args, **kwds)
                    else:
                        raise TypeError("Local/unknown orderings not supported by 'toy_buchberger' implementation.")

        elif algorithm.startswith('singular:'):
            gb = self._groebner_basis_singular(algorithm[9:], deg_bound=deg_bound, mult_bound=mult_bound, prot=prot, *args, **kwds)
        elif algorithm.startswith('libsingular:'):
            if prot == "sage":
                warn("The libsingular interface does not support prot='sage', reverting to 'prot=True'.")
            gb = self._groebner_basis_libsingular(algorithm[len('libsingular:'):], deg_bound=deg_bound, mult_bound=mult_bound, prot=prot, *args, **kwds)
        elif algorithm.startswith("macaulay2:"):
            gb = self._groebner_basis_macaulay2(strategy=algorithm.split(":")[1], *args, **kwds)
        elif algorithm == 'magma:GroebnerBasis':
            gb = self._groebner_basis_magma(prot=prot, deg_bound=deg_bound, *args, **kwds)
        elif algorithm == 'toy:buchberger':
            gb = toy_buchberger.buchberger(self, *args, **kwds)
        elif algorithm == 'toy:buchberger2':
            gb = toy_buchberger.buchberger_improved(self, *args, **kwds)
        elif algorithm == 'toy:d_basis':
            gb = toy_d_basis.d_basis(self, *args, **kwds)
        elif algorithm.startswith('ginv'):
            if algorithm == 'ginv':
                gb = self._groebner_basis_ginv(*args, **kwds)
            elif ":" in algorithm:
                ginv,alg = algorithm.split(":")
                gb = self._groebner_basis_ginv(algorithm=alg,*args, **kwds)
            else:
                raise NameError("Algorithm '%s' unknown." % algorithm)
        elif algorithm == 'giac:gbasis':
            from sage.libs.giac import groebner_basis as groebner_basis_libgiac
            gb = groebner_basis_libgiac(self, prot=prot, *args, **kwds)
        elif algorithm == 'msolve':
            if self.ring().term_order() != 'degrevlex':
                raise NotImplementedError("msolve only supports the degrevlex order "
                                          "(use transformed_basis())")
            if not (deg_bound is mult_bound is None) or prot:
                raise NotImplementedError("unsupported options for msolve")
            from . import msolve
            return msolve.groebner_basis_degrevlex(self, *args, **kwds)
        else:
            raise NameError("Algorithm '%s' unknown." % algorithm)

        gb = sorted(gb, reverse=True)
        if self.ring().base_ring().is_field():
            _gb = []
            for f in gb:
                if f.lc():
                    _gb.append(f*f.lc()**(-1))
                else:
                    _gb.append(f)
            gb = _gb
        elif self.ring().base_ring() is ZZ:
            if gb[-1].degree() == 0:
                gb = [f % gb[-1] for f in gb[:-1]] + [gb[-1]]

        gb = PolynomialSequence(self.ring(), gb, immutable=True)
        return gb

    def groebner_cover(self):
        r"""
        Compute the Gröbner cover of the ideal, over a field with parameters.

        The Groebner cover is a partition of the space of parameters,
        such that the Groebner basis in each part is given by the same
        expression.

        EXAMPLES::

            sage: F = PolynomialRing(QQ,'a').fraction_field()
            sage: F.inject_variables()
            Defining a
            sage: R.<x,y,z> = F[]
            sage: I = R.ideal([-x+3*y+z-5,2*x+a*z+4,4*x-3*z-1/a])
            sage: I.groebner_cover()
            {Quasi-affine subscheme X - Y of Affine Space of dimension 1 over Rational Field,
               where X is defined by:
                 0
               and Y is defined by:
                 2*a^2 + 3*a: [(2*a^2 + 3*a)*z + (8*a + 1),
                               (12*a^2 + 18*a)*y + (-20*a^2 - 35*a - 2), (4*a + 6)*x + 11],
             Quasi-affine subscheme X - Y of Affine Space of dimension 1 over Rational Field,
               where X is defined by:
                 ...
               and Y is defined by:
                 1: [1],
             Quasi-affine subscheme X - Y of Affine Space of dimension 1 over Rational Field,
               where X is defined by:
                 ...
               and Y is defined by:
                 1: [1]}
        """
        from sage.schemes.affine.affine_space import AffineSpace
        gc = self._groebner_cover()
        F = self.base_ring()
        A = AffineSpace(F.base_ring(), F.ngens(), list(F.gens_dict()))
        result = {}
        ring = F.ring()
        for segment in gc:
            for piece in segment[2]:
                X = A.subscheme([ring(c) for c in piece[0]])
                Y = A.subscheme([ring(c) for c in piece[1][0]])
                for pol in piece[1][1:]:
                    Y = Y.union(A.subscheme([ring(c) for c in pol]))
                result[Y.complement(X)] = segment[1]
        return result

    def change_ring(self, P):
        r"""
        Return the ideal ``I`` in ``P`` spanned by the generators
        `g_1, ..., g_n` of ``self`` as returned by ``self.gens()``.

        INPUT:

        - ``P`` -- a multivariate polynomial ring

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3,order='lex')
            sage: I = sage.rings.ideal.Cyclic(P)
            sage: I
            Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of
            Multivariate Polynomial Ring in x, y, z over Rational Field

        ::

            sage: I.groebner_basis()
            [x + y + z, y^2 + y*z + z^2, z^3 - 1]

        ::

            sage: Q.<x,y,z> = P.change_ring(order='degrevlex'); Q
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: Q.term_order()
            Degree reverse lexicographic term order

        ::

            sage: J = I.change_ring(Q); J
            Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of
            Multivariate Polynomial Ring in x, y, z over Rational Field

        ::

            sage: J.groebner_basis()
            [z^3 - 1, y^2 + y*z + z^2, x + y + z]
        """
        return P.ideal([P(f) for f in self.gens()])

    def subs(self, in_dict=None, **kwds):
        """
        Substitute variables.

        This method substitutes some variables in the polynomials that
        generate the ideal with given values. Variables that are not
        specified in the input remain unchanged.

        INPUT:

        - ``in_dict`` -- (optional) dictionary of inputs

        - ``**kwds`` -- named parameters

        OUTPUT:

        A new ideal with modified generators. If possible, in the same
        polynomial ring. Raises a :exc:`TypeError` if no common
        polynomial ring of the substituted generators can be found.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 2, 'xy')
            sage: I = R.ideal(x^5 + y^5, x^2 + y + x^2*y^2 + 5); I
            Ideal (x^5 + y^5, x^2*y^2 + x^2 + y + 5)
             of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: I.subs(x=y)
            Ideal (2*y^5, y^4 + y^2 + y + 5)
             of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: I.subs({x: y})    # same substitution but with dictionary
            Ideal (2*y^5, y^4 + y^2 + y + 5)
             of Multivariate Polynomial Ring in x, y over Integer Ring

        The new ideal can be in a different ring::

            sage: R.<a,b> = PolynomialRing(QQ, 2)
            sage: S.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal(a^2 + b^2 + a - b + 2); I
            Ideal (a^2 + b^2 + a - b + 2)
             of Multivariate Polynomial Ring in a, b over Rational Field
            sage: I.subs(a=x, b=y)
            Ideal (x^2 + y^2 + x - y + 2)
             of Multivariate Polynomial Ring in x, y over Rational Field

        The resulting ring need not be a multivariate polynomial ring::

            sage: T.<t> = PolynomialRing(QQ)
            sage: I.subs(a=t, b=t)
            Principal ideal (t^2 + 1) of Univariate Polynomial Ring in t over Rational Field
            sage: var("z")                                                              # needs sage.symbolic
            z
            sage: I.subs(a=z, b=z)                                                      # needs sage.symbolic
            Principal ideal (2*z^2 + 2) of Symbolic Ring

        Variables that are not substituted remain unchanged::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal(x^2 + y^2 + x - y + 2); I
            Ideal (x^2 + y^2 + x - y + 2)
             of Multivariate Polynomial Ring in x, y over Rational Field
            sage: I.subs(x=1)
            Ideal (y^2 - y + 4) of Multivariate Polynomial Ring in x, y over Rational Field
        """
        ring = self.ring()
        generators = [f.subs(in_dict, **kwds) for f in self.gens()]
        if not all(gen in ring for gen in generators):
            ring = Sequence(generators).universe()
        try:
            return ring.ideal(generators)
        except AttributeError:
            raise TypeError('Cannot construct an ideal from the substituted generators!')

    def reduce(self, f):
        """
        Reduce an element modulo the reduced Groebner basis for this ideal.
        This returns 0 if and only if the element is in this ideal. In any
        case, this reduction is unique up to monomial orders.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = (x^3 + y, y) * R
            sage: I.reduce(y)
            0
            sage: I.reduce(x^3)
            0
            sage: I.reduce(x - y)
            x

            sage: I = (y^2 - (x^3 + x)) * R
            sage: I.reduce(x^3)
            y^2 - x
            sage: I.reduce(x^6)
            y^4 - 2*x*y^2 + x^2
            sage: (y^2 - x)^2
            y^4 - 2*x*y^2 + x^2

        .. NOTE::

            Requires computation of a Groebner basis, which can be a
            very expensive operation.
        """
        try:
            strat = self._groebner_strategy()
            return strat.normal_form(f)
        except (TypeError, NotImplementedError, ValueError):
            pass

        gb = self.groebner_basis()
        return f.reduce(gb)

    def _contains_(self, f):
        r"""
        Return ``True`` if ``f`` is in this ideal,
        ``False`` otherwise.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = (x^3 + y, y) * R
            sage: x in I  # indirect doctest
            False
            sage: y in I
            True
            sage: x^3 + 2*y in I
            True

        .. NOTE::

            Requires computation of a Groebner basis, which can be a very
            expensive operation.
        """
        return self.reduce(f).is_zero()

    def homogenize(self, var='h'):
        """
        Return homogeneous ideal spanned by the homogeneous polynomials
        generated by homogenizing the generators of this ideal.

        INPUT:

        - ``h`` -- variable name or variable in cover ring (default: ``'h'``)

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(2))
            sage: I = Ideal([x^2*y + z + 1, x + y^2 + 1]); I
            Ideal (x^2*y + z + 1, y^2 + x + 1) of Multivariate
            Polynomial Ring in x, y, z over Finite Field of size 2

        ::

            sage: I.homogenize()
            Ideal (x^2*y + z*h^2 + h^3, y^2 + x*h + h^2) of
            Multivariate Polynomial Ring in x, y, z, h over Finite
            Field of size 2

        ::

            sage: I.homogenize(y)
            Ideal (x^2*y + y^3 + y^2*z, x*y) of Multivariate
            Polynomial Ring in x, y, z over Finite Field of size 2

        ::

            sage: I = Ideal([x^2*y + z^3 + y^2*x, x + y^2 + 1])
            sage: I.homogenize()
            Ideal (x^2*y + x*y^2 + z^3, y^2 + x*h + h^2) of
            Multivariate Polynomial Ring in x, y, z, h over Finite
            Field of size 2
        """
        I = [f.homogenize(var) for f in self.gens()]
        P = max(I, key=lambda x: x.parent().ngens()).parent()
        return P.ideal([P(f) for f in I])

    def is_homogeneous(self):
        r"""
        Return ``True`` if this ideal is spanned by homogeneous
        polynomials, i.e., if it is a homogeneous ideal.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I
            Ideal (x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y +
            2*y*z - y) of Multivariate Polynomial Ring in x, y, z over
            Rational Field

        ::

            sage: I.is_homogeneous()
            False

        ::

            sage: J = I.homogenize()
            sage: J
            Ideal (x + 2*y + 2*z - h, x^2 + 2*y^2 + 2*z^2 - x*h, 2*x*y
            + 2*y*z - y*h) of Multivariate Polynomial Ring in x, y, z,
            h over Rational Field

        ::

            sage: J.is_homogeneous()
            True
        """
        for f in self.gens():
            if not f.is_homogeneous():
                return False
        return True

    def degree_of_semi_regularity(self):
        r"""
        Return the degree of semi-regularity of this ideal under the
        assumption that it is semi-regular.

        Let `\{f_1, ... , f_m\} \subset K[x_1 , ... , x_n]` be
        homogeneous polynomials of degrees `d_1,... ,d_m`
        respectively. This sequence is semi-regular if:

         * `\{f_1, ... , f_m\} \neq K[x_1 , ... , x_n]`

         * for all `1 \leq i \leq m` and `g \in K[x_1,\dots,x_n]`:
           `deg(g \cdot pi ) < D` and
           `g \cdot f_i \in <f_1 , \dots , f_{i-1}>` implies that
           `g \in <f_1, ..., f_{i-1}>` where `D` is the degree of regularity.

        This notion can be extended to affine polynomials by
        considering their homogeneous components of highest degree.

        The degree of regularity of a semi-regular sequence
        `f_1, ...,f_m` of respective degrees `d_1,...,d_m` is given by the
        index of the first nonpositive coefficient of:

            `\sum c_k z^k = \frac{\prod (1 - z^{d_i})}{(1-z)^n}`

        EXAMPLES:

        We consider a homogeneous example::

            sage: n = 8
            sage: K = GF(127)
            sage: P = PolynomialRing(K, n, 'x')
            sage: s = [K.random_element() for _ in range(n)]
            sage: L = []
            sage: for i in range(2 * n):
            ....:     f = P.random_element(degree=2, terms=binomial(n, 2))
            ....:     f -= f(*s)
            ....:     L.append(f.homogenize())
            sage: I = Ideal(L)
            sage: I.degree_of_semi_regularity()
            4

        From this, we expect a Groebner basis computation to reach at
        most degree 4. For homogeneous systems this is equivalent to
        the largest degree in the Groebner basis::

            sage: max(f.degree() for f in I.groebner_basis())
            4

        We increase the number of polynomials and observe a decrease
        the degree of regularity::

            sage: for i in range(2 * n):
            ....:     f = P.random_element(degree=2, terms=binomial(n, 2))
            ....:     f -= f(*s)
            ....:     L.append(f.homogenize())
            sage: I = Ideal(L)
            sage: I.degree_of_semi_regularity()
            3

            sage: max(f.degree() for f in I.groebner_basis())
            3

        The degree of regularity approaches 2 for quadratic systems as
        the number of polynomials approaches `n^2`::

            sage: for i in range((n-4) * n):
            ....:     f = P.random_element(degree=2, terms=binomial(n, 2))
            ....:     f -= f(*s)
            ....:     L.append(f.homogenize())
            sage: I = Ideal(L)
            sage: I.degree_of_semi_regularity()
            2

            sage: max(f.degree() for f in I.groebner_basis())
            2

        .. NOTE::

            It is unknown whether semi-regular sequences
            exist. However, it is expected that random systems are
            semi-regular sequences. For more details about
            semi-regular sequences see [BFS2004]_.
        """
        degs = [f.degree() for f in self.gens() if f != 0]  # we ignore zeroes
        m, n = self.ngens(), len(set(sum([f.variables() for f in self.gens()],())))
        if m <= n:
            raise ValueError("This function requires an overdefined system of polynomials.")

        from sage.misc.misc_c import prod
        from sage.rings.power_series_ring import PowerSeriesRing
        from sage.rings.rational_field import QQ

        z = PowerSeriesRing(QQ, 'z', default_prec=sum(degs)).gen()
        s = prod([1-z**d for d in degs]) / (1-z)**n
        for dreg in range(sum(degs)):
            if s[dreg] <= 0:
                return ZZ(dreg)
        raise ValueError("BUG: Could not compute the degree of semi-regularity")

    def plot(self, *args, **kwds):
        """
        Plot the real zero locus of this principal ideal.

        INPUT:

        - ``self`` -- a principal ideal in 2 variables

        - ``algorithm`` -- set this to 'surf' if you want 'surf' to
          plot the ideal (default: ``None``)

        - ``*args`` -- (optional) tuples ``(variable, minimum, maximum)``
          for plotting dimensions

        - ``**kwds`` -- optional keyword arguments passed on to
          ``implicit_plot``

        EXAMPLES:

        Implicit plotting in 2-d::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal([y^3 - x^2])
            sage: I.plot()                         # cusp                               # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: I = R.ideal([y^2 - x^2 - 1])
            sage: I.plot((x,-3, 3), (y, -2, 2))    # hyperbola                          # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: I = R.ideal([y^2 + x^2*(1/4) - 1])
            sage: I.plot()                         # ellipse                            # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: I = R.ideal([y^2-(x^2-1)*(x-2)])
            sage: I.plot()                         # elliptic curve                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: f = ((x+3)^3 + 2*(x+3)^2 - y^2)*(x^3 - y^2)*((x-3)^3-2*(x-3)^2-y^2)
            sage: I = R.ideal(f)
            sage: I.plot()                         # the Singular logo                  # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal([x - 1])
            sage: I.plot((y, -2, 2))               # vertical line                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        ::

            sage: I = R.ideal([-x^2*y + 1])
            sage: I.plot()                         # blow up                            # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.all import implicit_plot
        from sage.rings.real_mpfr import RR

        K = self.base_ring()
        if not RR.has_coerce_map_from(K):
            raise NotImplementedError("plotting of curves over %s is not implemented yet" % K)

        if not self.is_principal():
            raise TypeError("ideal must be principal")

        f = self.gens()[0]

        variables = sorted(f.parent().gens(), reverse=True)

        if len(variables) == 2 and kwds.get('algorithm','') != 'surf':
            V = [(variables[0], None, None), (variables[1], None, None)]

            if len(args) > 2:
                raise TypeError("Expected up to 2 optional parameters but got %d." % len(args))

            # first check whether user supplied boundaries
            for e in args:
                if not isinstance(e, (tuple, list)) or len(e) != 3:
                    raise TypeError("Optional parameter must be list or tuple or length 3.")
                v,mi,ma = e

                if v not in variables:
                    raise TypeError("Optional parameter must contain variable of ideal generator.")

                vi = variables.index(v)
                V[vi] = v,mi,ma

            # now check whether we should find boundaries
            for var_index in range(2):
                if V[var_index][1] is None:
                    v, mi, ma = variables[var_index], -10, 10
                    for i in range(mi, ma):
                        poly = f.subs({v:i}).univariate_polynomial().change_ring(RR)
                        if not poly or len(poly.roots()) > 0:
                            mi = i - 1
                            break

                    for i in range(ma, mi, -1):
                        poly = f.subs({v:i}).univariate_polynomial().change_ring(RR)
                        if not poly or len(poly.roots()) > 0:
                            ma = i + 1
                            break
                    V[var_index] = variables[var_index], mi, ma

            kwds.setdefault("plot_points",200)
            kwds.pop('algorithm', '')
            return implicit_plot(f, V[0], V[1], **kwds)

        elif len(variables) == 3 or kwds.get('algorithm','') == 'surf':
            MPolynomialIdeal_singular_repr.plot(self, kwds.get("singular",singular_default))
        else:
            raise TypeError("Ideal generator may not have either 2 or 3 variables.")

    def random_element(self, degree, compute_gb=False, *args, **kwds):
        r"""
        Return a random element in this ideal as `r = \sum h_i·f_i`.

        INPUT:

        - ``compute_gb`` -- if ``True`` then a Gröbner basis is computed first
          and `f_i` are the elements in the Gröbner basis. Otherwise whatever
          basis is returned by ``self.gens()`` is used.

        - ``*args`` and ``**kwds`` are passed to ``R.random_element()`` with
          ``R = self.ring()``.

        EXAMPLES:

        We compute a uniformly random element up to the provided degree. ::

            sage: P.<x,y,z> = GF(127)[]
            sage: I = sage.rings.ideal.Katsura(P)
            sage: f = I.random_element(degree=4, compute_gb=True, terms=infinity)
            sage: f.degree() <= 4
            True
            sage: len(list(f)) <= 35
            True

        Note that sampling uniformly at random from the ideal at some large enough degree is
        equivalent to computing a Gröbner basis. We give an example showing how to compute a Gröbner
        basis if we can sample uniformly at random from an ideal::

            sage: n = 3; d = 4
            sage: P = PolynomialRing(GF(127), n, 'x')
            sage: I = sage.rings.ideal.Cyclic(P)

        1. We sample `n^d` uniformly random elements in the ideal::

               sage: F = Sequence(I.random_element(degree=d, compute_gb=True,
               ....:                               terms=infinity)
               ....:              for _ in range(n^d))

        2. We linearize and compute the echelon form::

               sage: A, v = F.coefficients_monomials()
               sage: A.echelonize()

        3. The result is the desired Gröbner basis::

               sage: G = Sequence((A * v).list())
               sage: G.is_groebner()
               True
               sage: Ideal(G) == I
               True

        We return some element in the ideal with no guarantee on the distribution::

            sage: # needs sage.rings.finite_rings
            sage: P = PolynomialRing(GF(127), 10, 'x')
            sage: I = sage.rings.ideal.Katsura(P)
            sage: f = I.random_element(degree=3)
            sage: f  # random
            -25*x0^2*x1 + 14*x1^3 + 57*x0*x1*x2 + ... + 19*x7*x9 + 40*x8*x9 + 49*x1
            sage: f.degree()
            3

        We show that the default method does not sample uniformly at random from the ideal::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,y,z> = GF(127)[]
            sage: G = Sequence([x + 7, y - 2, z + 110])
            sage: I = Ideal([sum(P.random_element() * g for g in G)
            ....:            for _ in range(4)])
            sage: all(I.random_element(degree=1) == 0 for _ in range(100))
            True

        If ``degree`` equals the degree of the generators, a random linear
        combination of the generators is returned::

            sage: P.<x,y> = QQ[]
            sage: I = P.ideal([x^2,y^2])
            sage: set_random_seed(5)
            sage: I.random_element(degree=2)
            -2*x^2 + 2*y^2
        """
        if compute_gb:
            gens = self.groebner_basis()
        else:
            gens = self.basis

        R = self.ring()

        r = R.zero()

        for f in gens:
            d = degree - f.degree()
            if d >= 0:
                h = R.random_element(degree=d, *args, **kwds)
                r += h*f
        return r

    @require_field
    def weil_restriction(self):
        r"""
        Compute the Weil restriction of this ideal over some extension
        field. If the field is a finite field, then this computes
        the Weil restriction to the prime subfield.

        A Weil restriction of scalars - denoted `Res_{L/k}` - is a
        functor which, for any finite extension of fields `L/k` and
        any algebraic variety `X` over `L`, produces another
        corresponding variety `Res_{L/k}(X)`, defined over `k`. It is
        useful for reducing questions about varieties over large
        fields to questions about more complicated varieties over
        smaller fields.

        This function does not compute this Weil restriction directly
        but computes on generating sets of polynomial ideals:

        Let `d` be the degree of the field extension `L/k`, let `a` a
        generator of `L/k` and `p` the minimal polynomial of
        `L/k`. Denote this ideal by `I`.

        Specifically, this function first maps each variable `x` to
        its representation over `k`: `\sum_{i=0}^{d-1} a^i x_i`. Then
        each generator of `I` is evaluated over these representations
        and reduced modulo the minimal polynomial `p`. The result is
        interpreted as a univariate polynomial in `a` and its
        coefficients are the new generators of the returned ideal.

        If the input and the output ideals are radical, this is
        equivalent to the statement about algebraic varieties above.

        OUTPUT: :class:`MPolynomialIdeal`

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(2^2)
            sage: P.<x,y> = PolynomialRing(k, 2)
            sage: I = Ideal([x*y + 1, a*x + 1])
            sage: I.variety()
            [{y: a, x: a + 1}]
            sage: J = I.weil_restriction()
            sage: J
            Ideal (x0*y0 + x1*y1 + 1, x1*y0 + x0*y1 + x1*y1, x1 + 1, x0 + x1) of
             Multivariate Polynomial Ring in x0, x1, y0, y1 over Finite Field of size 2
            sage: J += sage.rings.ideal.FieldIdeal(J.ring())  # ensure radical ideal
            sage: J.variety()
            [{y1: 1, y0: 0, x1: 1, x0: 1}]
            sage: J.weil_restriction()  # returns J
            Ideal (x0*y0 + x1*y1 + 1, x1*y0 + x0*y1 + x1*y1, x1 + 1, x0 + x1,
             x0^2 + x0, x1^2 + x1, y0^2 + y0, y1^2 + y1) of Multivariate
             Polynomial Ring in x0, x1, y0, y1 over Finite Field of size 2

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(3^5)
            sage: P.<x,y,z> = PolynomialRing(k)
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I.dimension()
            0
            sage: I.variety()
             [{z: 0, y: 0, x: 1}]
            sage: J = I.weil_restriction(); J
            Ideal (x0 - y0 - z0 - 1,
                   x1 - y1 - z1, x2 - y2 - z2, x3 - y3 - z3, x4 - y4 - z4,
                   x0^2 + x2*x3 + x1*x4 - y0^2 - y2*y3 - y1*y4 - z0^2 - z2*z3 - z1*z4 - x0,
                   -x0*x1 - x2*x3 - x3^2 - x1*x4 + x2*x4 + y0*y1 + y2*y3
                        + y3^2 + y1*y4 - y2*y4 + z0*z1 + z2*z3 + z3^2 + z1*z4 - z2*z4 - x1,
                   x1^2 - x0*x2 + x3^2 - x2*x4 + x3*x4 - y1^2 + y0*y2
                         - y3^2 + y2*y4 - y3*y4 - z1^2 + z0*z2 - z3^2 + z2*z4 - z3*z4 - x2,
                   -x1*x2 - x0*x3 - x3*x4 - x4^2
                        + y1*y2 + y0*y3 + y3*y4 + y4^2 + z1*z2 + z0*z3 + z3*z4 + z4^2 - x3,
                   x2^2 - x1*x3 - x0*x4 + x4^2 - y2^2
                                 + y1*y3 + y0*y4 - y4^2 - z2^2 + z1*z3 + z0*z4 - z4^2 - x4,
                   -x0*y0 + x4*y1 + x3*y2 + x2*y3
                                      + x1*y4 - y0*z0 + y4*z1 + y3*z2 + y2*z3 + y1*z4 - y0,
                   -x1*y0 - x0*y1 - x4*y1 - x3*y2 + x4*y2 - x2*y3 + x3*y3
                                      - x1*y4 + x2*y4 - y1*z0 - y0*z1 - y4*z1 - y3*z2
                                              + y4*z2 - y2*z3 + y3*z3 - y1*z4 + y2*z4 - y1,
                   -x2*y0 - x1*y1 - x0*y2 - x4*y2 - x3*y3 + x4*y3 - x2*y4 + x3*y4
                      - y2*z0 - y1*z1 - y0*z2 - y4*z2 - y3*z3 + y4*z3 - y2*z4 + y3*z4 - y2,
                   -x3*y0 - x2*y1 - x1*y2 - x0*y3 - x4*y3 - x3*y4 + x4*y4
                              - y3*z0 - y2*z1 - y1*z2 - y0*z3 - y4*z3 - y3*z4 + y4*z4 - y3,
                   -x4*y0 - x3*y1 - x2*y2 - x1*y3 - x0*y4 - x4*y4
                                      - y4*z0 - y3*z1 - y2*z2 - y1*z3 - y0*z4 - y4*z4 - y4)
             of Multivariate Polynomial Ring in x0, x1, x2, x3, x4, y0, y1, y2, y3, y4,
             z0, z1, z2, z3, z4 over Finite Field of size 3
            sage: J += sage.rings.ideal.FieldIdeal(J.ring())  # ensure radical ideal
            sage: from sage.doctest.fixtures import reproducible_repr
            sage: print(reproducible_repr(J.variety()))
            [{x0: 1, x1: 0, x2: 0, x3: 0, x4: 0,
              y0: 0, y1: 0, y2: 0, y3: 0, y4: 0,
              z0: 0, z1: 0, z2: 0, z3: 0, z4: 0}]


        Weil restrictions are often used to study elliptic curves over
        extension fields so we give a simple example involving those::

            sage: K.<a> = QuadraticField(1/3)                                           # needs sage.rings.number_field
            sage: E = EllipticCurve(K, [1,2,3,4,5])                                     # needs sage.rings.number_field

        We pick a point on ``E``::

            sage: p = E.lift_x(1); p                                                    # needs sage.rings.number_field
            (1 : -6 : 1)

            sage: I = E.defining_ideal(); I                                             # needs sage.rings.number_field
            Ideal (-x^3 - 2*x^2*z + x*y*z + y^2*z - 4*x*z^2 + 3*y*z^2 - 5*z^3)
             of Multivariate Polynomial Ring in x, y, z
              over Number Field in a with defining polynomial x^2 - 1/3
               with a = 0.5773502691896258?

        Of course, the point ``p`` is a root of all generators of ``I``::

            sage: I.subs(x=1, y=2, z=1)                                                 # needs sage.rings.number_field
            Ideal (0) of Multivariate Polynomial Ring in x, y, z
             over Number Field in a with defining polynomial x^2 - 1/3
              with a = 0.5773502691896258?

        ``I`` is also radical::

            sage: I.radical() == I                                                      # needs sage.rings.number_field
            True

        So we compute its Weil restriction::

            sage: J = I.weil_restriction(); J                                           # needs sage.rings.number_field
            Ideal (-x0^3 - x0*x1^2 - 2*x0^2*z0 - 2/3*x1^2*z0 + x0*y0*z0 + y0^2*z0
                     + 1/3*x1*y1*z0 + 1/3*y1^2*z0 - 4*x0*z0^2 + 3*y0*z0^2 - 5*z0^3
                     - 4/3*x0*x1*z1 + 1/3*x1*y0*z1 + 1/3*x0*y1*z1 + 2/3*y0*y1*z1
                     - 8/3*x1*z0*z1 + 2*y1*z0*z1 - 4/3*x0*z1^2 + y0*z1^2 - 5*z0*z1^2,
                   -3*x0^2*x1 - 1/3*x1^3 - 4*x0*x1*z0 + x1*y0*z0 + x0*y1*z0
                     + 2*y0*y1*z0 - 4*x1*z0^2 + 3*y1*z0^2 - 2*x0^2*z1 - 2/3*x1^2*z1
                     + x0*y0*z1 + y0^2*z1 + 1/3*x1*y1*z1 + 1/3*y1^2*z1 - 8*x0*z0*z1
                     + 6*y0*z0*z1 - 15*z0^2*z1 - 4/3*x1*z1^2 + y1*z1^2 - 5/3*z1^3)
            of Multivariate Polynomial Ring in x0, x1, y0, y1, z0, z1 over Rational Field

        We can check that the point ``p`` is still a root of all generators of ``J``::

            sage: J.subs(x0=1, y0=2, z0=1, x1=0, y1=0, z1=0)                            # needs sage.rings.number_field
            Ideal (0, 0) of Multivariate Polynomial Ring in x0, x1, y0, y1, z0, z1
             over Rational Field

        Example for relative number fields::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^5 - 2)
            sage: R.<x> = K[]
            sage: L.<v> = K.extension(x^2 + 1)
            sage: S.<x,y> = L[]
            sage: I = S.ideal([y^2 - x^3 - 1])
            sage: I.weil_restriction()
            Ideal (-x0^3 + 3*x0*x1^2 + y0^2 - y1^2 - 1, -3*x0^2*x1 + x1^3 + 2*y0*y1) of
             Multivariate Polynomial Ring in x0, x1, y0, y1
              over Number Field in w with defining polynomial x^5 - 2

        .. NOTE::

            Based on a Singular implementation by Michael Brickenstein
        """
        from sage.rings.polynomial.polynomial_ring_constructor import \
            PolynomialRing

        R = self.ring()
        nvars = R.ngens()
        L = R.base_ring()
        if L.is_finite():
            k = L.prime_subfield()
            d = L.degree()
            poly = L.polynomial()
        else:
            k = L.base_field()
            d = L.relative_degree()
            poly = L.relative_polynomial()

        if d == 1:
            return self

        helper = PolynomialRing(k, nvars + 1, (L.variable_name(),) + R.variable_names(), order='lex')
        myminpoly = poly.subs(helper.gen(0))

        l = [helper(str(f)) for f in self.gens()]
        r = myminpoly.degree()
        intermediate_ring = PolynomialRing(k, nvars * r + 1, 'x')
        a = intermediate_ring.gen(0)

        # map e.g. x -> a^2*x_2 + a*x_1 + x_0, where x_0,..,x_2
        # represent the components of x if viewed as a vector in k^r
        map_ideal = [a]

        variables = iter(intermediate_ring.gens()[1:])
        map_ideal.extend(sum([a**i * next(variables) for i in range(r)])
                         for _ in range(nvars))

        myminpoly = myminpoly(*map_ideal)
        l = [f(*map_ideal).reduce([myminpoly]) for f in l]

        result = []
        # split e.g. a^2*x2+a*x1+x0 to x0,x1,x2
        for f in l:
            t = []
            for i in reversed(range(r)):
                g = f.coefficient(a**i)
                f = f - a**i * g
                t.append(g)
            result += reversed(t)

        # eliminate parameter
        new_var_names = [str(var) + "%d" % i for var in R.gens()
                         for i in range(r)]

        result_ring = PolynomialRing(k, nvars * r, new_var_names)

        map_ideal = (0,) + result_ring.gens()
        result = [h(*map_ideal) for h in result]

        return result_ring.ideal(result)


class MPolynomialIdeal_quotient(MPolynomialIdeal):
    r"""
    An ideal in a quotient of a multivariate polynomial ring.

    EXAMPLES::

        sage: Q.<x,y,z,w> = QQ['x,y,z,w'].quotient(['x*y-z^2', 'y^2-w^2'])
        sage: I = ideal(x + y^2 + z - 1)
        sage: I
        Ideal (w^2 + x + z - 1) of Quotient
         of Multivariate Polynomial Ring in x, y, z, w over Rational Field
         by the ideal (x*y - z^2, y^2 - w^2)
    """
    def reduce(self, f):
        r"""
        Reduce an element modulo a Gröbner basis for this ideal.
        This returns 0 if and only if the element is in this ideal. In any
        case, this reduction is unique up to monomial orders.

        EXAMPLES::

            sage: R.<T,U,V,W,X,Y,Z> = PolynomialRing(QQ, order='lex')
            sage: I = R.ideal([T^2 + U^2 - 1, V^2 + W^2 - 1, X^2 + Y^2 + Z^2 - 1])
            sage: Q.<t,u,v,w,x,y,z> = R.quotient(I)
            sage: J = Q.ideal([u*v - x, u*w - y, t - z])
            sage: J.reduce(t^2 - z^2)
            0
            sage: J.reduce(u^2)
            -z^2 + 1
            sage: t^2 - z^2 in J
            True
            sage: u^2 in J
            False
        """
        Q = self.ring()
        gb = self.groebner_basis()
        # In quotient rings, gb is not a Gröbner basis of self, but gb0 is
        # a (possibly non-reduced) Gröbner basis of the preimage of self in
        # the cover ring (see :issue:`33217`). We only use Gröbner bases of
        # pre-existing ideals to potentially take advantage of caching.
        gb0 = Q.defining_ideal().groebner_basis() + [g.lift() for g in gb]
        f0 = f.lift().reduce(gb0)
        return Q._element_constructor_from_element_class(f0, reduce=False)

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- an ideal in a quotient of a polynomial ring

        OUTPUT: boolean

        TESTS::

            sage: R.<T,U,V,W,X,Y,Z> = PolynomialRing(QQ, order='lex')
            sage: Q.<t,u,v,w,x,y,z> = R.quotient([T^2+U^2-1, V^2+W^2-1, X^2+Y^2+Z^2-1])
            sage: J = Q.ideal([u*v-x, u*w-y, t-z])
            sage: I1 = J + Q.ideal(t^2-z^2)
            sage: I1 <= J, I1 < J, I1 == J, I1 != J
            (True, False, True, False)
            sage: I2 = J + Q.ideal(t-z^2)
            sage: J <= I2, J < I2, J > I2, J >= I2
            (True, True, False, False)

        The ideals must belong to the same quotient ring::

            sage: J0 = R.ideal([g.lift() for g in J.gens()])
            sage: J0 <= J
            Traceback (most recent call last):
            ...
            AttributeError:...
            sage: J <= J0
            Traceback (most recent call last):
            ...
            TypeError: '<=' not supported...
        """
        # The implementation largely follows the superclass, but for simplicity
        # does not deal with different rings or term orders.
        if not isinstance(other, MPolynomialIdeal_quotient):
            return NotImplemented

        if self is other:
            return rich_to_bool(op, 0)

        # comparison for >= and > : swap the arguments
        if op == op_GE:
            return other.__richcmp__(self, op_LE)
        elif op == op_GT:
            return other.__richcmp__(self, op_LT)

        if other.ring() != self.ring():
            return NotImplemented

        s_gens = self.gens()
        o_gens = other.gens()
        try:
            if (s_gens == o_gens) or (set(s_gens) == set(o_gens)):
                # the first clause works in the non-hashable case
                return rich_to_bool(op, 0)
        except TypeError:
            pass

        # comparison for <= and == and != and <
        if op in [op_LE, op_EQ, op_NE, op_LT]:
            contained = all(f in other for f in s_gens)
            if op == op_LE:
                return contained
            contains = all(g in self for g in o_gens)
            if op == op_EQ:
                return contained and contains
            elif op == op_NE:
                return not (contained and contains)
            else:  # remaining case <
                return contained and not contains
