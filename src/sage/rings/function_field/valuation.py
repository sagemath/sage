r"""
Discrete valuations on function fields

AUTHORS:

- Julian Rüth (2016-10-16): initial version

EXAMPLES:

We can create classical valuations that correspond to finite and infinite
places on a rational function field::

    sage: K.<x> = FunctionField(QQ)
    sage: v = K.valuation(1); v
    (x - 1)-adic valuation
    sage: v = K.valuation(x^2 + 1); v
    (x^2 + 1)-adic valuation
    sage: v = K.valuation(1/x); v
    Valuation at the infinite place

Note that we can also specify valuations which do not correspond to a place of
the function field::

    sage: R.<x> = QQ[]
    sage: w = valuations.GaussValuation(R, QQ.valuation(2))
    sage: v = K.valuation(w); v
    2-adic valuation

Valuations on a rational function field can then be extended to finite
extensions::

    sage: v = K.valuation(x - 1); v
    (x - 1)-adic valuation
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)                                                  # needs sage.rings.function_field
    sage: w = v.extensions(L); w                                                        # needs sage.rings.function_field
    [[ (x - 1)-adic valuation, v(y + 1) = 1 ]-adic valuation,
     [ (x - 1)-adic valuation, v(y - 1) = 1 ]-adic valuation]

TESTS:

Run test suite for classical places over rational function fields::

    sage: K.<x> = FunctionField(QQ)
    sage: v = K.valuation(1)
    sage: TestSuite(v).run(max_runs=100) # long time

    sage: v = K.valuation(x^2 + 1)
    sage: TestSuite(v).run(max_runs=100) # long time

    sage: v = K.valuation(1/x)
    sage: TestSuite(v).run(max_runs=100) # long time

Run test suite over classical places of finite extensions::

    sage: K.<x> = FunctionField(QQ)
    sage: v = K.valuation(x - 1)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)                                                  # needs sage.rings.function_field
    sage: ws = v.extensions(L)                                                          # needs sage.rings.function_field
    sage: for w in ws: TestSuite(w).run(max_runs=100)   # long time                     # needs sage.rings.function_field

Run test suite for valuations that do not correspond to a classical place::

    sage: K.<x> = FunctionField(QQ)
    sage: R.<x> = QQ[]
    sage: v = GaussValuation(R, QQ.valuation(2))
    sage: w = K.valuation(v)
    sage: TestSuite(w).run() # long time

Run test suite for a non-classical valuation that does not correspond to an
affinoid contained in the unit disk::

    sage: w = K.valuation((w, K.hom(K.gen()/2), K.hom(2*K.gen()))); w
    2-adic valuation (in Rational function field in x over Rational Field after x |--> 1/2*x)
    sage: TestSuite(w).run() # long time

Run test suite for some other classical places over large ground fields::

    sage: K.<t> = FunctionField(GF(3))
    sage: M.<x> = FunctionField(K)
    sage: v = M.valuation(x^3 - t)
    sage: TestSuite(v).run(max_runs=10)         # long time

Run test suite for extensions over the infinite place::

    sage: K.<x> = FunctionField(QQ)
    sage: v = K.valuation(1/x)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - 1/(x^2 + 1))                                        # needs sage.rings.function_field
    sage: w = v.extensions(L)                                                           # needs sage.rings.function_field
    sage: TestSuite(w).run()                    # long time                             # needs sage.rings.function_field

Run test suite for a valuation with `v(1/x) > 0` which does not come from a
classical valuation of the infinite place::

    sage: K.<x> = FunctionField(QQ)
    sage: R.<x> = QQ[]
    sage: w = GaussValuation(R, QQ.valuation(2)).augmentation(x, 1)
    sage: w = K.valuation(w)
    sage: v = K.valuation((w, K.hom([~K.gen()]), K.hom([~K.gen()])))
    sage: TestSuite(v).run() # long time

Run test suite for extensions which come from the splitting in the base field::

    sage: K.<x> = FunctionField(QQ)
    sage: v = K.valuation(x^2 + 1)
    sage: L.<x> = FunctionField(GaussianIntegers().fraction_field())
    sage: ws = v.extensions(L)                                                          # needs sage.rings.function_field
    sage: for w in ws: TestSuite(w).run(max_runs=100)   # long time                     # needs sage.rings.function_field

Run test suite for a finite place with residual degree and ramification::

    sage: K.<t> = FunctionField(GF(3))
    sage: L.<x> = FunctionField(K)
    sage: v = L.valuation(x^6 - t)
    sage: TestSuite(v).run(max_runs=10) # long time

Run test suite for a valuation which is backed by limit valuation::

    sage: K.<x> = FunctionField(QQ)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - (x^2 + x + 1))
    sage: v = K.valuation(x - 1)
    sage: w = v.extension(L)                                                            # needs sage.rings.function_field
    sage: TestSuite(w).run()                    # long time                             # needs sage.rings.function_field

Run test suite for a valuation which sends an element to `-\infty`::

    sage: R.<x> = QQ[]
    sage: v = GaussValuation(QQ['x'], QQ.valuation(2)).augmentation(x, infinity)
    sage: K.<x> = FunctionField(QQ)
    sage: w = K.valuation(v)
    sage: TestSuite(w).run() # long time

REFERENCES:

An overview of some computational tools relating to valuations on function
fields can be found in Section 4.6 of [Rüt2014]_. Most of this was originally
developed for number fields in [Mac1936I]_ and [Mac1936II]_.

"""

# ****************************************************************************
#       Copyright (C) 2016-2018 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.factory import UniqueFactory
from sage.rings.rational_field import QQ
from sage.misc.cachefunc import cached_method

from sage.rings.valuation.valuation import DiscreteValuation, DiscretePseudoValuation, InfiniteDiscretePseudoValuation, NegativeInfiniteDiscretePseudoValuation
from sage.rings.valuation.trivial_valuation import TrivialValuation
from sage.rings.valuation.mapped_valuation import FiniteExtensionFromLimitValuation, MappedValuation_base

class FunctionFieldValuationFactory(UniqueFactory):
    r"""
    Create a valuation on ``domain`` corresponding to ``prime``.

    INPUT:

    - ``domain`` -- a function field

    - ``prime`` -- a place of the function field, a valuation on a subring, or
      a valuation on another function field together with information for
      isomorphisms to and from that function field

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: v = K.valuation(1); v  # indirect doctest
        (x - 1)-adic valuation
        sage: v(x)
        0
        sage: v(x - 1)
        1

    See :meth:`sage.rings.function_field.function_field.FunctionField.valuation` for further examples.

    """
    def create_key_and_extra_args(self, domain, prime):
        r"""
        Create a unique key which identifies the valuation given by ``prime``
        on ``domain``.

        TESTS:

        We specify a valuation on a function field by two different means and
        get the same object::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x - 1)  # indirect doctest

            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, valuations.TrivialValuation(QQ)).augmentation(x - 1, 1)
            sage: K.valuation(w) is v
            True

        The normalization is, however, not smart enough, to unwrap
        substitutions that turn out to be trivial::

            sage: w = GaussValuation(R, QQ.valuation(2))
            sage: w = K.valuation(w)
            sage: w is K.valuation((w, K.hom([~K.gen()]), K.hom([~K.gen()])))
            False

        """
        from sage.categories.function_fields import FunctionFields
        if domain not in FunctionFields():
            raise ValueError("Domain must be a function field.")

        if isinstance(prime, tuple):
            if len(prime) == 3:
                # prime is a triple of a valuation on another function field with
                # isomorphism information
                return self.create_key_and_extra_args_from_valuation_on_isomorphic_field(domain, prime[0], prime[1], prime[2])

        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        if prime.parent() is DiscretePseudoValuationSpace(domain):
            # prime is already a valuation of the requested domain
            # if we returned (domain, prime), we would break caching
            # because this element has been created from a different key
            # Instead, we return the key that was used to create prime
            # so the caller gets back a correctly cached version of prime
            if not hasattr(prime, "_factory_data"):
                raise NotImplementedError("Valuations on function fields must be unique and come out of the FunctionFieldValuation factory but %r has been created by other means" % (prime,))
            return prime._factory_data[2], {}

        if prime in domain:
            # prime defines a place
            return self.create_key_and_extra_args_from_place(domain, prime)
        if prime.parent() is DiscretePseudoValuationSpace(domain._ring):
            # prime is a discrete (pseudo-)valuation on the polynomial ring
            # that the domain is constructed from
            return self.create_key_and_extra_args_from_valuation(domain, prime)
        if domain.base_field() is not domain:
            # prime might define a valuation on a subring of domain and have a
            # unique extension to domain
            base_valuation = domain.base_field().valuation(prime)
            return self.create_key_and_extra_args_from_valuation(domain, base_valuation)
        from sage.rings.ideal import is_Ideal
        if is_Ideal(prime):
            raise NotImplementedError("a place cannot be given by an ideal yet")

        raise NotImplementedError("argument must be a place or a pseudo-valuation on a supported subring but %r does not satisfy this for the domain %r" % (prime, domain))

    def create_key_and_extra_args_from_place(self, domain, generator):
        r"""
        Create a unique key which identifies the valuation at the place
        specified by ``generator``.

        TESTS:

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(1/x)  # indirect doctest

        """
        if generator not in domain.base_field():
            raise NotImplementedError("a place must be defined over a rational function field")

        if domain.base_field() is not domain:
            # if this is an extension field, construct the unique place over
            # the place on the subfield
            return self.create_key_and_extra_args(domain, domain.base_field().valuation(generator))

        if generator in domain.constant_base_field():
            # generator is a constant, we associate to it the place which
            # corresponds to the polynomial (x - generator)
            return self.create_key_and_extra_args(domain, domain.gen() - generator)

        if generator in domain._ring:
            # generator is a polynomial
            generator = domain._ring(generator)
            if not generator.is_monic():
                raise ValueError("place must be defined by a monic polynomial but %r is not monic" % (generator,))
            if not generator.is_irreducible():
                raise ValueError("place must be defined by an irreducible polynomial but %r factors over %r" % (generator, domain._ring))
            # we construct the corresponding valuation on the polynomial ring
            # with v(generator) = 1
            from sage.rings.valuation.gauss_valuation import GaussValuation
            valuation = GaussValuation(domain._ring, TrivialValuation(domain.constant_base_field())).augmentation(generator, 1)
            return self.create_key_and_extra_args(domain, valuation)
        elif generator == ~domain.gen():
            # generator is 1/x, the infinite place
            return (domain, (domain.valuation(domain.gen()), domain.hom(~domain.gen()), domain.hom(~domain.gen()))), {}
        else:
            raise ValueError("a place must be given by an irreducible polynomial or the inverse of the generator; %r does not define a place over %r" % (generator, domain))

    def create_key_and_extra_args_from_valuation(self, domain, valuation):
        r"""
        Create a unique key which identifies the valuation which extends
        ``valuation``.

        TESTS:

            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, valuations.TrivialValuation(QQ)).augmentation(x - 1, 1)
            sage: v = K.valuation(w) # indirect doctest

        Check that :trac:`25294` has been resolved::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 1/x^3*y + 2/x^4)                            # needs sage.rings.function_field
            sage: v = K.valuation(x)                                                    # needs sage.rings.function_field
            sage: v.extensions(L)                                                       # needs sage.rings.function_field
            [[ (x)-adic valuation, v(y) = 1 ]-adic valuation
               (in Function field in y defined by y^3 + x*y + 2*x^2 after y |--> 1/x^2*y),
             [ (x)-adic valuation, v(y) = 1/2 ]-adic valuation
               (in Function field in y defined by y^3 + x*y + 2*x^2 after y |--> 1/x^2*y)]

        """
        # this should have been handled by create_key already
        assert valuation.domain() is not domain

        if valuation.domain() is domain._ring:
            if domain.base_field() is not domain:
                vK = valuation.restriction(valuation.domain().base_ring())
                if vK.domain() is not domain.base_field():
                    raise ValueError("valuation must extend a valuation on the base field but %r extends %r whose domain is not %r" % (valuation, vK, domain.base_field()))
                # Valuation is an approximant that describes a single valuation
                # on domain.
                # For uniqueness of valuations (which provides better caching
                # and easier pickling) we need to find a normal form of
                # valuation, i.e., the smallest approximant that describes this
                # valuation
                approximants = vK.mac_lane_approximants(domain.polynomial(), require_incomparability=True)
                approximant = vK.mac_lane_approximant(domain.polynomial(), valuation, approximants)
                return (domain, approximant), {'approximants': approximants}
            else:
                # on a rational function field K(x), any valuation on K[x] that
                # does not have an element with valuation -infty extends to a
                # pseudo-valuation on K(x)
                if valuation.is_negative_pseudo_valuation():
                    raise ValueError("there must not be an element of valuation -Infinity in the domain of valuation %r" % (valuation,))
                return (domain, valuation), {}

        if valuation.domain().is_subring(domain.base_field()):
            # valuation is defined on a subring of this function field, try to lift it
            return self.create_key_and_extra_args(domain, valuation.extension(domain))

        raise NotImplementedError("extension of valuation from %r to %r not implemented yet" % (valuation.domain(), domain))

    def create_key_and_extra_args_from_valuation_on_isomorphic_field(self, domain, valuation, to_valuation_domain, from_valuation_domain):
        r"""
        Create a unique key which identifies the valuation which is
        ``valuation`` after mapping through ``to_valuation_domain``.

        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)                                    # needs sage.rings.function_field
            sage: v = K.valuation(1/x)                                                  # needs sage.rings.function_field
            sage: w = v.extension(L) # indirect doctest                                 # needs sage.rings.function_field

        """
        from sage.categories.function_fields import FunctionFields
        if valuation.domain() not in FunctionFields():
            raise ValueError("valuation must be defined over an isomorphic function field but %r is not a function field" % (valuation.domain(),))

        from sage.categories.homset import Hom
        if to_valuation_domain not in Hom(domain, valuation.domain()):
            raise ValueError("to_valuation_domain must map from %r to %r but %r maps from %r to %r" % (domain, valuation.domain(), to_valuation_domain, to_valuation_domain.domain(), to_valuation_domain.codomain()))
        if from_valuation_domain not in Hom(valuation.domain(), domain):
            raise ValueError("from_valuation_domain must map from %r to %r but %r maps from %r to %r" % (valuation.domain(), domain, from_valuation_domain, from_valuation_domain.domain(), from_valuation_domain.codomain()))

        if domain is domain.base():
            if valuation.domain() is not valuation.domain().base() or valuation.domain().constant_base_field() != domain.constant_base_field():
                raise NotImplementedError("maps must be isomorphisms with a rational function field over the same base field, not with %r" % (valuation.domain(),))
            if domain != valuation.domain():
                # make it harder to create different representations of the same valuation
                # (nothing bad happens if we did, but >= and <= are only implemented when this is the case.)
                raise NotImplementedError("domain and valuation.domain() must be the same rational function field but %r is not %r" % (domain, valuation.domain()))
        else:
            if domain.base() is not valuation.domain().base():
                raise NotImplementedError("domain and valuation.domain() must have the same base field but %r is not %r" % (domain.base(), valuation.domain().base()))
            if to_valuation_domain != domain.hom([to_valuation_domain(domain.gen())]):
                raise NotImplementedError("to_valuation_domain must be trivial on the base fields but %r is not %r" % (to_valuation_domain, domain.hom([to_valuation_domain(domain.gen())])))
            if from_valuation_domain != valuation.domain().hom([from_valuation_domain(valuation.domain().gen())]):
                raise NotImplementedError("from_valuation_domain must be trivial on the base fields but %r is not %r" % (from_valuation_domain, valuation.domain().hom([from_valuation_domain(valuation.domain().gen())])))
            if to_valuation_domain(domain.gen()) == valuation.domain().gen():
                raise NotImplementedError("to_valuation_domain seems to be trivial but trivial maps would currently break partial orders of valuations")

        if from_valuation_domain(to_valuation_domain(domain.gen())) != domain.gen():
            # only a necessary condition
            raise ValueError("to_valuation_domain and from_valuation_domain are not inverses of each other")

        return (domain, (valuation, to_valuation_domain, from_valuation_domain)), {}

    def create_object(self, version, key, **extra_args):
        r"""
        Create the valuation specified by ``key``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: w = valuations.GaussValuation(R, QQ.valuation(2))
            sage: v = K.valuation(w); v  # indirect doctest
            2-adic valuation

        """
        domain, valuation = key
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(domain)

        if isinstance(valuation, tuple) and len(valuation) == 3:
            valuation, to_valuation_domain, from_valuation_domain = valuation
            if domain is domain.base() and valuation.domain() is valuation.domain().base():
                if valuation == valuation.domain().valuation(valuation.domain().gen()):
                    if to_valuation_domain != domain.hom([~valuation.domain().gen()]) or from_valuation_domain != valuation.domain().hom([~domain.gen()]):
                        raise ValueError("the only allowed automorphism for classical valuations is the automorphism x |--> 1/x")
                    # valuation on the rational function field after x |--> 1/x,
                    # i.e., the classical valuation at infinity
                    return parent.__make_element_class__(InfiniteRationalFunctionFieldValuation)(parent)

                from sage.structure.dynamic_class import dynamic_class
                clazz = RationalFunctionFieldMappedValuation
                if valuation.is_discrete_valuation():
                    clazz = dynamic_class("RationalFunctionFieldMappedValuation_discrete", (clazz, DiscreteValuation))
                else:
                    clazz = dynamic_class("RationalFunctionFieldMappedValuation_infinite", (clazz, InfiniteDiscretePseudoValuation))
                return parent.__make_element_class__(clazz)(parent, valuation, to_valuation_domain, from_valuation_domain)
            return parent.__make_element_class__(FunctionFieldExtensionMappedValuation)(parent, valuation, to_valuation_domain, from_valuation_domain)

        if domain is valuation.domain():
            # we cannot just return valuation in this case
            # as this would break uniqueness and pickling
            raise ValueError("valuation must not be a valuation on domain yet but %r is a valuation on %r" % (valuation, domain))

        if domain.base_field() is domain:
            # valuation is a base valuation on K[x] that induces a valuation on K(x)
            if valuation.restriction(domain.constant_base_field()).is_trivial() and valuation.is_discrete_valuation():
                # valuation corresponds to a finite place
                return parent.__make_element_class__(FiniteRationalFunctionFieldValuation)(parent, valuation)
            else:
                from sage.structure.dynamic_class import dynamic_class
                clazz = NonClassicalRationalFunctionFieldValuation
                if valuation.is_discrete_valuation():
                    clazz = dynamic_class("NonClassicalRationalFunctionFieldValuation_discrete", (clazz, DiscreteFunctionFieldValuation_base))
                else:
                    clazz = dynamic_class("NonClassicalRationalFunctionFieldValuation_negative_infinite", (clazz, NegativeInfiniteDiscretePseudoValuation))
                return parent.__make_element_class__(clazz)(parent, valuation)
        else:
            # valuation is a limit valuation that singles out an extension
            return parent.__make_element_class__(FunctionFieldFromLimitValuation)(parent, valuation, domain.polynomial(), extra_args['approximants'])

        raise NotImplementedError("valuation on %r from %r on %r" % (domain, valuation, valuation.domain()))

FunctionFieldValuation = FunctionFieldValuationFactory("sage.rings.function_field.valuation.FunctionFieldValuation")


class FunctionFieldValuation_base(DiscretePseudoValuation):
    r"""
    Abstract base class for any discrete (pseudo-)valuation on a function
    field.

    TESTS::

        sage: K.<x> = FunctionField(QQ)
        sage: v = K.valuation(x)  # indirect doctest
        sage: from sage.rings.function_field.valuation import FunctionFieldValuation_base
        sage: isinstance(v, FunctionFieldValuation_base)
        True

    """


class DiscreteFunctionFieldValuation_base(DiscreteValuation):
    r"""
    Base class for discrete valuations on function fields.

    TESTS::

        sage: K.<x> = FunctionField(QQ)
        sage: v = K.valuation(x)  # indirect doctest
        sage: from sage.rings.function_field.valuation import DiscreteFunctionFieldValuation_base
        sage: isinstance(v, DiscreteFunctionFieldValuation_base)
        True

    """
    def extensions(self, L):
        r"""
        Return the extensions of this valuation to ``L``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)                                          # needs sage.rings.function_field
            sage: v.extensions(L)                                                       # needs sage.rings.function_field
            [(x)-adic valuation]

        TESTS:

        Valuations over the infinite place::

            sage: v = K.valuation(1/x)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - 1/(x^2 + 1))                                # needs sage.rings.function_field
            sage: sorted(v.extensions(L), key=str)                                      # needs sage.rings.function_field
            [[ Valuation at the infinite place, v(y + 1/x) = 3 ]-adic valuation,
             [ Valuation at the infinite place, v(y - 1/x) = 3 ]-adic valuation]

        Iterated extensions over the infinite place::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)
            sage: v = K.valuation(1/x)
            sage: w = v.extension(L)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: w.extension(M)                # not implemented
            Traceback (most recent call last):
            ...
            NotImplementedError

        A case that caused some trouble at some point::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(v)

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - x^4 - 1)                                    # needs sage.rings.function_field
            sage: v.extensions(L)                                                       # needs sage.rings.function_field
            [2-adic valuation]

        Test that this works in towers::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y - x)
            sage: R.<z> = L[]
            sage: L.<z> = L.extension(z - y)
            sage: v = K.valuation(x)
            sage: v.extensions(L)
            [(x)-adic valuation]
        """
        K = self.domain()
        from sage.categories.function_fields import FunctionFields
        if L is K:
            return [self]
        if L in FunctionFields():
            if K.is_subring(L):
                if L.base() is K:
                    # L = K[y]/(G) is a simple extension of the domain of this valuation
                    G = L.polynomial()
                    if not G.is_monic():
                        G = G / G.leading_coefficient()
                    if any(self(c) < 0 for c in G.coefficients()):
                        # rewrite L = K[u]/(H) with H integral and compute the extensions
                        from sage.rings.valuation.gauss_valuation import GaussValuation
                        g = GaussValuation(G.parent(), self)
                        y_to_u, u_to_y, H = g.monic_integral_model(G)
                        M = K.extension(H, names=L.variable_names())
                        H_extensions = self.extensions(M)

                        from sage.rings.morphism import RingHomomorphism_im_gens
                        if isinstance(y_to_u, RingHomomorphism_im_gens) and isinstance(u_to_y, RingHomomorphism_im_gens):
                            return [L.valuation((w, L.hom([M(y_to_u(y_to_u.domain().gen()))]), M.hom([L(u_to_y(u_to_y.domain().gen()))]))) for w in H_extensions]
                        raise NotImplementedError
                    return [L.valuation(w) for w in self.mac_lane_approximants(L.polynomial(), require_incomparability=True)]
                elif L.base() is not L and K.is_subring(L):
                    # recursively call this method for the tower of fields
                    from operator import add
                    from functools import reduce
                    A = [base_valuation.extensions(L) for base_valuation in self.extensions(L.base())]
                    return reduce(add, A, [])
                elif L.constant_base_field() is not K.constant_base_field() and K.constant_base_field().is_subring(L):
                    # subclasses should override this method and handle this case, so we never get here
                    raise NotImplementedError("Cannot compute the extensions of %r from %r to %r since the base ring changes." % (self, self.domain(), L))
        raise NotImplementedError("extension of %r from %r to %r not implemented" % (self, K, L))


class RationalFunctionFieldValuation_base(FunctionFieldValuation_base):
    r"""
    Base class for valuations on rational function fields.

    TESTS::

        sage: K.<x> = FunctionField(GF(2))
        sage: v = K.valuation(x)  # indirect doctest
        sage: from sage.rings.function_field.valuation import RationalFunctionFieldValuation_base
        sage: isinstance(v, RationalFunctionFieldValuation_base)
        True

    """
    @cached_method
    def element_with_valuation(self, s):
        r"""
        Return an element with valuation ``s``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 + 6)
            sage: v = K.valuation(2)
            sage: R.<x> = K[]
            sage: w = GaussValuation(R, v).augmentation(x, 1/123)
            sage: K.<x> = FunctionField(K)
            sage: w = w.extension(K)
            sage: w.element_with_valuation(122/123)
            2/x
            sage: w.element_with_valuation(1)
            2

        """
        constant_valuation = self.restriction(self.domain().constant_base_field())
        if constant_valuation.is_trivial():
            return super().element_with_valuation(s)

        a, b = self.value_group()._element_with_valuation(constant_valuation.value_group(), s)
        ret = self.uniformizer()**a * constant_valuation.element_with_valuation(constant_valuation.value_group().gen()*b)

        return self.simplify(ret, error=s)


class ClassicalFunctionFieldValuation_base(DiscreteFunctionFieldValuation_base):
    r"""
    Base class for discrete valuations on rational function fields that come
    from points on the projective line.

    TESTS::

        sage: K.<x> = FunctionField(GF(5))
        sage: v = K.valuation(x)  # indirect doctest
        sage: from sage.rings.function_field.valuation import ClassicalFunctionFieldValuation_base
        sage: isinstance(v, ClassicalFunctionFieldValuation_base)
        True

    """
    def _test_classical_residue_field(self, **options):
        r"""
        Check correctness of the residue field of a discrete valuation at a
        classical point.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x^2 + 1)
            sage: v._test_classical_residue_field()

        """
        tester = self._tester(**options)

        tester.assertTrue(self.domain().constant_base_field().is_subring(self.residue_field()))

    def _ge_(self, other):
        r"""
        Return whether ``self`` is greater or equal to ``other`` everywhere.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x^2 + 1)
            sage: w = K.valuation(x)
            sage: v >= w
            False
            sage: w >= v
            False

        """
        if other.is_trivial():
            return other.is_discrete_valuation()
        if isinstance(other, ClassicalFunctionFieldValuation_base):
            return self == other
        super()._ge_(other)


class InducedRationalFunctionFieldValuation_base(FunctionFieldValuation_base):
    r"""
    Base class for function field valuation induced by a valuation on the
    underlying polynomial ring.

    TESTS::

        sage: K.<x> = FunctionField(QQ)
        sage: v = K.valuation(x^2 + 1) # indirect doctest

    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x) # indirect doctest
            sage: from sage.rings.function_field.valuation import InducedRationalFunctionFieldValuation_base
            sage: isinstance(v, InducedRationalFunctionFieldValuation_base)
            True

        """
        FunctionFieldValuation_base.__init__(self, parent)

        domain = parent.domain()
        if base_valuation.domain() is not domain._ring:
            raise ValueError("base valuation must be defined on %r but %r is defined on %r" % (domain._ring, base_valuation, base_valuation.domain()))

        self._base_valuation = base_valuation

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.valuation(x).uniformizer()
            x

        """
        return self.domain()(self._base_valuation.uniformizer())

    def lift(self, F):
        r"""
        Return a lift of ``F`` to the domain of this valuation such
        that :meth:`reduce` returns the original element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x)
            sage: v.lift(0)
            0
            sage: v.lift(1)
            1

        """
        F = self.residue_ring().coerce(F)
        if F in self._base_valuation.residue_ring():
            num = self._base_valuation.residue_ring()(F)
            den = self._base_valuation.residue_ring()(1)
        elif F in self._base_valuation.residue_ring().fraction_field():
            num = self._base_valuation.residue_ring()(F.numerator())
            den = self._base_valuation.residue_ring()(F.denominator())
        else:
            raise NotImplementedError("lifting not implemented for this valuation")

        return self.domain()(self._base_valuation.lift(num)) / self.domain()(self._base_valuation.lift(den))

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.valuation(x).value_group()
            Additive Abelian Group generated by 1

        """
        return self._base_valuation.value_group()

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` in :meth:`residue_ring`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x^2 + 1)
            sage: v.reduce(x)
            u1

        """
        f = self.domain().coerce(f)

        if self(f) > 0:
            return self.residue_field().zero()
        if self(f) < 0:
            raise ValueError("cannot reduce element of negative valuation")

        base = self._base_valuation

        num = f.numerator()
        den = f.denominator()

        assert base(num) == base(den)
        shift = base.element_with_valuation(-base(num))
        num *= shift
        den *= shift
        ret = base.reduce(num) / base.reduce(den)
        assert not ret.is_zero()
        return self.residue_field()(ret)

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.valuation(x^2 + 1) # indirect doctest
            (x^2 + 1)-adic valuation

        """
        from sage.rings.valuation.augmented_valuation import AugmentedValuation_base
        from sage.rings.valuation.gauss_valuation import GaussValuation
        if isinstance(self._base_valuation, AugmentedValuation_base):
            if self._base_valuation._base_valuation == GaussValuation(self.domain()._ring, TrivialValuation(self.domain().constant_base_field())):
                if self._base_valuation._mu == 1:
                    return "(%r)-adic valuation" % (self._base_valuation.phi())
        vK = self._base_valuation.restriction(self._base_valuation.domain().base_ring())
        if self._base_valuation == GaussValuation(self.domain()._ring, vK):
            return repr(vK)
        return "Valuation on rational function field induced by %s" % self._base_valuation

    def extensions(self, L):
        r"""
        Return all extensions of this valuation to ``L`` which has a larger
        constant field than the domain of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x^2 + 1)
            sage: L.<x> = FunctionField(GaussianIntegers().fraction_field())
            sage: v.extensions(L)  # indirect doctest
            [(x - I)-adic valuation, (x + I)-adic valuation]

        """
        K = self.domain()
        if L is K:
            return [self]

        from sage.categories.function_fields import FunctionFields
        if (L in FunctionFields()
            and K.is_subring(L)
            and L.base() is L
            and L.constant_base_field() is not K.constant_base_field()
            and K.constant_base_field().is_subring(L.constant_base_field())):
            # The above condition checks whether L is an extension of K that
            # comes from an extension of the field of constants
            # Condition "L.base() is L" is important so we do not call this
            # code for extensions from K(x) to K(x)(y)

            # We extend the underlying valuation on the polynomial ring
            W = self._base_valuation.extensions(L._ring)
            return [L.valuation(w) for w in W]

        return super().extensions(L)

    def _call_(self, f):
        r"""
        Evaluate this valuation at the function ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x)  # indirect doctest
            sage: v((x+1)/x^2)
            -2

        """
        return self._base_valuation(f.numerator()) - self._base_valuation(f.denominator())

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.valuation(x).residue_ring()
            Rational Field

        """
        return self._base_valuation.residue_ring().fraction_field()

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.valuation(x).restriction(QQ)
            Trivial valuation on Rational Field

        """
        if ring.is_subring(self._base_valuation.domain()):
            return self._base_valuation.restriction(ring)
        return super().restriction(ring)

    def simplify(self, f, error=None, force=False):
        r"""
        Return a simplified version of ``f``.

        Produce an element which differs from ``f`` by an element of
        valuation strictly greater than the valuation of ``f`` (or strictly
        greater than ``error`` if set.)

        If ``force`` is not set, then expensive simplifications may be avoided.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(2)
            sage: f = (x + 1)/(x - 1)

        As the coefficients of this fraction are small, we do not simplify as
        this could be very costly in some cases::

            sage: v.simplify(f)
            (x + 1)/(x - 1)

        However, simplification can be forced::

            sage: v.simplify(f, force=True)
            3

        """
        f = self.domain().coerce(f)

        if error is None:
            # if the caller was sure that we should simplify, then we should try to do the best simplification possible
            error = self(f) if force else self.upper_bound(f)

        from sage.rings.infinity import infinity
        if error is infinity:
            return f

        numerator = f.numerator()
        denominator = f.denominator()

        v_numerator = self._base_valuation(numerator)
        v_denominator = self._base_valuation(denominator)

        if v_numerator - v_denominator > error:
            return self.domain().zero()

        if error == -infinity:
            # This case is not implemented yet, so we just return f which is always safe.
            return f

        numerator = self.domain()(self._base_valuation.simplify(numerator, error=error+v_denominator, force=force))
        denominator = self.domain()(self._base_valuation.simplify(denominator, error=max(v_denominator, error - v_numerator + 2*v_denominator), force=force))

        ret = numerator/denominator
        assert self(ret - f) > error
        return ret

    def _relative_size(self, f):
        r"""
        Return an estimate on the coefficient size of ``f``.

        The number returned is an estimate on the factor between the number of
        bits used by ``f`` and the minimal number of bits used by an element
        congruent to ``f``.

        This can be used by :meth:`simplify` to decide whether simplification
        of coefficients is going to lead to a significant shrinking of the
        coefficients of ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(0)
            sage: f = (x + 1024)/(x - 1024)

        Here we report a small size, as the numerator and the denominator
        independently cannot be simplified much::

            sage: v._relative_size(f)
            1

        However, a forced simplification, finds that we could have saved many
        more bits::

            sage: v.simplify(f, force=True)
            -1

        """
        return max(self._base_valuation._relative_size(f.numerator()), self._base_valuation._relative_size(f.denominator()))


class FiniteRationalFunctionFieldValuation(InducedRationalFunctionFieldValuation_base, ClassicalFunctionFieldValuation_base, RationalFunctionFieldValuation_base):
    r"""
    Valuation of a finite place of a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: v = K.valuation(x + 1); v  # indirect doctest
        (x + 1)-adic valuation

    A finite place with residual degree::

        sage: w = K.valuation(x^2 + 1); w
        (x^2 + 1)-adic valuation

    A finite place with ramification::

        sage: K.<t> = FunctionField(GF(3))
        sage: L.<x> = FunctionField(K)
        sage: u = L.valuation(x^3 - t); u
        (x^3 + 2*t)-adic valuation

    A finite place with residual degree and ramification::

        sage: q = L.valuation(x^6 - t); q
        (x^6 + 2*t)-adic valuation

    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(x + 1)
            sage: from sage.rings.function_field.valuation import FiniteRationalFunctionFieldValuation
            sage: isinstance(v, FiniteRationalFunctionFieldValuation)
            True

        """
        InducedRationalFunctionFieldValuation_base.__init__(self, parent, base_valuation)
        ClassicalFunctionFieldValuation_base.__init__(self, parent)
        RationalFunctionFieldValuation_base.__init__(self, parent)


class NonClassicalRationalFunctionFieldValuation(InducedRationalFunctionFieldValuation_base, RationalFunctionFieldValuation_base):
    r"""
    Valuation induced by a valuation on the underlying polynomial ring which is
    non-classical.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: v = GaussValuation(QQ['x'], QQ.valuation(2))
        sage: w = K.valuation(v); w  # indirect doctest
        2-adic valuation

    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS:

        There is some support for discrete pseudo-valuations on rational
        function fields in the code. However, since these valuations must send
        elements to `-\infty`, they are not supported yet::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(QQ['x'], QQ.valuation(2)).augmentation(x, infinity)
            sage: K.<x> = FunctionField(QQ)
            sage: w = K.valuation(v)
            sage: from sage.rings.function_field.valuation import NonClassicalRationalFunctionFieldValuation
            sage: isinstance(w, NonClassicalRationalFunctionFieldValuation)
            True

        """
        InducedRationalFunctionFieldValuation_base.__init__(self, parent, base_valuation)
        RationalFunctionFieldValuation_base.__init__(self, parent)

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: v = valuations.GaussValuation(QQ['x'], QQ.valuation(2))
            sage: w = K.valuation(v)
            sage: w.residue_ring()
            Rational function field in x over Finite Field of size 2

            sage: R.<x> = QQ[]
            sage: vv = v.augmentation(x, 1)
            sage: w = K.valuation(vv)
            sage: w.residue_ring()
            Rational function field in x over Finite Field of size 2

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + 2*x)                                        # needs sage.rings.function_field
            sage: w.extension(L).residue_ring()                                         # needs sage.rings.function_field
            Function field in u2 defined by u2^2 + x

        TESTS:

        This still works for pseudo-valuations::

            sage: R.<x> = QQ[]
            sage: v = valuations.GaussValuation(R, QQ.valuation(2))
            sage: vv = v.augmentation(x, infinity)
            sage: K.<x> = FunctionField(QQ)
            sage: w = K.valuation(vv)
            sage: w.residue_ring()
            Finite Field of size 2

        """
        if not self.is_discrete_valuation():
            # A pseudo valuation attaining negative infinity does typically not have a function field as its residue ring
            return super().residue_ring()
        return self._base_valuation.residue_ring().fraction_field().function_field()


class FunctionFieldFromLimitValuation(FiniteExtensionFromLimitValuation, DiscreteFunctionFieldValuation_base):
    r"""
    A valuation on a finite extensions of function fields `L=K[y]/(G)` where `K` is
    another function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - (x^2 + x + 1))                                  # needs sage.rings.function_field
        sage: v = K.valuation(x - 1) # indirect doctest                                 # needs sage.rings.function_field
        sage: w = v.extension(L); w                                                     # needs sage.rings.function_field
        (x - 1)-adic valuation

    """
    def __init__(self, parent, approximant, G, approximants):
        r"""
        TESTS::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x^2 + x + 1))
            sage: v = K.valuation(x - 1) # indirect doctest
            sage: w = v.extension(L)
            sage: from sage.rings.function_field.valuation import FunctionFieldFromLimitValuation
            sage: isinstance(w, FunctionFieldFromLimitValuation)
            True

        """
        FiniteExtensionFromLimitValuation.__init__(self, parent, approximant, G, approximants)
        DiscreteFunctionFieldValuation_base.__init__(self, parent)

    def _to_base_domain(self, f):
        r"""
        Return ``f`` as an element of the domain of the underlying limit valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x^2 + x + 1))
            sage: v = K.valuation(x - 1) # indirect doctest
            sage: w = v.extension(L)
            sage: w._to_base_domain(y).parent()
            Univariate Polynomial Ring in y over Rational function field in x over Rational Field

        """
        return f.element()

    def scale(self, scalar):
        r"""
        Return this valuation scaled by ``scalar``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x^2 + x + 1))
            sage: v = K.valuation(x - 1) # indirect doctest
            sage: w = v.extension(L)
            sage: 3*w
            3 * (x - 1)-adic valuation

        """
        if scalar in QQ and scalar > 0 and scalar != 1:
            return self.domain().valuation(self._base_valuation._initial_approximation.scale(scalar))
        return super().scale(scalar)


class FunctionFieldMappedValuation_base(FunctionFieldValuation_base, MappedValuation_base):
    r"""
    A valuation on a function field which relies on a ``base_valuation`` on an
    isomorphic function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2))
        sage: v = K.valuation(1/x); v
        Valuation at the infinite place

    """
    def __init__(self, parent, base_valuation, to_base_valuation_domain, from_base_valuation_domain):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: v = K.valuation(1/x)
            sage: from sage.rings.function_field.valuation import FunctionFieldMappedValuation_base
            sage: isinstance(v, FunctionFieldMappedValuation_base)
            True

        """
        FunctionFieldValuation_base.__init__(self, parent)
        MappedValuation_base.__init__(self, parent, base_valuation)

        self._to_base = to_base_valuation_domain
        self._from_base = from_base_valuation_domain

    def _to_base_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of ``_base_valuation``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)                                    # needs sage.rings.function_field
            sage: v = K.valuation(1/x)
            sage: w = v.extension(L)                                                    # needs sage.rings.function_field
            sage: w._to_base_domain(y)                                                  # needs sage.rings.function_field
            x^2*y

        """
        return self._to_base(f)

    def _from_base_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)                                    # needs sage.rings.function_field
            sage: v = K.valuation(1/x)
            sage: w = v.extension(L)                                                    # needs sage.rings.function_field
            sage: w._from_base_domain(w._to_base_domain(y))                             # needs sage.rings.function_field
            y

        r"""
        return self._from_base(f)

    def scale(self, scalar):
        r"""
        Return this valuation scaled by ``scalar``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)                                    # needs sage.rings.function_field
            sage: v = K.valuation(1/x)
            sage: w = v.extension(L)                                                    # needs sage.rings.function_field
            sage: 3*w                                                                   # needs sage.rings.function_field
            3 * (x)-adic valuation (in Rational function field in x over Finite Field of size 2 after x |--> 1/x)

        """
        from sage.rings.rational_field import QQ
        if scalar in QQ and scalar > 0 and scalar != 1:
            return self.domain().valuation((self._base_valuation.scale(scalar), self._to_base, self._from_base))
        return super().scale(scalar)

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)                                    # needs sage.rings.function_field
            sage: v = K.valuation(1/x)
            sage: v.extension(L)  # indirect doctest                                    # needs sage.rings.function_field
            Valuation at the infinite place

        """
        to_base = repr(self._to_base)
        if hasattr(self._to_base, '_repr_defn'):
            to_base = self._to_base._repr_defn().replace('\n', ', ')
        return "%r (in %r after %s)" % (self._base_valuation, self._base_valuation.domain(), to_base)

    def is_discrete_valuation(self):
        r"""
        Return whether this is a discrete valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^4 - 1)
            sage: v = K.valuation(1/x)
            sage: w0,w1 = v.extensions(L)
            sage: w0.is_discrete_valuation()
            True

        """
        return self._base_valuation.is_discrete_valuation()


class FunctionFieldMappedValuationRelative_base(FunctionFieldMappedValuation_base):
    r"""
    A valuation on a function field which relies on a ``base_valuation`` on an
    isomorphic function field and which is such that the map from and to the
    other function field is the identity on the constant field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2))
        sage: v = K.valuation(1/x); v
        Valuation at the infinite place

    """
    def __init__(self, parent, base_valuation, to_base_valuation_domain, from_base_valuation_domain):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: v = K.valuation(1/x)
            sage: from sage.rings.function_field.valuation import FunctionFieldMappedValuationRelative_base
            sage: isinstance(v, FunctionFieldMappedValuationRelative_base)
            True

        """
        FunctionFieldMappedValuation_base.__init__(self, parent, base_valuation, to_base_valuation_domain, from_base_valuation_domain)
        if self.domain().constant_base_field() is not base_valuation.domain().constant_base_field():
            raise ValueError("constant fields must be identical but they differ for %r and %r" % (self.domain(), base_valuation.domain()))

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: K.valuation(1/x).restriction(GF(2))
            Trivial valuation on Finite Field of size 2

        """
        if ring.is_subring(self.domain().constant_base_field()):
            return self._base_valuation.restriction(ring)
        return super().restriction(ring)


class RationalFunctionFieldMappedValuation(FunctionFieldMappedValuationRelative_base, RationalFunctionFieldValuation_base):
    r"""
    Valuation on a rational function field that is implemented after a map to
    an isomorphic rational function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<x> = QQ[]
        sage: w = GaussValuation(R, QQ.valuation(2)).augmentation(x, 1)
        sage: w = K.valuation(w)
        sage: v = K.valuation((w, K.hom([~K.gen()]), K.hom([~K.gen()]))); v
        Valuation on rational function field induced by
        [ Gauss valuation induced by 2-adic valuation, v(x) = 1 ]
        (in Rational function field in x over Rational Field after x |--> 1/x)

    """
    def __init__(self, parent, base_valuation, to_base_valuation_doain, from_base_valuation_domain):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, QQ.valuation(2)).augmentation(x, 1)
            sage: w = K.valuation(w)
            sage: v = K.valuation((w, K.hom([~K.gen()]), K.hom([~K.gen()])))
            sage: from sage.rings.function_field.valuation import RationalFunctionFieldMappedValuation
            sage: isinstance(v, RationalFunctionFieldMappedValuation)
            True

        """
        FunctionFieldMappedValuationRelative_base.__init__(self, parent, base_valuation, to_base_valuation_doain, from_base_valuation_domain)
        RationalFunctionFieldValuation_base.__init__(self, parent)


class InfiniteRationalFunctionFieldValuation(FunctionFieldMappedValuationRelative_base, RationalFunctionFieldValuation_base, ClassicalFunctionFieldValuation_base):
    r"""
    Valuation of the infinite place of a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: v = K.valuation(1/x)  # indirect doctest

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(1/x)  # indirect doctest
            sage: from sage.rings.function_field.valuation import InfiniteRationalFunctionFieldValuation
            sage: isinstance(v, InfiniteRationalFunctionFieldValuation)
            True

        """
        x = parent.domain().gen()
        FunctionFieldMappedValuationRelative_base.__init__(self, parent, FunctionFieldValuation(parent.domain(), x), parent.domain().hom([1/x]), parent.domain().hom([1/x]))
        RationalFunctionFieldValuation_base.__init__(self, parent)
        ClassicalFunctionFieldValuation_base.__init__(self, parent)

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.valuation(1/x)  # indirect doctest
            Valuation at the infinite place

        """
        return "Valuation at the infinite place"


class FunctionFieldExtensionMappedValuation(FunctionFieldMappedValuationRelative_base):
    r"""
    A valuation on a finite extensions of function fields `L=K[y]/(G)` where `K` is
    another function field which redirects to another ``base_valuation`` on an
    isomorphism function field `M=K[y]/(H)`.

    The isomorphisms must be trivial on ``K``.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2))
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 + y + x^3)                                        # needs sage.rings.function_field
        sage: v = K.valuation(1/x)
        sage: w = v.extension(L)                                                        # needs sage.rings.function_field

        sage: w(x)                                                                      # needs sage.rings.function_field
        -1
        sage: w(y)                                                                      # needs sage.rings.function_field
        -3/2
        sage: w.uniformizer()                                                           # needs sage.rings.function_field
        1/x^2*y

    TESTS::

        sage: from sage.rings.function_field.valuation import FunctionFieldExtensionMappedValuation
        sage: isinstance(w, FunctionFieldExtensionMappedValuation)                      # needs sage.rings.function_field
        True

    """
    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)
            sage: v = K.valuation(1/x)
            sage: w = v.extension(L); w
            Valuation at the infinite place

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - 1/x^2 - 1)
            sage: v = K.valuation(1/x)
            sage: w = v.extensions(L); w
            [[ Valuation at the infinite place, v(y + 1) = 2 ]-adic valuation,
             [ Valuation at the infinite place, v(y - 1) = 2 ]-adic valuation]

        """
        assert(self.domain().base() is not self.domain())
        if repr(self._base_valuation) == repr(self.restriction(self.domain().base())):
            return repr(self._base_valuation)
        return super()._repr_()

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + y + x^3)
            sage: v = K.valuation(1/x)
            sage: w = v.extension(L)
            sage: w.restriction(K) is v
            True
        """
        if ring.is_subring(self.domain().base()):
            return self._base_valuation.restriction(ring)
        return super().restriction(ring)
