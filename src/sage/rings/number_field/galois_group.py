"""
Galois Groups of Number Fields

AUTHORS:

- William Stein (2004, 2005): initial version
- David Loeffler (2009): rewrite to give explicit homomorphism groups

TESTS:

Standard test of pickleability::

    sage: G = NumberField(x^3 + 2, 'alpha').galois_group(type="pari"); G
    Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in alpha with defining polynomial x^3 + 2
    sage: G == loads(dumps(G))
    True

    sage: G = NumberField(x^3 + 2, 'alpha').galois_group(names='beta'); G
    Galois group of Galois closure in beta of Number Field in alpha with defining polynomial x^3 + 2
    sage: G == loads(dumps(G))
    True
"""

from sage.structure.sage_object import SageObject
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.misc.cachefunc import cached_method
from sage.libs.pari.all import pari
from sage.rings.infinity import infinity
from sage.rings.number_field.number_field import refine_embedding
from sage.rings.number_field.morphism import NumberFieldHomomorphism_im_gens

class GaloisGroup_v1(SageObject):
    r"""
    A wrapper around a class representing an abstract transitive group.

    This is just a fairly minimal object at present.  To get the underlying
    group, do ``G.group()``, and to get the corresponding number field do
    ``G.number_field()``. For a more sophisticated interface use the
    ``type=None`` option.

    EXAMPLES::

        sage: K = QQ[2^(1/3)]
        sage: G = K.galois_group(type="pari"); G
        Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in a with defining polynomial x^3 - 2
        sage: G.order()
        6
        sage: G.group()
        PARI group [6, -1, 2, "S3"] of degree 3
        sage: G.number_field()
        Number Field in a with defining polynomial x^3 - 2
    """

    def __init__(self, group, number_field):
        """
        Create a Galois group.

        EXAMPLES::

            sage: NumberField([x^2 + 1, x^2 + 2],'a').galois_group(type="pari")
            Galois group PARI group [4, 1, 2, "E(4) = 2[x]2"] of degree 4 of the Number Field in a0 with defining polynomial x^2 + 1 over its base field
        """
        self.__group = group
        self.__number_field = number_field

    def __cmp__(self, other):
        """
        Compare two number field Galois groups.  First the number
        fields are compared, then the Galois groups if the number
        fields are equal.  (Of course, if the number fields are the
        same, the Galois groups are automatically equal.)

        EXAMPLES::

            sage: G = NumberField(x^3 + 2, 'alpha').galois_group(type="pari")
            sage: H = QQ[sqrt(2)].galois_group(type="pari")
            sage: cmp(G,H)
            -1
            sage: H == H
            True
            sage: G == G
            True
        """
        if not isinstance(other, GaloisGroup_v1):
            return cmp(type(self), type(other))
        return cmp( (self.__number_field, self.__group),
                    (other.__number_field, other.__group) )

    def __repr__(self):
        """
        Display print representation of a Galois group.

        EXAMPLES::

            sage: G = NumberField(x^4 + 2*x + 2, 'a').galois_group(type="pari")
            sage: G.__repr__()
            'Galois group PARI group [24, -1, 5, "S4"] of degree 4 of the Number Field in a with defining polynomial x^4 + 2*x + 2'
        """
        return "Galois group %s of the %s"%(
            self.__group, self.__number_field)

    def group(self):
        """
        Return the underlying abstract group.

        EXAMPLES::

            sage: G = NumberField(x^3 + 2*x + 2, 'theta').galois_group(type="pari")
            sage: H = G.group(); H
            PARI group [6, -1, 2, "S3"] of degree 3
            sage: P = H.permutation_group(); P  # optional - database_gap
            Transitive group number 2 of degree 3
            sage: list(P)                       # optional
            [(), (2,3), (1,2), (1,2,3), (1,3,2), (1,3)]
        """
        return self.__group

    def order(self):
        """
        Return the order of this Galois group.

        EXAMPLES::

            sage: G = NumberField(x^5 + 2, 'theta_1').galois_group(type="pari"); G
            Galois group PARI group [20, -1, 3, "F(5) = 5:4"] of degree 5 of the Number Field in theta_1 with defining polynomial x^5 + 2
            sage: G.order()
            20
        """
        return self.__group.order()

    def number_field(self):
        """
        Return the number field of which this is the Galois group.

        EXAMPLES::

            sage: G = NumberField(x^6 + 2, 't').galois_group(type="pari"); G
            Galois group PARI group [12, -1, 3, "D(6) = S(3)[x]2"] of degree 6 of the Number Field in t with defining polynomial x^6 + 2
            sage: G.number_field()
            Number Field in t with defining polynomial x^6 + 2
        """
        return self.__number_field


class GaloisGroup_v2(PermutationGroup_generic):

    r"""
    The Galois group of an (absolute) number field.

    .. note::

        We define the Galois group of a non-normal field K to be the
        Galois group of its Galois closure L, and elements are stored as
        permutations of the roots of the defining polynomial of L, *not* as
        permutations of the roots (in L) of the defining polynomial of K. The
        latter would probably be preferable, but is harder to implement. Thus
        the permutation group that is returned is always simply-transitive.

        The 'arithmetical' features (decomposition and ramification groups,
        Artin symbols etc) are only available for Galois fields.
    """

    def __init__(self, number_field, names=None):
        r"""
        Create a Galois group.

        EXAMPLES::

            sage: QuadraticField(-23,'a').galois_group()
            Galois group of Number Field in a with defining polynomial x^2 + 23
            sage: NumberField(x^3 - 2, 'b').galois_group()
            Traceback (most recent call last):
            ...
            TypeError: You must specify the name of the generator.
            sage: NumberField(x^3 - 2, 'b').galois_group(names="c")
            Galois group of Galois closure in c of Number Field in b with defining polynomial x^3 - 2
        """

        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        self._number_field = number_field

        if not number_field.is_galois():
            self._galois_closure, self._gc_map = number_field.galois_closure(names=names, map=True)
        else:
            self._galois_closure, self._gc_map = (number_field, number_field.hom(number_field.gen(), number_field))

        self._pari_gc = self._galois_closure._pari_()

        g = self._pari_gc.galoisinit()
        self._pari_data = g

        PermutationGroup_generic.__init__(self, sorted(g[6]))

        # PARI computes all the elements of self anyway, so we might as well store them
        self._elts = sorted([self(x, check=False) for x in g[5]])

    def _element_class(self):
        r"""
        Return the class to be used for creating elements of this group, which
        is GaloisGroupElement.

        EXAMPLE::

            sage: F.<z> = CyclotomicField(7)
            sage: G = F.galois_group()
            sage: G._element_class()
            <class 'sage.rings.number_field.galois_group.GaloisGroupElement'>

        We test that a method inherited from PermutationGroup_generic returns
        the right type of element (see trac #133)::

            sage: phi = G.random_element()
            sage: type(phi)
            <class 'sage.rings.number_field.galois_group.GaloisGroupElement'>
            sage: phi(z) # random
            z^3
        """
        return GaloisGroupElement

    def __call__(self, x, check=True):
        r""" Create an element of self from x. Here x had better be one of:
        -- the integer 1, denoting the identity of G
        -- an element of G
        -- a permutation of the right length which defines an element of G, or anything that
            coerces into a permutation of the right length
        -- an abstract automorphism of the underlying number field.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-23)
            sage: G = K.galois_group()
            sage: G(1)
            ()
            sage: G(G.gens()[0])
            (1,2)
            sage: G([(1,2)])
            (1,2)
            sage: G(K.hom(-a, K))
            (1,2)
        """
        if x == 1:
            return self.identity()

        from sage.rings.number_field.morphism import NumberFieldHomomorphism_im_gens
        if isinstance(x, NumberFieldHomomorphism_im_gens) and x.parent() == self.number_field().Hom(self.number_field()):
            l = [g for g in self if g.as_hom() == x]
            if len(l) != 1: raise ArithmeticError
            return l[0]
        return GaloisGroupElement(x, parent=self, check=check)

    def is_galois(self):
        r"""
        Return True if the underlying number field of self is actually Galois.

        EXAMPLE::

            sage: NumberField(x^3 - x + 1,'a').galois_group(names='b').is_galois()
            False
            sage: NumberField(x^2 - x + 1,'a').galois_group().is_galois()
            True
        """
        if self._number_field == self._galois_closure:
            return True
        else:
            return False

    def ngens(self):
        r""" Number of generators of self.

        EXAMPLE::

            sage: QuadraticField(-23, 'a').galois_group().ngens()
            1
        """
        return len(self._gens)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: G = QuadraticField(-23, 'a').galois_group()
            sage: G._repr_()
            'Galois group of Number Field in a with defining polynomial x^2 + 23'
            sage: G = NumberField(x^3 - 2, 'a').galois_group(names='b')
            sage: G._repr_()
            'Galois group of Galois closure in b of Number Field in a with defining polynomial x^3 - 2'
        """
        if self.is_galois():
            return "Galois group of %s" % self.number_field()
        else:
            return "Galois group of Galois closure in %s of %s" % (self.splitting_field().gen(), self.number_field())

    def number_field(self):
        r"""
        The ambient number field.

        EXAMPLE::

            sage: K = NumberField(x^3 - x + 1, 'a')
            sage: K.galois_group(names='b').number_field() is K
            True
        """
        return self._number_field

    def splitting_field(self):
        r"""
        The Galois closure of the ambient number field.

        EXAMPLE::

            sage: K = NumberField(x^3 - x + 1, 'a')
            sage: K.galois_group(names='b').splitting_field()
            Number Field in b with defining polynomial x^6 - 14*x^4 - 20*x^3 + 49*x^2 + 140*x + 307
            sage: L = QuadraticField(-23, 'c'); L.galois_group().splitting_field() is L
            True
        """
        return self._galois_closure

    def list(self):
        r"""
        List of the elements of self.

        EXAMPLE::

            sage: NumberField(x^3 - 3*x + 1,'a').galois_group().list()
            [(), (1,2,3), (1,3,2)]
        """
        return self._elts

    def subgroup(self, elts):
        r"""
        Return the subgroup of self with the given elements. Mostly for internal use.

        EXAMPLE::

            sage: G = NumberField(x^3 - x - 1,'a').galois_closure('b').galois_group()
            sage: G.subgroup([ G(1), G([(1,5,2),(3,4,6)]), G([(1,2,5),(3,6,4)])])
            Subgroup [(), (1,5,2)(3,4,6), (1,2,5)(3,6,4)] of Galois group of Number Field in b with defining polynomial x^6 - 14*x^4 + 20*x^3 + 49*x^2 - 140*x + 307
        """
        if len(elts) == self.order():
            return self
        else:
            return GaloisGroup_subgroup(self, elts)

    # Proper number theory starts here. All the functions below make no sense
    # unless the field is Galois.

    def decomposition_group(self, P):
        """
        Decomposition group of a prime ideal P, i.e. the subgroup of elements
        that map P to itself. This is the same as the Galois group of the
        extension of local fields obtained by completing at P.

        This function will raise an error if P is not prime or the given number
        field is not Galois.

        P can also be an infinite prime, i.e. an embedding into `\RR` or `\CC`.

        EXAMPLE::

            sage: K.<a> = NumberField(x^4 - 2*x^2 + 2,'b').galois_closure()
            sage: P = K.ideal([17, a^2])
            sage: G = K.galois_group()
            sage: G.decomposition_group(P)
            Subgroup [(), (1,8)(2,7)(3,6)(4,5)] of Galois group of Number Field in a with defining polynomial x^8 - 20*x^6 + 104*x^4 - 40*x^2 + 1156
            sage: G.decomposition_group(P^2)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (...) is not prime
            sage: G.decomposition_group(17)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (17) is not prime

        An example with an infinite place::

            sage: L.<b> = NumberField(x^3 - 2,'a').galois_closure(); G=L.galois_group()
            sage: x = L.places()[0]
            sage: G.decomposition_group(x).order()
            2
        """
        if not self.is_galois():
            raise TypeError, "Decomposition groups only defined for Galois extensions"

        if isinstance(P, NumberFieldHomomorphism_im_gens):
            if self.number_field().is_totally_real():
                return self.subgroup([self.identity()])
            else:
                return self.subgroup([self.identity(), self.complex_conjugation(P)])
        else:
            P = self.number_field().ideal_monoid()(P)
            if not P.is_prime():
                raise ValueError, "%s is not prime" % P
            return self.subgroup([s for s in self if s(P) == P])

    def complex_conjugation(self, P=None):
        """
        Return the unique element of self corresponding to complex conjugation,
        for a specified embedding P into the complex numbers. If P is not
        specified, use the "standard" embedding, whenever that is well-defined.

        EXAMPLE::

            sage: L = CyclotomicField(7)
            sage: G = L.galois_group()
            sage: G.complex_conjugation()
            (1,6)(2,3)(4,5)

        An example where the field is not CM, so complex conjugation really
        depends on the choice of embedding::

            sage: L = NumberField(x^6 + 40*x^3 + 1372,'a')
            sage: G = L.galois_group()
            sage: [G.complex_conjugation(x) for x in L.places()]
            [(1,3)(2,6)(4,5), (1,5)(2,4)(3,6), (1,2)(3,4)(5,6)]
        """
        if P is None:
            Q = self.number_field().specified_complex_embedding()
            if Q is None:
                raise ValueError, "No default complex embedding specified"
            P = Q

        P = refine_embedding(P, infinity)

        if not self.number_field().is_galois():
            raise TypeError, "Extension is not Galois"
        if self.number_field().is_totally_real():
            raise TypeError, "No complex conjugation (field is real)"

        g = self.number_field().gen()
        gconj = P(g).conjugate()
        elts = [s for s in self if P(s(g)) == gconj]
        if len(elts) != 1: raise ArithmeticError, "Something has gone very wrong here"
        return elts[0]

    def ramification_group(self, P, v):
        """
        Return the vth ramification group of self for the prime P, i.e. the set
        of elements s of self such that s acts trivially modulo P^(v+1). This
        is only defined for Galois fields.

        EXAMPLE::

            sage: K.<b> = NumberField(x^3 - 3,'a').galois_closure()
            sage: G=K.galois_group()
            sage: P = K.primes_above(3)[0]
            sage: G.ramification_group(P, 3)
            Subgroup [(), (1,3,6)(2,4,5), (1,6,3)(2,5,4)] of Galois group of Number Field in b with defining polynomial x^6 + 60*x^3 + 3087
            sage: G.ramification_group(P, 5)
            Subgroup [()] of Galois group of Number Field in b with defining polynomial x^6 + 60*x^3 + 3087
        """
        if not self.is_galois():
            raise TypeError, "Ramification groups only defined for Galois extensions"
        P = self.number_field().ideal_monoid()(P)
        if not P.is_prime():
            raise ValueError, "%s is not prime"
        return self.subgroup([g for g in self if g(P) == P and g.ramification_degree(P) >= v + 1])

    def inertia_group(self, P):
        """
        Return the inertia group of the prime P, i.e. the group of elements acting
        trivially modulo P. This is just the 0th ramification group of P.

        EXAMPLE::

            sage: K.<b> = NumberField(x^2 - 3,'a')
            sage: G = K.galois_group()
            sage: G.inertia_group(K.primes_above(2)[0])
            Galois group of Number Field in b with defining polynomial x^2 - 3
            sage: G.inertia_group(K.primes_above(5)[0])
            Subgroup [()] of Galois group of Number Field in b with defining polynomial x^2 - 3
        """
        if not self.is_galois():
            raise TypeError, "Inertia groups only defined for Galois extensions"
        return self.ramification_group(P, 0)

    def ramification_breaks(self, P):
        r"""
        Return the set of ramification breaks of the prime ideal P, i.e. the
        set of indices i such that the ramification group `G_{i+1} \ne G_{i}`.
        This is only defined for Galois fields.

        EXAMPLE::

            sage: K.<b> = NumberField(x^8 - 20*x^6 + 104*x^4 - 40*x^2 + 1156)
            sage: G = K.galois_group()
            sage: P = K.primes_above(2)[0]
            sage: G.ramification_breaks(P)
            {1, 3, 5}
            sage: min( [ G.ramification_group(P, i).order() / G.ramification_group(P, i+1).order() for i in G.ramification_breaks(P)] )
            2
        """
        if not self.is_galois():
            raise TypeError, "Ramification breaks only defined for Galois extensions"
        from sage.rings.infinity import infinity
        from sage.sets.set import Set
        i = [g.ramification_degree(P) - 1 for g in self.decomposition_group(P)]
        i.remove(infinity)
        return Set(i)

    def artin_symbol(self, P):
        r"""
        Return the Artin symbol `\left(\frac{K /
        \QQ}{\mathfrak{P}}\right)`, where K is the number field of self,
        and `\mathfrak{P}` is an unramified prime ideal. This is the unique
        element s of the decomposition group of `\mathfrak{P}` such that `s(x) = x^p \bmod
        \mathfrak{P}`, where p is the residue characteristic of `\mathfrak{P}`.

        EXAMPLES::

            sage: K.<b> = NumberField(x^4 - 2*x^2 + 2, 'a').galois_closure()
            sage: G = K.galois_group()
            sage: [G.artin_symbol(P) for P in K.primes_above(7)]
            [(1,5)(2,6)(3,7)(4,8), (1,5)(2,6)(3,7)(4,8), (1,4)(2,3)(5,8)(6,7), (1,4)(2,3)(5,8)(6,7)]
            sage: G.artin_symbol(17)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (17) is not prime
            sage: QuadraticField(-7,'c').galois_group().artin_symbol(13)
            (1,2)
            sage: G.artin_symbol(K.primes_above(2)[0])
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (...) is ramified
        """
        if not self.is_galois():
            raise TypeError, "Artin symbols only defined for Galois extensions"

        P = self.number_field().ideal_monoid()(P)
        if not P.is_prime():
            raise ValueError, "%s is not prime" % P
        p = P.smallest_integer()
        t = []
        gens = self.number_field().ring_of_integers().ring_generators()
        for s in self.decomposition_group(P):
            w = [ (s(g) - g**p).valuation(P) for g in gens]
            if min(w) >= 1:
                t.append(s)
        if len(t) > 1: raise ValueError, "%s is ramified" % P
        return t[0]

class GaloisGroup_subgroup(GaloisGroup_v2):
    r"""
    A subgroup of a Galois group, as returned by functions such as ``decomposition_group``.
    """

    def __init__(self, ambient, elts):
        r"""
        Create a subgroup of a Galois group with the given elements. It is generally better to
        use the subgroup() method of the parent group.

        EXAMPLE::

            sage: from sage.rings.number_field.galois_group import GaloisGroup_subgroup
            sage: G = NumberField(x^3 - x - 1,'a').galois_closure('b').galois_group()
            sage: GaloisGroup_subgroup( G, [ G(1), G([(1,5,2),(3,4,6)]), G([(1,2,5),(3,6,4)])])
            Subgroup [(), (1,5,2)(3,4,6), (1,2,5)(3,6,4)] of Galois group of Number Field in b with defining polynomial x^6 - 14*x^4 + 20*x^3 + 49*x^2 - 140*x + 307
        """
        #XXX: This should be fixed so that this can use GaloisGroup_v2.__init__
        PermutationGroup_generic.__init__(self, elts, canonicalize = True)
        self._ambient = ambient
        self._number_field = ambient.number_field()
        self._galois_closure = ambient._galois_closure
        self._pari_data = ambient._pari_data
        self._pari_gc = ambient._pari_gc
        self._gc_map = ambient._gc_map
        self._elts = elts

    def fixed_field(self):
        r"""
        Return the fixed field of this subgroup (as a subfield of the Galois
        closure of the number field associated to the ambient Galois group).

        EXAMPLE::

            sage: L.<a> = NumberField(x^4 + 1)
            sage: G = L.galois_group()
            sage: H = G.decomposition_group(L.primes_above(3)[0])
            sage: H.fixed_field()
            (Number Field in a0 with defining polynomial x^2 + 2, Ring morphism:
            From: Number Field in a0 with defining polynomial x^2 + 2
            To:   Number Field in a with defining polynomial x^4 + 1
            Defn: a0 |--> a^3 + a)

        """
        if self.order() == 1:
            return self._galois_closure # work around a silly error

        vecs = [pari(g.domain()).Vecsmall() for g in self._elts]
        v = self._ambient._pari_data.galoisfixedfield(vecs)
        x = self._galois_closure(v[1])
        return self._galois_closure.subfield(x)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: G = NumberField(x^3 - x - 1,'a').galois_closure('b').galois_group()
            sage: H = G.subgroup([ G(1), G([(1,5,2),(3,4,6)]), G([(1,2,5),(3,6,4)])])
            sage: H._repr_()
            'Subgroup [(), (1,5,2)(3,4,6), (1,2,5)(3,6,4)] of Galois group of Number Field in b with defining polynomial x^6 - 14*x^4 + 20*x^3 + 49*x^2 - 140*x + 307'
        """
        return "Subgroup %s of %s" % (self._elts, self._ambient)

class GaloisGroupElement(PermutationGroupElement):
    r"""
    An element of a Galois group. This is stored as a permutation, but may also
    be made to act on elements of the field (generally returning elements of
    its Galois closure).

    EXAMPLE::

        sage: K.<w> = QuadraticField(-7); G = K.galois_group()
        sage: G[1]
        (1,2)
        sage: G[1](w + 2)
        -w + 2

        sage: L.<v> = NumberField(x^3 - 2); G = L.galois_group(names='y')
        sage: G[1]
        (1,2)(3,4)(5,6)
        sage: G[1](v)
        1/84*y^4 + 13/42*y
        sage: G[1]( G[1](v) )
        -1/252*y^4 - 55/126*y
    """

    @cached_method
    def as_hom(self):
        r"""
        Return the homomorphism L -> L corresponding to self, where L is the
        Galois closure of the ambient number field.

        EXAMPLE::

            sage: G = QuadraticField(-7,'w').galois_group()
            sage: G[1].as_hom()
            Ring endomorphism of Number Field in w with defining polynomial x^2 + 7
            Defn: w |--> -w
        """
        L = self.parent().splitting_field()
        a = L(self.parent()._pari_data.galoispermtopol(pari(self.domain()).Vecsmall()))
        return L.hom(a, L)

    def __call__(self, x):
        r"""
        Return the action of self on an element x in the number field of self
        (or its Galois closure).

        EXAMPLE::

            sage: K.<w> = QuadraticField(-7)
            sage: f = K.galois_group()[1]
            sage: f(w)
            -w
        """
        if x.parent() == self.parent().splitting_field():
            return self.as_hom()(x)
        else:
            return self.as_hom()( self.parent()._gc_map(x) )

    def ramification_degree(self, P):
        """
        Return the greatest value of v such that s acts trivially modulo P^v.
        Should only be used if P is prime and s is in the decomposition group of P.

        EXAMPLE::

            sage: K.<b> = NumberField(x^3 - 3,'a').galois_closure()
            sage: G=K.galois_group()
            sage: P = K.primes_above(3)[0]
            sage: s = hom(K, K, 1/54*b^4 + 1/18*b)
            sage: G(s).ramification_degree(P)
            4
        """
        if not self.parent().is_galois():
            raise TypeError, "Ramification degree only defined for Galois extensions"
        gens = self.parent().number_field().ring_of_integers().ring_generators()
        w = [ (self(g) - g).valuation(P) for g in gens]
        return min(w)

    def __cmp__(self, other):
        r"""
        Compare self to other. For some bizarre reason, if you just let it
        inherit the cmp routine from PermutationGroupElement, cmp(x, y) works
        but sorting lists doesn't.

        TEST::

            sage: K.<a> = NumberField(x^6 + 40*x^3 + 1372);G = K.galois_group()
            sage: sorted([G.artin_symbol(Q) for Q in K.primes_above(5)])
            [(1,3)(2,6)(4,5), (1,2)(3,4)(5,6), (1,5)(2,4)(3,6)]
        """
        return PermutationGroupElement.__cmp__(self, other)



# For unpickling purposes we rebind GaloisGroup as GaloisGroup_v1.

GaloisGroup = GaloisGroup_v1

