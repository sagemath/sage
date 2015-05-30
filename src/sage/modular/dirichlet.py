# -*- coding: utf-8 -*-
r"""
Dirichlet characters

A ``DirichletCharacter`` is the extension of a
homomorphism

.. math::

       (\ZZ/N\ZZ)^* \to R^*,


for some ring `R`, to the map
`\ZZ/N\ZZ \to R` obtained by sending those
`x\in\ZZ/N\ZZ` with `\gcd(N,x)>1` to
`0`.

EXAMPLES::

    sage: G = DirichletGroup(35)
    sage: x = G.gens()
    sage: e = x[0]*x[1]^2; e
    Dirichlet character modulo 35 of conductor 35 mapping 22 |--> zeta12^3, 31 |--> zeta12^2 - 1
    sage: e.order()
    12

This illustrates a canonical coercion.

::

    sage: e = DirichletGroup(5, QQ).0
    sage: f = DirichletGroup(5,CyclotomicField(4)).0
    sage: e*f
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -zeta4

AUTHORS:

-  William Stein (2005-09-02): Fixed bug in comparison of Dirichlet
   characters. It was checking that their values were the same, but
   not checking that they had the same level!

-  William Stein (2006-01-07): added more examples

-  William Stein (2006-05-21): added examples of everything; fix a
   *lot* of tiny bugs and design problem that became clear when
   creating examples.

-  Craig Citro (2008-02-16): speed up __call__ method for
   Dirichlet characters, miscellaneous fixes

-  Julian Rueth (2014-03-06): use UniqueFactory to cache DirichletGroups

"""

########################################################################
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
########################################################################

import weakref

import sage.categories.all                  as cat
from sage.misc.all import prod
import sage.misc.prandom                    as random
import sage.modules.free_module             as free_module
import sage.modules.free_module_element     as free_module_element
import sage.rings.all                       as rings
import sage.rings.arith                     as arith
import sage.rings.number_field.number_field as number_field
import sage.structure.parent_gens           as parent_gens

from sage.rings.rational_field import is_RationalField
from sage.rings.complex_field import is_ComplexField
from sage.rings.ring import is_Ring

from sage.misc.cachefunc                    import cached_method
from sage.rings.arith                       import binomial, bernoulli
from sage.structure.element                 import MultiplicativeGroupElement
from sage.structure.sequence                import Sequence
from sage.structure.factory                 import UniqueFactory

def trivial_character(N, base_ring=rings.RationalField()):
    r"""
    Return the trivial character of the given modulus, with values in the given
    base ring.

    EXAMPLE::

        sage: t = trivial_character(7)
        sage: [t(x) for x in [0..20]]
        [0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1]
        sage: t(1).parent()
        Rational Field
        sage: trivial_character(7, Integers(3))(1).parent()
        Ring of integers modulo 3
    """
    return DirichletGroup(N, base_ring)(1)

TrivialCharacter = trivial_character

def kronecker_character(d):
    """
    Returns the quadratic Dirichlet character (d/.) of minimal
    conductor.

    EXAMPLES::

        sage: kronecker_character(97*389*997^2)
        Dirichlet character modulo 37733 of conductor 37733 mapping 1557 |--> -1, 37346 |--> -1

    ::

        sage: a = kronecker_character(1)
        sage: b = DirichletGroup(2401,QQ)(a)    # NOTE -- over QQ!
        sage: b.modulus()
        2401

    AUTHORS:

    - Jon Hanke (2006-08-06)
    """
    d = rings.Integer(d)
    if d == 0:
        raise ValueError("d must be nonzero")

    D = arith.fundamental_discriminant(d)
    G = DirichletGroup(abs(D), rings.RationalField())
    return G([arith.kronecker(D,u) for u in G.unit_gens()])


def kronecker_character_upside_down(d):
    """
    Returns the quadratic Dirichlet character (./d) of conductor d, for
    d0.

    EXAMPLES::

        sage: kronecker_character_upside_down(97*389*997^2)
        Dirichlet character modulo 37506941597 of conductor 37733 mapping 13533432536 |--> -1, 22369178537 |--> -1, 14266017175 |--> 1

    AUTHORS:

    - Jon Hanke (2006-08-06)
    """
    d = rings.Integer(d)
    if d <= 0:
        raise ValueError("d must be positive")

    G = DirichletGroup(d, rings.RationalField())
    return G([arith.kronecker(u.lift(),d) for u in G.unit_gens()])


def is_DirichletCharacter(x):
    r"""
    Return True if x is of type DirichletCharacter.

    EXAMPLES::

        sage: from sage.modular.dirichlet import is_DirichletCharacter
        sage: is_DirichletCharacter(trivial_character(3))
        True
        sage: is_DirichletCharacter([1])
        False
    """
    return isinstance(x, DirichletCharacter)

class DirichletCharacter(MultiplicativeGroupElement):
    """
    A Dirichlet character
    """
    def __init__(self, parent, x, check=True):
        r"""
        Create with ``DirichletCharacter(parent, values_on_gens)``

        INPUT:


        -  ``parent`` - DirichletGroup, a group of Dirichlet
           characters

        -  ``x``

           - tuple (or list) of ring elements, the values of the
             Dirichlet character on the chosen generators of
             `(\ZZ/N\ZZ)^*`.

           - Vector over Z/eZ, where e is the order of the root of
             unity.


        OUTPUT:


        -  ``DirichletCharacter`` - a Dirichlet character


        EXAMPLES::

            sage: G.<e> = DirichletGroup(13)
            sage: G
            Group of Dirichlet characters of modulus 13 over Cyclotomic Field of order 12 and degree 4
            sage: e
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12
            sage: loads(e.dumps()) == e
            True

        ::

            sage: G, x = DirichletGroup(35).objgens()
            sage: e = x[0]*x[1]; e
            Dirichlet character modulo 35 of conductor 35 mapping 22 |--> zeta12^3, 31 |--> zeta12^2
            sage: e.order()
            12
            sage: loads(e.dumps()) == e
            True
        """
        MultiplicativeGroupElement.__init__(self, parent)
        self.__modulus = parent.modulus()
        if check:
            if len(x) != len(parent.unit_gens()):
                raise ValueError("wrong number of values(=%s) on unit gens (want %s)"%( \
                       x,len(parent.unit_gens())))
            if free_module_element.is_FreeModuleElement(x):
                self.__element = parent._module(x)
            else:
                R = parent.base_ring()
                self.__values_on_gens = tuple([R(z) for z in x])
        else:
            if free_module_element.is_FreeModuleElement(x):
                self.__element = x
            else:
                self.__values_on_gens = x


    @cached_method
    def __eval_at_minus_one(self):
        r"""
        Efficiently evaluate the character at -1 using knowledge of its
        order. This is potentially much more efficient than computing the
        value of -1 directly using dlog and a large power of the image root
        of unity.

        We use the following. Proposition: Suppose eps is a character mod
        `p^n`, where `p` is a prime. Then
        `\varepsilon(-1) = -1` if and only if `p = 2` and
        the factor of eps at 4 is nontrivial or `p > 2` and 2 does
        not divide `\phi(p^n)/\mbox{\rm ord}(\varepsilon)`.

        EXAMPLES::

            sage: chi = DirichletGroup(20).0; chi._DirichletCharacter__eval_at_minus_one()
            -1
        """
        D = self.decomposition()
        val = self.base_ring()(1)
        for e in D:
            if e.modulus() % 2 == 0:
                if e.modulus() % 4 == 0:
                    val *= e.values_on_gens()[0] # first gen is -1 for 2-power modulus
            elif (arith.euler_phi(e.parent().modulus()) / e.order()) % 2 != 0:
                val *= -1
        return val

    def __call__(self, m):
        """
        Return the value of this character at the integer `m`.

        .. warning::

           A table of values of the character is made the first time
           you call this (unless `m` equals -1)

        EXAMPLES::

            sage: G = DirichletGroup(60)
            sage: e = prod(G.gens(), G(1))
            sage: e
            Dirichlet character modulo 60 of conductor 60 mapping 31 |--> -1, 41 |--> -1, 37 |--> zeta4
            sage: e(-1)
            -1
            sage: e(2)
            0
            sage: e(7)
            -zeta4
            sage: Integers(60).unit_gens()
            (31, 41, 37)
            sage: e(31)
            -1
            sage: e(41)
            -1
            sage: e(37)
            zeta4
            sage: e(31*37)
            -zeta4
            sage: parent(e(31*37))
            Cyclotomic Field of order 4 and degree 2
        """
        m = int(m % self.__modulus)
        if self.values.is_in_cache() or m != self.__modulus - 1:
            return self.values()[m]
        else:
            return self.__eval_at_minus_one()

    def change_ring(self, R):
        """
        Returns the base extension of self to the ring R.

        EXAMPLE::

            sage: e = DirichletGroup(7, QQ).0
            sage: f = e.change_ring(QuadraticField(3, 'a'))
            sage: f.parent()
            Group of Dirichlet characters of modulus 7 over Number Field in a with defining polynomial x^2 - 3

        ::

            sage: e = DirichletGroup(13).0
            sage: e.change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: cannot coerce element of order 12 into self
        """
        if self.base_ring() is R:
            return self
        G = self.parent().change_ring(R)
        return G(self)

    def __cmp__(self, other):
        """
        Compare self to other. Note that this only gets called when the parents
        of self and other are identical, via a canonical coercion map; this
        means that characters of different moduli compare as unequal, even if
        they define identical functions on ZZ.

        EXAMPLES::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: f = e.restrict(8)
            sage: e == e
            True
            sage: f == f
            True
            sage: e == f
            False
            sage: k = DirichletGroup(7)([-1])
            sage: k == e
            False
        """
        return cmp(self.element(), other.element())

    def __hash__(self):
        """
        EXAMPLES::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: hash(e)
            1498523633                  # 32-bit
            3713082714464823281         # 64-bit
        """
        return self.element()._hash()

    def __invert__(self):
        """
        Return the multiplicative inverse of self. The notation is self.

        EXAMPLES::

            sage: e = DirichletGroup(13).0
            sage: f = ~e
            sage: f*e
            Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1
        """
        return DirichletCharacter(self.parent(), -self.element(), check=False)

    def _mul_(self,  other):
        """
        Return the product of self and other.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: b
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4
            sage: a*b # indirect doctest
            Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4

        Multiplying elements whose parents have different zeta orders works::

            sage: a = DirichletGroup(3, QQ, zeta=1, zeta_order=1)(1)
            sage: b = DirichletGroup(3, QQ, zeta=-1, zeta_order=2)([-1])
            sage: a * b # indirect doctest
            Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1
        """
        x = self.element() + other.element()
        return DirichletCharacter(self.parent(), x, check=False)

    #values_on_gens = [self.__values_on_gens[i]*other.__values_on_gens[i]
    #for i in range(len(self.__values_on_gens))]
    #    return DirichletCharacter(self.parent(), values_on_gens)

##     def x_mul_(self,  other):
##         """
##         Return the product of self and other.

##         EXAMPLES:
##             sage: G.<a,b> = DirichletGroup(20)
##             sage: a
##             [-1, 1]
##             sage: b
##             [1, zeta4]
##             sage: a*b
##             [-1, zeta4]
##         """
##         values_on_gens = [self.__values_on_gens[i]*other.__values_on_gens[i]
##                           for i in range(len(self.__values_on_gens))]
##         return DirichletCharacter(self.parent(), values_on_gens)
##         P = self.parent()
##         dlog = P._zeta_dlog
##         pows = P._zeta_powers
##         n = len(pows)
##         values_on_gens = [None]*len(self.__values_on_gens)
##         for i in range(len(self.__values_on_gens)):
##             k = (dlog[self.__values_on_gens[i]] + dlog[other.__values_on_gens[i]]) % n
##             values_on_gens[i] = pows[k]
##         return DirichletCharacter(self.parent(), values_on_gens)

    def __copy__(self):
        """
        Return a (shallow) copy of this Dirichlet character.

        EXAMPLE::

            sage: G.<a> = DirichletGroup(11)
            sage: b = copy(a)
            sage: a is b
            False
            sage: a.element() is b.element()
            True
        """
        # This method exists solely because of a bug in the cPickle module --
        # see modsym/manin_symbols.py.
        return DirichletCharacter(self.parent(), self.element(), check=False)

    def __pow__(self, n):
        """
        Return self raised to the power of n

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a^2
            Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1
            sage: b^2
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1
        """
        return DirichletCharacter(self.parent(), n * self.element(), check=False)

    def _repr_short_(self):
        r"""
        A short string representation of self, often used in string representations of modular forms

        EXAMPLES::

            sage: chi = DirichletGroup(24).0
            sage: chi._repr_short_()
            '[-1, 1, 1]'

        """
        return str(list(self.values_on_gens()))

    def _repr_(self):
        """
        String representation of self.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: repr(a) # indirect doctest
            'Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1'

        TESTS:

        Dirichlet characters modulo 1 and 2 are printed correctly (see
        :trac:`17338`)::

            sage: DirichletGroup(1)[0]
            Dirichlet character modulo 1 of conductor 1
            sage: DirichletGroup(2)[0]
            Dirichlet character modulo 2 of conductor 1
        """
        s = 'Dirichlet character modulo %s of conductor %s' % (self.modulus(), self.conductor())
        r = len(self.values_on_gens())
        if r != 0:
            s += ' mapping '
        for i in range(r):
            if i != 0:
                s += ', '
            s += str(self.parent().unit_gens()[i]) + ' |--> ' + str(self.values_on_gens()[i])
        return s

    def _latex_(self):
        r"""
        LaTeX representation of self.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(16)
            sage: latex(b)  # indirect doctest
            \hbox{Dirichlet character modulo } 16 \hbox{ of conductor } 16 \hbox{ mapping } 15 \mapsto 1,\ 5 \mapsto \zeta_{4}

        TESTS:

        Dirichlet characters modulo 1 and 2 are printed correctly (see
        :trac:`17338`)::

            sage: latex(DirichletGroup(1)[0])
            \hbox{Dirichlet character modulo } 1 \hbox{ of conductor } 1
            sage: latex(DirichletGroup(2)[0])
            \hbox{Dirichlet character modulo } 2 \hbox{ of conductor } 1
        """
        from sage.misc.latex import latex
        s = r'\hbox{Dirichlet character modulo } %s \hbox{ of conductor } %s' % (self.modulus(), self.conductor())
        r = len(self.values_on_gens())
        if r != 0:
            s += r' \hbox{ mapping } '
        for i in range(r):
            if i != 0:
                s += r',\ '
            s += self.parent().unit_gens()[i]._latex_() + r' \mapsto ' + self.values_on_gens()[i]._latex_()
        return s

    def base_ring(self):
        """
        Returns the base ring of this Dirichlet character.

        EXAMPLES::

            sage: G = DirichletGroup(11)
            sage: G.gen(0).base_ring()
            Cyclotomic Field of order 10 and degree 4
            sage: G = DirichletGroup(11, RationalField())
            sage: G.gen(0).base_ring()
            Rational Field
        """
        return self.parent().base_ring()

    def bar(self):
        """
        Return the complex conjugate of this Dirichlet character.

        EXAMPLES::

            sage: e = DirichletGroup(5).0
            sage: e
            Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4
            sage: e.bar()
            Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -zeta4
        """
        return ~self

    def bernoulli(self, k, algorithm='recurrence', cache=True, **opts):
        r"""
        Returns the generalized Bernoulli number `B_{k,eps}`.

        INPUT:


        -  ``k`` - an integer

        -  ``algorithm`` - string (default: 'recurrence');
           either 'recurrence' or 'definition'. The 'recurrence' algorithm
           expresses generalized Bernoulli numbers in terms of classical
           Bernoulli numbers using a recurrence formula and is usually
           optimal. In this case ``**opts`` is passed onto the bernoulli
           function.

        -  ``cache`` - if True, cache answers


        Let eps be this character (not necessarily primitive), and let
        `k \geq 0` be an integer weight. This function computes the
        (generalized) Bernoulli number `B_{k,eps}`, e.g., as
        defined on page 44 of Diamond-Im:

        .. math::

                       \sum_{a=1}^{N} \varepsilon(a) t*e^{at} / (e^{Nt}-1) = sum_{k=0}^{\infty} B_{k,eps}/{k!} t^k.


        where `N` is the modulus of `\varepsilon`.

        The default algorithm is the recurrence on page 656 of Cohen's GTM
        'Number Theory and Diophantine Equations', section 9.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.bernoulli(5)
            7430/13*zeta12^3 - 34750/13*zeta12^2 - 11380/13*zeta12 + 9110/13
            sage: eps = DirichletGroup(9).0
            sage: eps.bernoulli(3)
            10*zeta6 + 4
            sage: eps.bernoulli(3, algorithm="definition")
            10*zeta6 + 4
        """
        if cache:
            try:
                self.__bernoulli
            except AttributeError:
                self.__bernoulli = {}
            if k in self.__bernoulli:
                return self.__bernoulli[k]
        N = self.modulus()
        K = self.base_ring()

        if N != 1 and self(-1) != K((-1)**k):
            ans = K(0)
            if cache:
                self.__bernoulli[k] = ans
            return ans

        if algorithm == "recurrence":
            # The following code is pretty fast, at least compared to
            # the other algorithm below.  That said, I'm sure it could
            # be sped up by a factor of 10 or more in many cases,
            # especially since we end up computing all the bernoulli
            # numbers up to k, which should be done with power series
            # instead of calls to the bernoulli function.  Likewise
            # computing all binomial coefficients can be done much
            # more efficiently.
            if self.modulus() == 1:
                ber = K(bernoulli(k))
            else:
                v = self.values()
                S = lambda n: sum(v[r] * r**n for r in range(1, N))
                ber = K(sum(binomial(k,j) * bernoulli(j, **opts) *
                                    N**(j-1) * S(k-j) for j in range(k+1)))

        elif algorithm == "definition":
            # This is better since it computes the same thing, but requires
            # no arith in a poly ring over a number field.
            prec = k+2
            R = rings.PowerSeriesRing(rings.QQ, 't')
            t = R.gen()
            # g(t) = t/(e^{Nt}-1)
            g = t/((N*t).exp(prec) - 1)

            # h(n) = g(t)*e^{nt}
            h = [0] + [g * ((n*t).exp(prec)) for n in range(1,N+1)]
            ber = sum([self(a)*h[a][k] for a in range(1,N+1)]) * arith.factorial(k)

        else:
            raise ValueError("algorithm = '%s' unknown"%algorithm)

        if cache:
            self.__bernoulli[k] = ber
        return ber

    @cached_method
    def conductor(self):
        """
        Computes and returns the conductor of this character.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.conductor()
            4
            sage: b.conductor()
            5
            sage: (a*b).conductor()
            20

        TESTS::

            sage: G.<a, b> = DirichletGroup(20)
            sage: type(G(1).conductor())
            <type 'sage.rings.integer.Integer'>
        """
        if self.modulus() == 1 or self.is_trivial():
            return rings.Integer(1)
        F = arith.factor(self.modulus())
        if len(F) > 1:
            return prod([d.conductor() for d in self.decomposition()])
        p = F[0][0]
        # When p is odd, and x =/= 1, the conductor is the smallest p**r such that
        #   Order(x) divides EulerPhi(p**r) = p**(r-1)*(p-1).
        # For a given r, whether or not the above divisibility holds
        # depends only on the factor of p**(r-1) on the right hand side.
        # Since p-1 is coprime to p, this smallest r such that the
        # divisibility holds equals Valuation(Order(x),p)+1.
        cond = p**(arith.valuation(self.order(),p) + 1)
        if p == 2 and F[0][1] > 2 and self.values_on_gens()[1].multiplicative_order() != 1:
            cond *= 2;
        return rings.Integer(cond)

    @cached_method
    def decomposition(self):
        """
        Return the decomposition of self as a product of Dirichlet
        characters of prime power modulus, where the prime powers exactly
        divide the modulus of this character.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: c = a*b
            sage: d = c.decomposition(); d
            [Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4]
            sage: d[0].parent()
            Group of Dirichlet characters of modulus 4 over Cyclotomic Field of order 4 and degree 2
            sage: d[1].parent()
            Group of Dirichlet characters of modulus 5 over Cyclotomic Field of order 4 and degree 2

        We can't multiply directly, since coercion of one element into the
        other parent fails in both cases::

            sage: d[0]*d[1] == c
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Group of Dirichlet characters of modulus 4 over Cyclotomic Field of order 4 and degree 2' and 'Group of Dirichlet characters of modulus 5 over Cyclotomic Field of order 4 and degree 2'

        We can multiply if we're explicit about where we want the
        multiplication to take place.

        ::

            sage: G(d[0])*G(d[1]) == c
            True

        Conductors that are divisible by various powers of 2 present
        some problems as the multiplicative group modulo `2^k` is
        trivial for `k = 1` and non-cyclic for `k \ge 3`::

            sage: (DirichletGroup(18).0).decomposition()
            [Dirichlet character modulo 2 of conductor 1, Dirichlet character modulo 9 of conductor 9 mapping 2 |--> zeta6]
            sage: (DirichletGroup(36).0).decomposition()
            [Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1, Dirichlet character modulo 9 of conductor 1 mapping 2 |--> 1]
            sage: (DirichletGroup(72).0).decomposition()
            [Dirichlet character modulo 8 of conductor 4 mapping 7 |--> -1, 5 |--> 1, Dirichlet character modulo 9 of conductor 1 mapping 2 |--> 1]
        """
        D = self.parent().decomposition()
        vals = [[z] for z in self.values_on_gens()]
        R = self.base_ring()
        if self.modulus() % 8 == 0:   # 2 factors at 2.
            vals[0].append(vals[1][0])
            del vals[1]
        elif self.modulus() % 4 == 2: # 0 factors at 2.
            vals = [1] + vals
        return [D[i](vals[i]) for i in range(len(D))]

    def extend(self, M):
        """
        Returns the extension of this character to a Dirichlet character
        modulo the multiple M of the modulus.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: H.<c> = DirichletGroup(4)
            sage: c.extend(20)
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: a
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: c.extend(20) == a
            True
        """
        if M % self.modulus() != 0:
            raise ArithmeticError("M(=%s) must be a multiple of the modulus(=%s)"%(M,self.modulus()))
        H = DirichletGroup(M, self.base_ring())
        return H(self)

    def galois_orbit(self, sort=True):
        r"""
        Return the orbit of this character under the action of the absolute
        Galois group of the prime subfield of the base ring.

        EXAMPLES::

            sage: G = DirichletGroup(30); e = G.1
            sage: e.galois_orbit()
            [Dirichlet character modulo 30 of conductor 5 mapping 11 |--> 1, 7 |--> zeta4, Dirichlet character modulo 30 of conductor 5 mapping 11 |--> 1, 7 |--> -zeta4]

        Another example::

            sage: G = DirichletGroup(13)
            sage: G.galois_orbits()
            [
            [Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1], ... [Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1]
            ]
            sage: e = G.0
            sage: e
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12
            sage: e.galois_orbit()
            [Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12^3 - zeta12, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta12, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta12^3 + zeta12]
            sage: e = G.0^2; e
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12^2
            sage: e.galois_orbit()
            [Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12^2, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta12^2 + 1]

        A non-example::

            sage: chi = DirichletGroup(7, Integers(9), zeta = Integers(9)(2)).0
            sage: chi.galois_orbit()
            Traceback (most recent call last):
            ...
            TypeError: Galois orbits only defined if base ring is an integral domain
        """
        if not self.base_ring().is_integral_domain():
            raise TypeError("Galois orbits only defined if base ring is an integral domain")
        k = self.order()
        if k <= 2:
            return [self]
        P = self.parent()
        z = self.element()
        o = int(z.additive_order())
        Auts = set([m % o for m in P._automorphisms()])
        v = [DirichletCharacter(P, m * z, check=False) for m in Auts]
        if sort:
            v.sort()
        return v


    def gauss_sum(self, a=1):
        r"""
        Return a Gauss sum associated to this Dirichlet character.

        The Gauss sum associated to `\chi` is

        .. math::

                           g_a(\chi) = \sum_{r \in \ZZ/m\ZZ} \chi(r)\,\zeta^{ar},


        where `m` is the modulus of `\chi` and
        `\zeta` is a primitive `m^{th}` root of unity, i.e.,
        `\zeta` is ``self.parent().zeta()``.

        FACTS: If the modulus is a prime `p` and the character is
        nontrivial, then the Gauss sum has absolute value
        `\sqrt{p}`.

        CACHING: Computed Gauss sums are *not* cached with this character.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G([-1])
            sage: e.gauss_sum(1)
            2*zeta6 - 1
            sage: e.gauss_sum(2)
            -2*zeta6 + 1
            sage: norm(e.gauss_sum())
            3

        ::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.gauss_sum()
            -zeta156^46 + zeta156^45 + zeta156^42 + zeta156^41 + 2*zeta156^40 + zeta156^37 - zeta156^36 - zeta156^34 - zeta156^33 - zeta156^31 + 2*zeta156^30 + zeta156^28 - zeta156^24 - zeta156^22 + zeta156^21 + zeta156^20 - zeta156^19 + zeta156^18 - zeta156^16 - zeta156^15 - 2*zeta156^14 - zeta156^10 + zeta156^8 + zeta156^7 + zeta156^6 + zeta156^5 - zeta156^4 - zeta156^2 - 1
            sage: factor(norm(e.gauss_sum()))
            13^24
        """
        G = self.parent()
        K = G.base_ring()
        if is_ComplexField(K):
            return self.gauss_sum_numerical()
        if not (number_field.is_CyclotomicField(K) or is_RationalField(K)):
            raise NotImplementedError("Gauss sums only currently implemented when the base ring is a cyclotomic field, QQ, or a complex field.")
        g = 0
        m = G.modulus()
        L = rings.CyclotomicField(arith.lcm(m,G.zeta_order()))
        zeta = L.gen(0)
        n = zeta.multiplicative_order()
        zeta = zeta ** (n // m)
        if a != 1:
            zeta = zeta**a
        z = 1
        for c in self.values()[1:]:
            z *= zeta
            g += L(c)*z
        return g

    def gauss_sum_numerical(self, prec=53, a=1):
        r"""
        Return a Gauss sum associated to this Dirichlet character as an
        approximate complex number with prec bits of precision.

        INPUT:


        - ``prec`` - integer (default: 53), *bits* of precision

        - ``a`` - integer, as for gauss_sum.


        The Gauss sum associated to `\chi` is

        .. math::

           g_a(\chi) = \sum_{r \in \ZZ/m\ZZ} \chi(r)\,\zeta^{ar},


        where `m` is the modulus of `\chi` and
        `\zeta` is a primitive `m^{th}` root of unity, i.e.,
        `\zeta` is ``self.parent().zeta()``.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G.0
            sage: abs(e.gauss_sum_numerical())
            1.7320508075...
            sage: sqrt(3.0)
            1.73205080756888
            sage: e.gauss_sum_numerical(a=2)
            -...e-15 - 1.7320508075...*I
            sage: e.gauss_sum_numerical(a=2, prec=100)
            4.7331654313260708324703713917e-30 - 1.7320508075688772935274463415*I
            sage: G = DirichletGroup(13)
            sage: H = DirichletGroup(13, CC)
            sage: e = G.0
            sage: f = H.0
            sage: e.gauss_sum_numerical()
            -3.07497205... + 1.8826966926...*I
            sage: f.gauss_sum_numerical()
            -3.07497205... + 1.8826966926...*I
            sage: abs(e.gauss_sum_numerical())
            3.60555127546...
            sage: abs(f.gauss_sum_numerical())
            3.60555127546...
            sage: sqrt(13.0)
            3.60555127546399
        """
        G = self.parent()
        K = G.base_ring()
        if not (number_field.is_CyclotomicField(K) or is_RationalField(K)
                or is_ComplexField(K)):
            raise NotImplementedError("Gauss sums only currently implemented when the base ring is a cyclotomic field, QQ, or a complex field.")

        if is_ComplexField(K):
            phi = lambda t : t
            CC = K
        else:
            phi = K.complex_embedding(prec)
            CC = phi.codomain()

        g = 0
        m = G.modulus()
        zeta = CC.zeta(m)
        if a != 1:
            zeta = zeta**a
        z = 1
        for c in self.values()[1:]:
            z *= zeta
            g += phi(c)*z
        return g

    def jacobi_sum(self, char, check=True):
        """
        Return the Jacobi sum associated to these Dirichlet characters
        (i.e., J(self,char)). This is defined as

        .. math::

            J(\chi, \psi) = \sum_{a \in \ZZ / N\ZZ} \chi(a) \psi(1-a)

        where `\chi` and `\psi` are both characters modulo `N`.

        EXAMPLES::

            sage: D = DirichletGroup(13)
            sage: e = D.0
            sage: f = D[-2]
            sage: e.jacobi_sum(f)
            3*zeta12^2 + 2*zeta12 - 3
            sage: f.jacobi_sum(e)
            3*zeta12^2 + 2*zeta12 - 3
            sage: p = 7
            sage: DP = DirichletGroup(p)
            sage: f = DP.0
            sage: e.jacobi_sum(f)
            Traceback (most recent call last):
            ...
            NotImplementedError: Characters must be from the same Dirichlet Group.

            sage: all_jacobi_sums = [(DP[i].values_on_gens(),DP[j].values_on_gens(),DP[i].jacobi_sum(DP[j])) \
            ...                       for i in range(p-1) for j in range(p-1)[i:]]
            ...
            sage: for s in all_jacobi_sums:
            ...       print s
            ((1,), (1,), 5)
            ((1,), (zeta6,), -1)
            ((1,), (zeta6 - 1,), -1)
            ((1,), (-1,), -1)
            ((1,), (-zeta6,), -1)
            ((1,), (-zeta6 + 1,), -1)
            ((zeta6,), (zeta6,), -zeta6 + 3)
            ((zeta6,), (zeta6 - 1,), 2*zeta6 + 1)
            ((zeta6,), (-1,), -2*zeta6 - 1)
            ((zeta6,), (-zeta6,), zeta6 - 3)
            ((zeta6,), (-zeta6 + 1,), 1)
            ((zeta6 - 1,), (zeta6 - 1,), -3*zeta6 + 2)
            ((zeta6 - 1,), (-1,), 2*zeta6 + 1)
            ((zeta6 - 1,), (-zeta6,), -1)
            ((zeta6 - 1,), (-zeta6 + 1,), -zeta6 - 2)
            ((-1,), (-1,), 1)
            ((-1,), (-zeta6,), -2*zeta6 + 3)
            ((-1,), (-zeta6 + 1,), 2*zeta6 - 3)
            ((-zeta6,), (-zeta6,), 3*zeta6 - 1)
            ((-zeta6,), (-zeta6 + 1,), -2*zeta6 + 3)
            ((-zeta6 + 1,), (-zeta6 + 1,), zeta6 + 2)

        Let's check that trivial sums are being calculated correctly::

            sage: N = 13
            sage: D = DirichletGroup(N)
            sage: g = D(1)
            sage: g.jacobi_sum(g)
            11
            sage: sum([g(x)*g(1-x) for x in IntegerModRing(N)])
            11

        And sums where exactly one character is nontrivial (see trac #6393)::

            sage: G=DirichletGroup(5); X=G.list(); Y=X[0]; Z=X[1]
            sage: Y.jacobi_sum(Z)
            -1
            sage: Z.jacobi_sum(Y)
            -1

        Now let's take a look at a non-prime modulus::

            sage: N = 9
            sage: D = DirichletGroup(N)
            sage: g = D(1)
            sage: g.jacobi_sum(g)
            3

        We consider a sum with values in a finite field::

            sage: g = DirichletGroup(17, GF(9,'a')).0
            sage: g.jacobi_sum(g**2)
            2*a

        TESTS:

        This shows that ticket #6393 has been fixed::

            sage: G = DirichletGroup(5); X = G.list(); Y = X[0]; Z = X[1]
            sage: # Y is trivial and Z is quartic
            sage: sum([Y(x)*Z(1-x) for x in IntegerModRing(5)])
            -1
            sage: # The value -1 above is the correct value of the Jacobi sum J(Y, Z).
            sage: Y.jacobi_sum(Z); Z.jacobi_sum(Y)
            -1
            -1
        """
        if check:
            if self.parent() != char.parent():
                raise NotImplementedError("Characters must be from the same Dirichlet Group.")

        return sum([self(x) * char(1-x) for x in rings.IntegerModRing(self.modulus())])

    def kloosterman_sum(self, a=1,b=0):
        r"""
        Return the "twisted" Kloosterman sum associated to this Dirichlet character.
        This includes Gauss sums, classical Kloosterman sums, Salie sums, etc.

        The Kloosterman sum associated to `\chi` and the integers a,b is

        .. math::

            K(a,b,\chi) = \sum_{r \in (\ZZ/m\ZZ)^\times} \chi(r)\,\zeta^{ar+br^{-1}},

        where `m` is the modulus of `\chi` and `\zeta` is a primitive
        `m` th root of unity. This reduces to to the Gauss sum if `b=0`.

        This method performs an exact calculation and returns an element of a
        suitable cyclotomic field; see also :meth:`.kloosterman_sum_numerical`,
        which gives an inexact answer (but is generally much quicker).

        CACHING: Computed Kloosterman sums are *not* cached with this
        character.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G([-1])
            sage: e.kloosterman_sum(3,5)
            -2*zeta6 + 1
            sage: G = DirichletGroup(20)
            sage: e = G([1 for  u in G.unit_gens()])
            sage: e.kloosterman_sum(7,17)
            -2*zeta20^6 + 2*zeta20^4 + 4

        """
        G = self.parent()
        K = G.base_ring()
        if not (number_field.is_CyclotomicField(K) or is_RationalField(K)):
            raise NotImplementedError("Kloosterman sums only currently implemented when the base ring is a cyclotomic field or QQ.")
        g = 0
        m = G.modulus()
        L = rings.CyclotomicField(arith.lcm(m,G.zeta_order()))
        zeta = L.gen(0)
        n = zeta.multiplicative_order()
        zeta = zeta ** (n // m)
        for c in range(1,m):
            if arith.gcd(c,m)==1:
                e = rings.Mod(c,m)
                z = zeta ** int(a*e + b*(e**(-1)))
                g += self.__call__(c)*z
        return g

    def kloosterman_sum_numerical(self, prec=53, a=1,b=0):
        r"""
        Return the Kloosterman sum associated to this Dirichlet character as
        an approximate complex number with prec bits of precision. See also
        :meth:`.kloosterman_sum`, which calculates the sum exactly (which is
        generally slower).

        INPUT:

        - ``prec`` -- integer (default: 53), *bits* of precision
        - ``a`` -- integer, as for :meth:`.kloosterman_sum`
        - ``b`` -- integer, as for :meth:`.kloosterman_sum`.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G.0

        The real component of the numerical value of e is near zero::

            sage: v=e.kloosterman_sum_numerical()
            sage: v.real() < 1.0e15
            True
            sage: v.imag()
            1.73205080756888
            sage: G = DirichletGroup(20)
            sage: e = G.1
            sage: e.kloosterman_sum_numerical(53,3,11)
            3.80422606518061 - 3.80422606518061*I
        """
        G = self.parent()
        K = G.base_ring()
        if not (number_field.is_CyclotomicField(K) or is_RationalField(K)):
            raise NotImplementedError("Kloosterman sums only currently implemented when the base ring is a cyclotomic field or QQ.")
        phi = K.complex_embedding(prec)
        CC = phi.codomain()
        g = 0
        m = G.modulus()
        zeta = CC.zeta(m)

        for c in range(1,m):
            if arith.gcd(c,m)==1:
                e = rings.Mod(c,m)
                z = zeta ** int(a*e + b*(e**(-1)))
                g += phi(self.__call__(c))*z
        return g

    @cached_method
    def is_even(self):
        r"""
        Return ``True`` if and only if `\varepsilon(-1) = 1`.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.is_even()
            False
            sage: e(-1)
            -1
            sage: [e.is_even() for e in G]
            [True, False, True, False, True, False, True, False, True, False, True, False]

            sage: G = DirichletGroup(13, CC)
            sage: e = G.0
            sage: e.is_even()
            False
            sage: e(-1)
            -1.000000...
            sage: [e.is_even() for e in G]
            [True, False, True, False, True, False, True, False, True, False, True, False]

            sage: G = DirichletGroup(100000, CC)
            sage: G.1.is_even()
            True

        Note that ``is_even`` need not be the negation of
        is_odd, e.g., in characteristic 2::

            sage: G.<e> = DirichletGroup(13, GF(4,'a'))
            sage: e.is_even()
            True
            sage: e.is_odd()
            True
        """
        R = self.base_ring()
        # self(-1) is either +1 or -1
        if not R.is_exact():
            return abs(self(-1) - R(1)) < 0.5
        return self(-1) == R(1)

    @cached_method
    def is_odd(self):
        r"""
        Return ``True`` if and only if
        `\varepsilon(-1) = -1`.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.is_odd()
            True
            sage: [e.is_odd() for e in G]
            [False, True, False, True, False, True, False, True, False, True, False, True]

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.is_odd()
            True
            sage: [e.is_odd() for e in G]
            [False, True, False, True, False, True, False, True, False, True, False, True]

            sage: G = DirichletGroup(100000, CC)
            sage: G.0.is_odd()
            True

        Note that ``is_even`` need not be the negation of
        is_odd, e.g., in characteristic 2::

            sage: G.<e> = DirichletGroup(13, GF(4,'a'))
            sage: e.is_even()
            True
            sage: e.is_odd()
            True
        """
        R = self.base_ring()
        # self(-1) is either +1 or -1
        if not R.is_exact():
            return abs(self(-1) - R(-1)) < 0.5
        return self(-1) == R(-1)

    @cached_method
    def is_primitive(self):
        """
        Return ``True`` if and only if this character is
        primitive, i.e., its conductor equals its modulus.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.is_primitive()
            False
            sage: b.is_primitive()
            False
            sage: (a*b).is_primitive()
            True
            sage: G.<a,b> = DirichletGroup(20, CC)
            sage: a.is_primitive()
            False
            sage: b.is_primitive()
            False
            sage: (a*b).is_primitive()
            True
        """
        return (self.conductor() == self.modulus())

    @cached_method
    def is_trivial(self):
        r"""
        Returns ``True`` if this is the trivial character,
        i.e., has order 1.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.is_trivial()
            False
            sage: (a^2).is_trivial()
            True
        """
        return (self.element() == 0)

    def kernel(self):
        r"""
        Return the kernel of this character.

        OUTPUT: Currently the kernel is returned as a list. This may
        change.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.kernel()
            [1, 9, 13, 17]
            sage: b.kernel()
            [1, 11]
        """
        one = self.base_ring()(1)
        return [x for x in range(self.modulus()) if self(x) == one]

    def maximize_base_ring(self):
        r"""
        Let

        .. math::

           \varepsilon : (\ZZ/N\ZZ)^* \to \QQ(\zeta_n)


        be a Dirichlet character. This function returns an equal Dirichlet
        character

        .. math::

           \chi : (\ZZ/N\ZZ)^* \to \QQ(\zeta_m)


        where `m` is the least common multiple of `n` and
        the exponent of `(\ZZ/N\ZZ)^*`.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20,QQ)
            sage: b.maximize_base_ring()
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1
            sage: b.maximize_base_ring().base_ring()
            Cyclotomic Field of order 4 and degree 2
            sage: DirichletGroup(20).base_ring()
            Cyclotomic Field of order 4 and degree 2
        """
        g = rings.IntegerModRing(self.modulus()).unit_group_exponent()
        if g == 1:
            g = 2
        z = self.base_ring().zeta()
        n = z.multiplicative_order()
        m = arith.LCM(g,n)
        if n == m:
            return self
        K = rings.CyclotomicField(m)
        return self.change_ring(K)

    def minimize_base_ring(self):
        r"""
        Return a Dirichlet character that equals this one, but over as
        small a subfield (or subring) of the base ring as possible.

        .. note::

           This function is currently only implemented when the base
           ring is a number field. It's the identity function in
           characteristic p.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = DirichletGroup(13).0
            sage: e.base_ring()
            Cyclotomic Field of order 12 and degree 4
            sage: e.minimize_base_ring().base_ring()
            Cyclotomic Field of order 12 and degree 4
            sage: (e^2).minimize_base_ring().base_ring()
            Cyclotomic Field of order 6 and degree 2
            sage: (e^3).minimize_base_ring().base_ring()
            Cyclotomic Field of order 4 and degree 2
            sage: (e^12).minimize_base_ring().base_ring()
            Rational Field
        """
        if isinstance(self.base_ring(),rings.RationalField):
            return self

        if self.is_trivial() and self.base_ring().characteristic() == 0:
            return self.change_ring(rings.QQ)

        if isinstance(self.base_ring(),number_field.NumberField_generic):
            if self.order() <= 2:
                return self.change_ring(rings.RationalField())
            if arith.euler_phi(self.order()) == self.base_ring().degree():
                return self
            K = rings.CyclotomicField(self.order())
            return self.change_ring(K)

        #raise NotImplementedError, "minimize_base_ring is currently " + \
        #      "only implemented when the base ring is a number field."
        return self

    def modulus(self):
        """
        The modulus of this character.

        EXAMPLES::

            sage: e = DirichletGroup(100, QQ).0
            sage: e.modulus()
            100
            sage: e.conductor()
            4
        """
        return self.__modulus

    def level(self):
        """
        Synonym for modulus.

        EXAMPLE::

            sage: e = DirichletGroup(100, QQ).0
            sage: e.level()
            100
        """
        return self.modulus()

    @cached_method
    def multiplicative_order(self):
        """
        The order of this character.

        EXAMPLES::

            sage: e = DirichletGroup(100).1
            sage: e.order()    # same as multiplicative_order, since group is multiplicative
            20
            sage: e.multiplicative_order()
            20
            sage: e = DirichletGroup(100).0
            sage: e.multiplicative_order()
            2
        """
        return self.element().additive_order()

    def primitive_character(self):
        """
        Returns the primitive character associated to self.

        EXAMPLES::

            sage: e = DirichletGroup(100).0; e
            Dirichlet character modulo 100 of conductor 4 mapping 51 |--> -1, 77 |--> 1
            sage: e.conductor()
            4
            sage: f = e.primitive_character(); f
            Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1
            sage: f.modulus()
            4
        """
        return self.restrict(self.conductor())

    def restrict(self, M):
        """
        Returns the restriction of this character to a Dirichlet character
        modulo the divisor M of the modulus, which must also be a multiple
        of the conductor of this character.

        EXAMPLES::

            sage: e = DirichletGroup(100).0
            sage: e.modulus()
            100
            sage: e.conductor()
            4
            sage: e.restrict(20)
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: e.restrict(4)
            Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1
            sage: e.restrict(50)
            Traceback (most recent call last):
            ...
            ValueError: conductor(=4) must divide M(=50)
        """
        M = int(M)
        if self.modulus()%M != 0:
            raise ValueError("M(=%s) must divide the modulus(=%s)"%(M,self.modulus()))
        if M%self.conductor() != 0:
            raise ValueError("conductor(=%s) must divide M(=%s)"%(self.conductor(),M))
        H = DirichletGroup(M, self.base_ring())
        return H(self)

    @cached_method
    def values(self):
        """
        Return a list of the values of this character on each integer
        between 0 and the modulus.

        EXAMPLES::

            sage: e = DirichletGroup(20)(1)
            sage: e.values()
            [0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1]
            sage: e = DirichletGroup(20).gen(0)
            sage: print e.values()
            [0, 1, 0, -1, 0, 0, 0, -1, 0, 1, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1]
            sage: e = DirichletGroup(20).gen(1)
            sage: e.values()
            [0, 1, 0, -zeta4, 0, 0, 0, zeta4, 0, -1, 0, 1, 0, -zeta4, 0, 0, 0, zeta4, 0, -1]
            sage: e = DirichletGroup(21).gen(0) ; e.values()
            [0, 1, -1, 0, 1, -1, 0, 0, -1, 0, 1, -1, 0, 1, 0, 0, 1, -1, 0, 1, -1]
            sage: e = DirichletGroup(21, base_ring=GF(37)).gen(0) ; e.values()
            [0, 1, 36, 0, 1, 36, 0, 0, 36, 0, 1, 36, 0, 1, 0, 0, 1, 36, 0, 1, 36]
            sage: e = DirichletGroup(21, base_ring=GF(3)).gen(0) ; e.values()
            [0, 1, 2, 0, 1, 2, 0, 0, 2, 0, 1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 2]

        ::

            sage: chi = DirichletGroup(100151, CyclotomicField(10)).0
            sage: ls = chi.values() ; ls[0:10]
            [0,
            1,
            -zeta10^3,
            -zeta10,
            -zeta10,
            1,
            zeta10^3 - zeta10^2 + zeta10 - 1,
            zeta10,
            zeta10^3 - zeta10^2 + zeta10 - 1,
            zeta10^2]

        TESTS:

        Test that :trac:`11783` and :trac:`14368` are fixed::

            sage: chi = DirichletGroup(1).list()[0]
            sage: chi.values()
            [1]
            sage: chi(1)
            1
        """
        G = self.parent()
        R = G.base_ring()

        mod = self.__modulus
        if mod == 1:
            return [R.one()]
        elif mod == 2:
            return [R.zero(), R.one()]

        result_list = [R.zero()] * mod
        gens = G.unit_gens()
        orders = G.integers_mod().unit_group().gens_orders()

        R_values = G._zeta_powers
        val_on_gen = self.element()

        exponents = [0] * len(orders)
        n = G.integers_mod().one()
        value = val_on_gen.base_ring().zero()

        while True:
            # record character value on n
            result_list[n] = R_values[value]
            # iterate:
            #   increase the exponent vector by 1,
            #   increase n accordingly, and increase value
            i = 0
            while True:
                try:
                    exponents[i] += 1
                except IndexError:  # Done!
                    return result_list
                value += val_on_gen[i]
                n *= gens[i]
                if exponents[i] < orders[i]:
                    break
                exponents[i] = 0
                i += 1


    def values_on_gens(self):
        r"""
        Return a tuple of the values of ``self`` on the standard
        generators of `(\ZZ/N\ZZ)^*`, where `N` is the modulus.

        EXAMPLES::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: e.values_on_gens ()
            (-1, 1)
        """
        try:
            return self.__values_on_gens
        except AttributeError:
            pows = self.parent()._zeta_powers
            v = tuple([pows[i] for i in self.element()])
            self.__values_on_gens = v
            return v

    def element(self):
        r"""
        Return the underlying `\ZZ/n\ZZ`-module
        vector of exponents.

        .. warning::

           Please do not change the entries of the returned vector;
           this vector is mutable *only* because immutable vectors are
           implemented yet.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.element()
            (2, 0)
            sage: b.element()
            (0, 1)
        """
        try:
            return self.__element
        except AttributeError:
            P = self.parent()
            M = P._module
            C = P.base_ring()
            if is_ComplexField(C):
                R = C._real_field()
                zeta = P._zeta
                zeta_argument = zeta.argument()
                v = M([int(round(x.argument()/zeta_argument))
                       for x in self.values_on_gens()])
            else:
                dlog = P._zeta_dlog
                v = M([dlog[x] for x in self.values_on_gens()])
            self.__element = v
            return v


class DirichletGroupFactory(UniqueFactory):
    r"""
    The group of Dirichlet characters modulo `N` with values in the subgroup
    `\langle \zeta_n\rangle` of the multiplicative group of the ``base_ring``.
    If ``base_ring`` is omitted then we use `\QQ(\zeta_n)`, where `n` is the
    exponent of `(\ZZ/N\ZZ)^*`. If `\zeta` is omitted then we compute and use a
    maximal-order zeta in ``base_ring``, if possible.

    INPUT:

    - ``N`` - an integer

    - ``base_ring`` - a ring (optional), where characters take their values
      (should be an integral domain)

    - ``zeta`` - an element (optional), a root of unity in ``base_ring``

    - ``zeta_order`` - an integer (optional), the order of zeta

    - ``names`` - ignored (needed so ``G.... = DirichletGroup(...)`` notation
      works)

    - ``integral`` - a boolean (default: ``False``). If ``True``, return the
      group with ``base_ring`` the ring of integers in the smallest choice of
      :meth:`sage.rings.number_field.number_field.CyclotomicField`. Ignored if
      ``base_ring`` is not ``None``.

    OUTPUT:

    a group of Dirichlet characters

    EXAMPLES:

    The default base ring is a cyclotomic field of order the exponent
    of `(\ZZ/N\ZZ)^*`::

        sage: DirichletGroup(20)
        Group of Dirichlet characters of modulus 20 over Cyclotomic Field of order 4 and degree 2

    We create the group of Dirichlet character mod 20 with values in
    the rational numbers::

        sage: G = DirichletGroup(20, QQ); G
        Group of Dirichlet characters of modulus 20 over Rational Field
        sage: G.order()
        4
        sage: G.base_ring()
        Rational Field

    The elements of G print as lists giving the values of the character
    on the generators of `(Z/NZ)^*`::

        sage: list(G)
        [Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1, Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1, Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1, Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -1]

    Next we construct the group of Dirichlet character mod 20, but with
    values in `\QQ(\zeta_n)`::

        sage: G = DirichletGroup(20)
        sage: G.1
        Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4

    We next compute several invariants of ``G``::

        sage: G.gens()
        (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1, Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4)
        sage: G.unit_gens()
        (11, 17)
        sage: G.zeta()
        zeta4
        sage: G.zeta_order()
        4

    In this example we create a Dirichlet character with values in a
    number field. We have to give ``zeta``, but not its order::

        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(x^4 + 1)
        sage: G = DirichletGroup(5, K, a); G
        Group of Dirichlet characters of modulus 5 over Number Field in a with defining polynomial x^4 + 1
        sage: G.list()
        [Dirichlet character modulo 5 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -a^2]

    ::

        sage: G.<e> = DirichletGroup(13)
        sage: loads(G.dumps()) == G
        True

    ::

        sage: G = DirichletGroup(19, GF(5))
        sage: loads(G.dumps()) == G
        True

    We compute a Dirichlet group over a large prime field::

        sage: p = next_prime(10^40)
        sage: g = DirichletGroup(19, GF(p)); g
        Group of Dirichlet characters of modulus 19 over Finite Field of size 10000000000000000000000000000000000000121

    Note that the root of unity has small order, i.e., it is not the
    largest order root of unity in the field::

        sage: g.zeta_order()
        2

    ::

        sage: r4 = CyclotomicField(4).ring_of_integers()
        sage: G = DirichletGroup(60, r4)
        sage: G.gens()
        (Dirichlet character modulo 60 of conductor 4 mapping 31 |--> -1, 41 |--> 1, 37 |--> 1, Dirichlet character modulo 60 of conductor 3 mapping 31 |--> 1, 41 |--> -1, 37 |--> 1, Dirichlet character modulo 60 of conductor 5 mapping 31 |--> 1, 41 |--> 1, 37 |--> zeta4)
        sage: val = G.gens()[2].values_on_gens()[2] ; val
        zeta4
        sage: parent(val)
        Maximal Order in Cyclotomic Field of order 4 and degree 2
        sage: r4.residue_field(r4.ideal(29).factor()[0][0])(val)
        17
        sage: r4.residue_field(r4.ideal(29).factor()[0][0])(val) * GF(29)(3)
        22
        sage: r4.residue_field(r4.ideal(29).factor()[0][0])(G.gens()[2].values_on_gens()[2]) * 3
        22
        sage: parent(r4.residue_field(r4.ideal(29).factor()[0][0])(G.gens()[2].values_on_gens()[2]) * 3)
        Residue field of Fractional ideal (-2*zeta4 + 5)

    ::

        sage: DirichletGroup(60, integral=True)
        Group of Dirichlet characters of modulus 60 over Maximal Order in Cyclotomic Field of order 4 and degree 2
        sage: parent(DirichletGroup(60, integral=True).gens()[2].values_on_gens()[2])
        Maximal Order in Cyclotomic Field of order 4 and degree 2

    TESTS:

    Dirichlet groups are cached, creating two groups with the same parameters
    yields the same object::

        sage: DirichletGroup(60) is DirichletGroup(60)
        True

    """
    def create_key(self, N, base_ring=None, zeta=None, zeta_order=None, names=None, integral=False):
        """
        Create a key that uniquely determines a Dirichlet group.

        TESTS::

            sage: DirichletGroup.create_key(60)
            (Cyclotomic Field of order 4 and degree 2, 60, zeta4, 4)

        An example to illustrate that ``base_ring`` is a part of the key::

            sage: k = DirichletGroup.create_key(2, base_ring=QQ); k
            (Rational Field, 2, 1, 1)
            sage: l = DirichletGroup.create_key(2, base_ring=CC); l
            (Complex Field with 53 bits of precision, 2, 1.00000000000000, 1)
            sage: k == l
            False
            sage: G = DirichletGroup.create_object(None, k); G
            Group of Dirichlet characters of modulus 2 over Rational Field
            sage: H = DirichletGroup.create_object(None, l); H
            Group of Dirichlet characters of modulus 2 over Complex Field with 53 bits of precision
            sage: G == H
            False

        If ``base_ring`` was not be a part of the key, the keys would compare
        equal and the caching would be broken::

            sage: k = k[1:]; k
            (2, 1, 1)
            sage: l = l[1:]; l
            (2, 1.00000000000000, 1)
            sage: k == l
            True
            sage: DirichletGroup(2, base_ring=QQ) is DirichletGroup(2, base_ring=CC)
            False

        """
        modulus = rings.Integer(N)

        if base_ring is None:
            if not (zeta is None and zeta_order is None):
                raise ValueError, "zeta and zeta_order must be None if base_ring not specified."
            e = rings.IntegerModRing(modulus).unit_group_exponent()
            base_ring = rings.CyclotomicField(e)
            if integral:
                base_ring = base_ring.ring_of_integers()

        if not is_Ring(base_ring):
            raise TypeError, "base_ring (=%s) must be a ring"%base_ring

        if zeta is None:
            e = rings.IntegerModRing(modulus).unit_group_exponent()
            for d in reversed(e.divisors()):
                try:
                    zeta = base_ring.zeta(d)
                    zeta_order = d
                    break
                except ValueError:
                    pass
        elif zeta_order is None:
            zeta_order = zeta.multiplicative_order()

        return (base_ring, modulus, zeta, zeta_order)

    def create_object(self, version, key, **extra_args):
        """
        Create the object from the key (extra arguments are ignored). This is
        only called if the object was not found in the cache.

        TESTS::

            sage: K = CyclotomicField(4)
            sage: DirichletGroup.create_object(None, (K, 60, K.gen(), 4))
            Group of Dirichlet characters of modulus 60 over Cyclotomic Field of order 4 and degree 2

        """
        base_ring, modulus, zeta, zeta_order = key
        return DirichletGroup_class(modulus, zeta, zeta_order)

DirichletGroup = DirichletGroupFactory("DirichletGroup")

def is_DirichletGroup(x):
    """
    Returns True if x is a Dirichlet group.

    EXAMPLES::

        sage: from sage.modular.dirichlet import is_DirichletGroup
        sage: is_DirichletGroup(DirichletGroup(11))
        True
        sage: is_DirichletGroup(11)
        False
        sage: is_DirichletGroup(DirichletGroup(11).0)
        False
    """
    return isinstance(x, DirichletGroup_class)

class DirichletGroup_class(parent_gens.ParentWithMultiplicativeAbelianGens):
    """
    Group of Dirichlet characters modulo `N` over a given base
    ring `R`.
    """
    def __init__(self, modulus, zeta, zeta_order):
        r"""
        Create a Dirichlet group. Not to be called directly (use the factory function ``DirichletGroup``.)

        EXAMPLES::

            sage: G = DirichletGroup(7, base_ring = Integers(9), zeta = Integers(9)(2)) # indirect doctest
            sage: G.base() # check that ParentWithBase.__init__ has been called
            Ring of integers modulo 9
        """
        parent_gens.ParentWithMultiplicativeAbelianGens.__init__(self, zeta.parent())
        self._zeta = zeta
        self._zeta_order = int(zeta_order)
        self._modulus = modulus
        self._integers = rings.IntegerModRing(modulus)
        a = zeta.parent()(1)
        v = {a:0}
        w = [a]
        if is_ComplexField(zeta.parent()):
            for i in range(1, self._zeta_order):
                a = a * zeta
                a._set_multiplicative_order(zeta_order/arith.GCD(zeta_order, i))
                v[a] = i
                w.append(a)
        else:
            for i in range(1, self._zeta_order):
                a = a * zeta
                v[a] = i
                w.append(a)
        self._zeta_powers = w  # gives quickly the ith power of zeta
        self._zeta_dlog = v    # dictionary that computes log_{zeta}(power of zeta).
        self._module = free_module.FreeModule(rings.IntegerModRing(zeta_order),
                                              len(self._integers.unit_gens()))

    def change_ring(self, R, zeta=None, zeta_order=None):
        """
        Returns the Dirichlet group over R with the same modulus as self.

        EXAMPLES::

            sage: G = DirichletGroup(7,QQ); G
            Group of Dirichlet characters of modulus 7 over Rational Field
            sage: G.change_ring(CyclotomicField(6))
            Group of Dirichlet characters of modulus 7 over Cyclotomic Field of order 6 and degree 2
        """
        return DirichletGroup(self.modulus(), R,
                              zeta=zeta,
                              zeta_order=zeta_order)

    def base_extend(self, R):
        """
        Returns the Dirichlet group over R obtained by extending scalars, with the same modulus and root of unity as self.

        EXAMPLES::

            sage: G = DirichletGroup(7,QQ); G
            Group of Dirichlet characters of modulus 7 over Rational Field
            sage: H = G.base_extend(CyclotomicField(6)); H
            Group of Dirichlet characters of modulus 7 over Cyclotomic Field of order 6 and degree 2
            sage: H.zeta()
            -1
            sage: G.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: No coercion map from 'Rational Field' to 'Integer Ring' is defined.

        """
        if not R.has_coerce_map_from(self.base_ring()):
            raise TypeError("No coercion map from '%s' to '%s' is defined." % (self.base_ring(), R))
        return DirichletGroup(self.modulus(), R, zeta=R(self.zeta()), zeta_order=self.zeta_order())

    def __call__(self, x):
        """
        Coerce x into this Dirichlet group.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: K = G.base_ring()
            sage: G(1)
            Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1
            sage: G([-1])
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1
            sage: G([K.0])
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12
            sage: G(0)
            Traceback (most recent call last):
            ...
            TypeError: No coercion of 0 into Group of Dirichlet characters of modulus 13 over Cyclotomic Field of order 12 and degree 4 defined.
        """
        if isinstance(x, (int,rings.Integer)) and x==1:
            R = self.base_ring()
            x = [R(1) for _ in self.unit_gens()]
        if isinstance(x, list):  # list of values on each unit generator
            return DirichletCharacter(self, x)
        elif isinstance(x, DirichletCharacter):  # coercion
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return DirichletCharacter(self, x.values_on_gens())
            return self._coerce_in_dirichlet_character(x)
        raise TypeError("No coercion of %s into %s defined."%(x, self))

    def _coerce_in_dirichlet_character(self, x):
        r"""
        Create an element of this group from a Dirichlet character, whose
        conductor must divide the modulus of self.

        EXAMPLES::

            sage: G = DirichletGroup(6)
            sage: G._coerce_in_dirichlet_character(DirichletGroup(3).0)
            Dirichlet character modulo 6 of conductor 3 mapping 5 |--> -1
            sage: G._coerce_in_dirichlet_character(DirichletGroup(15).0)
            Dirichlet character modulo 6 of conductor 3 mapping 5 |--> -1
            sage: G._coerce_in_dirichlet_character(DirichletGroup(15).1)
            Traceback (most recent call last):
            ...
            TypeError: conductor must divide modulus
            sage: H = DirichletGroup(16, QQ); H._coerce_in_dirichlet_character(DirichletGroup(16).1)
            Traceback (most recent call last):
            ...
            ValueError: cannot coerce element of order 4 into self
        """

        if self.modulus() % x.conductor() != 0:
            raise TypeError("conductor must divide modulus")
        elif not x.order().divides(self._zeta_order):
            raise ValueError("cannot coerce element of order %s into self"%x.order())
        a = []
        R = self.base_ring()
        for u in self.unit_gens():
            v = u.lift()
            # have to do this, since e.g., unit gens mod 11 are not units mod 22.
            while x.modulus().gcd(v) != 1:
                v += self.modulus()
            a.append(R(x(v)))
        return self(a)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of x into self.

        Note that although there is conversion between Dirichlet groups of
        different moduli, there is never canonical coercion. The main result of
        this is that Dirichlet characters of different moduli never compare as
        equal.

        TESTS::

            sage: trivial_character(6) == trivial_character(3) # indirect doctest
            False
            sage: trivial_character(3) == trivial_character(9)
            False
            sage: trivial_character(3) == DirichletGroup(3, QQ).0^2
            True
        """
        if isinstance(x, DirichletCharacter) and self.modulus() == x.modulus() and \
               self.base_ring().has_coerce_map_from(x.base_ring()) and \
               self.zeta_order() % x.parent().zeta_order() == 0:
            return self._coerce_in_dirichlet_character(x)
        raise TypeError

    def __cmp__(self, other):
        """
        Compare two Dirichlet groups. They are equal if they have the same
        modulus, are over the same base ring, and have the same chosen root
        of unity. Otherwise we compare first on the modulus, then the base
        ring, and finally the root of unity.

        EXAMPLES::

            sage: DirichletGroup(13) == DirichletGroup(13)
            True
            sage: DirichletGroup(13) == DirichletGroup(13,QQ)
            False
            sage: DirichletGroup(11) < DirichletGroup(13,QQ)
            True
            sage: DirichletGroup(17) > DirichletGroup(13)
            True
        """
        if not is_DirichletGroup(other):
            return -1
        c = cmp(self.modulus(), other.modulus())
        if c:
            return c
        c = cmp(self.base_ring(), other.base_ring())
        if c:
            return c
        c = cmp(self._zeta, other._zeta)
        if c:
            return c
        return 0

    def __len__(self):
        """
        Return the number of elements of this Dirichlet group. This is the
        same as self.order().

        EXAMPLES::

            sage: len(DirichletGroup(20))
            8
            sage: len(DirichletGroup(20, QQ))
            4
            sage: len(DirichletGroup(20, GF(5)))
            8
            sage: len(DirichletGroup(20, GF(2)))
            1
            sage: len(DirichletGroup(20, GF(3)))
            4
        """
        return self.order()

    def _repr_(self):
        """
        Return a print representation of this group, which can be renamed.

        EXAMPLES::

            sage: G = DirichletGroup(11)
            sage: repr(G) # indirect doctest
            'Group of Dirichlet characters of modulus 11 over Cyclotomic Field of order 10 and degree 4'
            sage: G.rename('Dir(11)')
            sage: G
            Dir(11)
        """
        return "Group of Dirichlet characters of modulus %s over %s"%\
               (self.modulus(),self.base_ring())

# Commented out base_ring since this class should inherit base_ring from superclass ParentWithBase.
#   -- David Loeffler, 2009-04-09.
#
#    def base_ring(self):
#        """
#        Returns the base ring of self.
#
#        EXAMPLES::
#
#            sage: DirichletGroup(11).base_ring()
#            Cyclotomic Field of order 10 and degree 4
#            sage: DirichletGroup(11,QQ).base_ring()
#            Rational Field
#            sage: DirichletGroup(11,GF(7)).base_ring()
#            Finite Field of size 7
#            sage: DirichletGroup(20).base_ring()
#            Cyclotomic Field of order 4 and degree 2
#        """
#        return self._zeta.parent()

    @cached_method
    def decomposition(self):
        r"""
        Returns the Dirichlet groups of prime power modulus corresponding
        to primes dividing modulus.

        (Note that if the modulus is 2 mod 4, there will be a "factor" of
        `(\ZZ/2\ZZ)^*`, which is the trivial group.)

        EXAMPLES::

            sage: DirichletGroup(20).decomposition()
            [
            Group of Dirichlet characters of modulus 4 over Cyclotomic Field of order 4 and degree 2,
            Group of Dirichlet characters of modulus 5 over Cyclotomic Field of order 4 and degree 2
            ]
            sage: DirichletGroup(20,GF(5)).decomposition()
            [
            Group of Dirichlet characters of modulus 4 over Finite Field of size 5,
            Group of Dirichlet characters of modulus 5 over Finite Field of size 5
            ]
        """
        R = self.base_ring()
        return Sequence([DirichletGroup(p**r,R) for p, r \
                           in arith.factor(self.modulus())],
                                cr=True,
                                universe = cat.Objects())

    def exponent(self):
        """
        Return the exponent of this group.

        EXAMPLES::

            sage: DirichletGroup(20).exponent()
            4
            sage: DirichletGroup(20,GF(3)).exponent()
            2
            sage: DirichletGroup(20,GF(2)).exponent()
            1
            sage: DirichletGroup(37).exponent()
            36
        """
        return self._zeta_order

    @cached_method
    def _automorphisms(self):
        """
        Compute the automorphisms of self. These are always given by raising to
        a power, so the return value is a list of integers.

        At present this is only implemented if the base ring has characteristic 0 or a prime.

        EXAMPLES::

            sage: DirichletGroup(17)._automorphisms()
            [1, 3, 5, 7, 9, 11, 13, 15]
            sage: DirichletGroup(17, GF(11^4, 'a'))._automorphisms()
            [1, 11, 121, 1331]
            sage: DirichletGroup(17, Integers(6), zeta=Integers(6)(5))._automorphisms()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automorphisms for finite non-field base rings not implemented
            sage: DirichletGroup(17, Integers(9), zeta=Integers(9)(2))._automorphisms()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automorphisms for finite non-field base rings not implemented
        """
        n = self.zeta_order()
        R = self.base_ring()
        p = R.characteristic()
        if p == 0:
            Auts = [e for e in xrange(1,n) if arith.GCD(e,n) == 1]
        else:
            if not rings.ZZ(p).is_prime():
                raise NotImplementedError("Automorphisms for finite non-field base rings not implemented")
            # The automorphisms in characteristic p are
            # k-th powering for
            #         k = 1, p, p^2, ..., p^(r-1),
            # where p^r = 1 (mod n), so r is the mult order of p modulo n.
            r = rings.IntegerModRing(n)(p).multiplicative_order()
            Auts = [p**m for m in xrange(0,r)]
        return Auts

    def galois_orbits(self, v=None, reps_only=False, sort=True, check=True):
        """
        Return a list of the Galois orbits of Dirichlet characters in self,
        or in v if v is not None.

        INPUT:


        -  ``v`` - (optional) list of elements of self

        -  ``reps_only`` - (optional: default False) if True
           only returns representatives for the orbits.

        -  ``sort`` - (optional: default True) whether to sort
           the list of orbits and the orbits themselves (slightly faster if
           False).

        -  ``check`` - (optional, default: True) whether or not
           to explicitly coerce each element of v into self.


        The Galois group is the absolute Galois group of the prime subfield
        of Frac(R). If R is not a domain, an error will be raised.

        EXAMPLES::

            sage: DirichletGroup(20).galois_orbits()
            [
            [Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1], ... [Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -1]
            ]
            sage: DirichletGroup(17, Integers(6), zeta=Integers(6)(5)).galois_orbits()
            Traceback (most recent call last):
            ...
            TypeError: Galois orbits only defined if base ring is an integral domain
            sage: DirichletGroup(17, Integers(9), zeta=Integers(9)(2)).galois_orbits()
            Traceback (most recent call last):
            ...
            TypeError: Galois orbits only defined if base ring is an integral domain
        """

        if v is None:
            v = self.list()
        else:
            if check:
                v = [self(x) for x in v]

        G = []
        n = self.zeta_order()
        R = self.base_ring()
        p = R.characteristic()
        seen_so_far = set([])
        for x in v:
            z = x.element()
            e = tuple(z)   # change when there are immutable vectors (and below)
            if e in seen_so_far:
                continue
            orbit = x.galois_orbit(sort=sort)
            if reps_only:
                G.append(x)
            else:
                G.append(orbit)
            for z in orbit:
                seen_so_far.add(tuple(z.element()))
        G = Sequence(G, cr=True)
        if sort:
            G.sort()
        return G

    def gen(self, n=0):
        """
        Return the n-th generator of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.gen(0)
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: G.gen(1)
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4
            sage: G.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: n(=2) must be between 0 and 1

        ::

            sage: G.gen(-1)
            Traceback (most recent call last):
            ...
            IndexError: n(=-1) must be between 0 and 1
        """
        n = int(n)
        g = self.gens()
        if n<0 or n>=len(g):
            raise IndexError("n(=%s) must be between 0 and %s"%(n,len(g)-1))
        return g[n]

    def gens(self):
        """
        Returns generators of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.gens()
            (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1, Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4)
        """
        if not (self._gens is None):
            return self._gens
        self._gens = []
        ug = self.unit_gens()
        R = self.base_ring()
        one = [R(1) for i in range(len(ug))]
        zeta = self.zeta()
        ord = self.zeta_order()
        M = self._module
        zero = M(0)
        orders = self.integers_mod().unit_group().gens_orders()
        for i in range(len(ug)):
            z = zero.__copy__()
            z[i] = ord//arith.GCD(ord, orders[i])
            self._gens.append(DirichletCharacter(self, z, check=False))
        self._gens = tuple(self._gens)
        return self._gens

    def integers_mod(self):
        r"""
        Returns the group of integers `\ZZ/N\ZZ`
        where `N` is the modulus of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.integers_mod()
            Ring of integers modulo 20
        """
        return self._integers

    def modulus(self):
        """
        Returns the modulus of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.modulus()
            20
        """
        return self._modulus

    def ngens(self):
        """
        Returns the number of generators of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.ngens()
            2
        """
        return len(self.gens())

    @cached_method
    def order(self):
        """
        Return the number of elements of self. This is the same as
        len(self).

        EXAMPLES::

            sage: DirichletGroup(20).order()
            8
            sage: DirichletGroup(37).order()
            36
        """
        ord = rings.Integer(1)
        for g in self.gens():
            ord *= int(g.order())
        return ord

    def random_element(self):
        """
        Return a random element of self.

        The element is computed by multiplying a random power of each
        generator together, where the power is between 0 and the order of
        the generator minus 1, inclusive.

        EXAMPLES::

            sage: DirichletGroup(37).random_element()
            Dirichlet character modulo 37 of conductor 37 mapping 2 |--> zeta36^4
            sage: DirichletGroup(20).random_element()
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: DirichletGroup(60).random_element()
            Dirichlet character modulo 60 of conductor 3 mapping 31 |--> 1, 41 |--> -1, 37 |--> 1
        """
        e = self(1)
        for i in range(self.ngens()):
            g = self.gen(i)
            n = random.randrange(g.order())
            e *= g**n
        return e

    def unit_gens(self):
        r"""
        Returns the minimal generators for the units of
        `(\ZZ/N\ZZ)^*`, where `N` is the
        modulus of self.

        EXAMPLES::

            sage: DirichletGroup(37).unit_gens()
            (2,)
            sage: DirichletGroup(20).unit_gens()
            (11, 17)
            sage: DirichletGroup(60).unit_gens()
            (31, 41, 37)
            sage: DirichletGroup(20,QQ).unit_gens()
            (11, 17)
        """
        return self._integers.unit_gens()

    def zeta(self):
        """
        Returns the chosen root zeta of unity in the base ring
        `R`.

        EXAMPLES::

            sage: DirichletGroup(37).zeta()
            zeta36
            sage: DirichletGroup(20).zeta()
            zeta4
            sage: DirichletGroup(60).zeta()
            zeta4
            sage: DirichletGroup(60,QQ).zeta()
            -1
            sage: DirichletGroup(60, GF(25,'a')).zeta()
            2
        """
        return self._zeta

    def zeta_order(self):
        """
        Returns the order of the chosen root zeta of unity in the base ring
        `R`.

        EXAMPLES::

            sage: DirichletGroup(20).zeta_order()
            4
            sage: DirichletGroup(60).zeta_order()
            4
            sage: DirichletGroup(60, GF(25,'a')).zeta_order()
            4
            sage: DirichletGroup(19).zeta_order()
            18
        """
        return self._zeta_order




