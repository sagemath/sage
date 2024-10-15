# sage_setup: distribution = sagemath-objects
"""
Functors

AUTHORS:

- David Kohel and William Stein

- David Joyner (2005-12-17): examples

- Robert Bradshaw (2007-06-23): Pyrexify

- Simon King (2010-04-30): more examples, several bug fixes,
  re-implementation of the default call method,
  making functors applicable to morphisms (not only to objects)

- Simon King (2010-12): Pickling of functors without losing domain and codomain
"""

# ****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories import category


def _Functor_unpickle(Cl, D, domain, codomain):
    """
    Generic unpickling function for functors.

    AUTHOR:

    - Simon King (2010-12): :issue:`10460`

    EXAMPLES::

        sage: R.<x,y> = InfinitePolynomialRing(QQ)
        sage: F = R.construction()[0]
        sage: F == loads(dumps(F))
        True
        sage: F.domain(), loads(dumps(F)).domain()
        (Category of rings, Category of rings)
    """
    F = Functor.__new__(Cl)
    Functor.__init__(F, domain, codomain)
    for s, v in D:
        setattr(F, s, v)
    return F


cdef class Functor(SageObject):
    """
    A class for functors between two categories.

    NOTE:

    - In the first place, a functor is given by its domain and codomain,
      which are both categories.
    - When defining a sub-class, the user should not implement a call method.
      Instead, one should implement three methods, which are composed in the
      default call method:

      - ``_coerce_into_domain(self, x)`` -- return an object of ``self``'s
        domain, corresponding to ``x``, or raise a :exc:`TypeError`.
        - Default: Raise :exc:`TypeError` if ``x`` is not in ``self``'s domain.

      - ``_apply_functor(self, x)`` -- apply ``self`` to an object ``x`` of
        ``self``'s domain.
        - Default: Conversion into ``self``'s codomain.

      - ``_apply_functor_to_morphism(self, f)`` -- apply ``self`` to a morphism
        ``f`` in ``self``'s domain.
        - Default: Return ``self(f.domain()).hom(f,self(f.codomain()))``.

    EXAMPLES::

        sage: rings  = Rings()
        sage: abgrps = CommutativeAdditiveGroups()
        sage: F = ForgetfulFunctor(rings, abgrps)
        sage: F.domain()
        Category of rings
        sage: F.codomain()
        Category of commutative additive groups
        sage: from sage.categories.functor import Functor
        sage: isinstance(F, Functor)
        True
        sage: I = IdentityFunctor(abgrps)
        sage: I
        The identity functor on Category of commutative additive groups
        sage: I.domain()
        Category of commutative additive groups
        sage: isinstance(I, Functor)
        True

    Note that by default, an instance of the class Functor is coercion
    from the domain into the codomain. The above subclasses overloaded
    this behaviour. Here we illustrate the default::

        sage: from sage.categories.functor import Functor
        sage: F = Functor(Rings(), Fields())
        sage: F
        Functor from Category of rings to Category of fields
        sage: F(ZZ)
        Rational Field
        sage: F(GF(2))
        Finite Field of size 2

    Functors are not only about the objects of a category, but also about
    their morphisms. We illustrate it, again, with the coercion functor
    from rings to fields.

    ::

        sage: R1.<x> = ZZ[]
        sage: R2.<a,b> = QQ[]
        sage: f = R1.hom([a + b], R2)
        sage: f
        Ring morphism:
          From: Univariate Polynomial Ring in x over Integer Ring
          To:   Multivariate Polynomial Ring in a, b over Rational Field
          Defn: x |--> a + b
        sage: F(f)                                                                      # needs sage.libs.singular
        Ring morphism:
          From: Fraction Field of Univariate Polynomial Ring in x over Integer Ring
          To:   Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
          Defn: x |--> a + b
        sage: F(f)(1/x)                                                                 # needs sage.libs.singular
        1/(a + b)

    We can also apply a polynomial ring construction functor to our homomorphism. The
    result is a homomorphism that is defined on the base ring::

        sage: F = QQ['t'].construction()[0]
        sage: F
        Poly[t]
        sage: F(f)
        Ring morphism:
          From: Univariate Polynomial Ring in t
                 over Univariate Polynomial Ring in x over Integer Ring
          To:   Univariate Polynomial Ring in t
                 over Multivariate Polynomial Ring in a, b over Rational Field
          Defn: Induced from base ring by
                Ring morphism:
                  From: Univariate Polynomial Ring in x over Integer Ring
                  To:   Multivariate Polynomial Ring in a, b over Rational Field
                  Defn: x |--> a + b
        sage: p = R1['t']('(-x^2 + x)*t^2 + (x^2 - x)*t - 4*x^2 - x + 1')
        sage: F(f)(p)
        (-a^2 - 2*a*b - b^2 + a + b)*t^2 + (a^2 + 2*a*b + b^2 - a - b)*t
        - 4*a^2 - 8*a*b - 4*b^2 - a - b + 1
    """
    def __init__(self, domain, codomain):
        """
        TESTS::

            sage: from sage.categories.functor import Functor
            sage: F = Functor(Rings(), Fields())
            sage: F
            Functor from Category of rings to Category of fields
            sage: F(ZZ)
            Rational Field
            sage: F(GF(2))
            Finite Field of size 2
        """
        if not isinstance(domain, category.Category):
            raise TypeError("domain (=%s) must be a category" % domain)
        if not isinstance(codomain, category.Category):
            raise TypeError("codomain (=%s) must be a category" % codomain)
        self.__domain = domain
        self.__codomain = codomain

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.categories.functor import Functor
            sage: F = Functor(Rings(), Fields())
            sage: hash(F)  # random
            42
        """
        return hash(self.__domain) ^ hash(self.__codomain)

    def __reduce__(self):
        """
        Generic pickling of functors.

        AUTHOR:

        - Simon King (2010-12):  :issue:`10460`

        TESTS::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: F == loads(dumps(F))
            True
            sage: F.codomain()
            Category of rings
        """
        return (_Functor_unpickle,
                (self.__class__, list(self.__dict__.items()),
                 self.__domain, self.__codomain))

    def _apply_functor(self, x):
        """
        Apply the functor to an object of ``self``'s domain.

        NOTE:

        Each subclass of :class:`Functor` should overload this method. By default,
        this method coerces into the codomain, without checking whether the
        argument belongs to the domain.

        TESTS::

            sage: from sage.categories.functor import Functor
            sage: F = Functor(FiniteFields(), Fields())
            sage: F._apply_functor(ZZ)
            Rational Field
        """
        return self.__codomain(x)

    def _apply_functor_to_morphism(self, f):
        """
        Apply the functor to a morphism between two objects of ``self``'s domain.

        NOTE:

        Each subclass of :class:`Functor` should overload this method. By
        default, this method coerces into the codomain, without checking
        whether the argument belongs to the domain.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.categories.functor import Functor
            sage: F = Functor(Rings(), Fields())
            sage: k.<a> = GF(25)
            sage: f = k.hom([-a - 4])
            sage: R.<t> = k[]
            sage: fR = R.hom(f, R)
            sage: fF = F(fR)         # indirect doctest
            sage: fF
            Ring endomorphism of Fraction Field of
             Univariate Polynomial Ring in t over Finite Field in a of size 5^2
              Defn: Induced from base ring by
                    Ring endomorphism of Univariate Polynomial Ring in t
                     over Finite Field in a of size 5^2
                      Defn: Induced from base ring by
                            Ring endomorphism of Finite Field in a of size 5^2
                              Defn: a |--> 4*a + 1
            sage: fF((a^2+a)*t^2/(a*t - a^2))
            ((4*a + 2)*t^2)/(t + a + 4)
        """
        try:
            return self(f.domain()).hom(f, self(f.codomain()))
        except Exception:
            raise TypeError(f'unable to transform {f} into a morphism in {self.codomain()}')

    def _coerce_into_domain(self, x):
        """
        Interpret the argument as an object of ``self``'s domain.

        NOTE:

        A subclass of :class:`Functor` may overload this method. It should
        return an object of ``self``'s domain, and should raise a
        :exc:`TypeError` if this is impossible.

        By default, the argument will not be changed, but a
        :exc:`TypeError` will be raised if the argument does not
        belong to the domain.

        TESTS::

            sage: from sage.categories.functor import Functor
            sage: F = Functor(Fields(),Fields())
            sage: F(QQ)
            Rational Field
            sage: F(ZZ) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: x (=Integer Ring) is not in Category of fields
        """
        if x not in self.__domain:
            raise TypeError("x (=%s) is not in %s" % (x, self.__domain))
        return x

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.functor import Functor
            sage: F = Functor(Rings(),Fields())
            sage: F #indirect doctest
            Functor from Category of rings to Category of fields

        A functor can be renamed if its type is a Python class
        (see :issue:`16156`)::

            sage: I = IdentityFunctor(Rings()); I
            The identity functor on Category of rings
            sage: I.rename('Id'); I
            Id
        """
        return "Functor from %s to %s" % (self.__domain, self.__codomain)

    def __call__(self, x):
        """
        NOTE:

        Implement _coerce_into_domain, _apply_functor and
        _apply_functor_to_morphism when subclassing Functor.

        TESTS:

        The default::

            sage: from sage.categories.functor import Functor
            sage: F = Functor(Rings(),Fields())
            sage: F
            Functor from Category of rings to Category of fields
            sage: F(ZZ)
            Rational Field
            sage: F(GF(2))
            Finite Field of size 2

        Two subclasses::

            sage: F1 = ForgetfulFunctor(FiniteFields(), Fields())
            sage: F1(GF(5))  # indirect doctest
            Finite Field of size 5
            sage: F1(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: x (=Integer Ring) is not in Category of finite enumerated fields
            sage: F2 = IdentityFunctor(Fields())
            sage: F2(RR) is RR #indirect doctest
            True
            sage: F2(ZZ['x','y'])
            Traceback (most recent call last):
            ...
            TypeError: x (=Multivariate Polynomial Ring in x, y over Integer Ring)
            is not in Category of fields

        The last example shows that it is tested whether the result of
        applying the functor lies in the functor's codomain. Note that
        the matrix functor used to be defined similar to this example,
        which was fixed in :issue:`8807`::

            sage: class IllFunctor(Functor):
            ....:   def __init__(self, m, n):
            ....:       self._m = m
            ....:       self._n = n
            ....:       Functor.__init__(self, Rings(), Rings())
            ....:   def _apply_functor(self, R):
            ....:       return MatrixSpace(R, self._m, self._n)
            sage: F = IllFunctor(2, 2)
            sage: F(QQ)                                                                 # needs sage.modules
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: F = IllFunctor(2, 3)
            sage: F(QQ)                                                                 # needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: Functor from Category of rings to Category of rings
            is ill-defined, since it sends x (=Rational Field) to something
            that is not in Category of rings.
        """
        from sage.categories.morphism import Morphism
        if isinstance(x, Morphism):
            return self._apply_functor_to_morphism(x)
        y = self._apply_functor(self._coerce_into_domain(x))
        if not ((y in self.__codomain) or (y in self.__codomain.Homsets())):
            raise TypeError("%s is ill-defined, since it sends x (=%s) to something that is not in %s." % (repr(self), x, self.__codomain))
        return y

    def domain(self):
        """
        The domain of ``self``.

        EXAMPLES::

            sage: F = ForgetfulFunctor(FiniteFields(), Fields())
            sage: F.domain()
            Category of finite enumerated fields
        """
        return self.__domain

    def codomain(self):
        """
        The codomain of ``self``.

        EXAMPLES::

            sage: F = ForgetfulFunctor(FiniteFields(), Fields())
            sage: F.codomain()
            Category of fields
        """
        return self.__codomain


def is_Functor(x):
    """
    Test whether the argument is a functor.

    This function is deprecated.

    EXAMPLES::

        sage: from sage.categories.functor import is_Functor
        sage: F1 = QQ.construction()[0]
        sage: F1
        FractionField
        sage: is_Functor(F1)
        doctest:warning...
        DeprecationWarning: The function is_Functor is deprecated;
        use 'isinstance(..., Functor)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        True
        sage: is_Functor(FractionField)
        False
        sage: F2 = ForgetfulFunctor(Fields(), Rings())
        sage: F2
        The forgetful functor from Category of fields to Category of rings
        sage: is_Functor(F2)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_Functor is deprecated; "
                "use 'isinstance(..., Functor)' instead.")
    return isinstance(x, Functor)


###########################################
# The natural functors in Sage
###########################################

class ForgetfulFunctor_generic(Functor):
    """
    The forgetful functor, i.e., embedding of a subcategory.

    NOTE:

    Forgetful functors should be created using :func:`ForgetfulFunctor`,
    since the init method of this class does not check whether the
    domain is a subcategory of the codomain.

    EXAMPLES::

        sage: F = ForgetfulFunctor(FiniteFields(), Fields())  # indirect doctest
        sage: F
        The forgetful functor
         from Category of finite enumerated fields
           to Category of fields
        sage: F(GF(3))
        Finite Field of size 3
    """
    def __reduce__(self):
        """
        EXAMPLES::

            sage: F = ForgetfulFunctor(Groups(), Sets())
            sage: loads(F.dumps()) == F
            True
        """
        return ForgetfulFunctor, (self.domain(), self.codomain())

    def _repr_(self):
        """
        TESTS::

            sage: F = ForgetfulFunctor(FiniteFields(), Fields())
            sage: F  # indirect doctest
            The forgetful functor from Category of finite enumerated fields
             to Category of fields
        """
        return "The forgetful functor from %s to %s" % (self.domain(),
                                                        self.codomain())

    def __eq__(self, other):
        """
        NOTE:

        It is tested whether the second argument belongs to the class
        of forgetful functors and has the same domain and codomain as
        ``self``. If the second argument is a functor of a different class
        but happens to be a forgetful functor, both arguments will
        still be considered as being *different*.

        TESTS::

            sage: F1 = ForgetfulFunctor(FiniteFields(), Fields())

        This is to test against a bug occurring in a previous version
        (see :issue:`8800`)::

            sage: F1 == QQ #indirect doctest
            False

        We now compare with the fraction field functor, that has a
        different domain::

            sage: F2 = QQ.construction()[0]
            sage: F1 == F2 #indirect doctest
            False
        """
        from sage.categories.pushout import IdentityConstructionFunctor
        if not isinstance(other, (self.__class__, IdentityConstructionFunctor)):
            return False
        return (self.domain() == other.domain() and
                self.codomain() == other.codomain())

    def __ne__(self, other):
        """
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: F1 = ForgetfulFunctor(FiniteFields(), Fields())
            sage: F1 != F1
            False
            sage: F1 != QQ
            True
        """
        return not self == other


class IdentityFunctor_generic(ForgetfulFunctor_generic):
    """
    Generic identity functor on any category.

    NOTE:

    This usually is created using :func:`IdentityFunctor`.

    EXAMPLES::

        sage: F = IdentityFunctor(Fields()) #indirect doctest
        sage: F
        The identity functor on Category of fields
        sage: F(RR) is RR
        True
        sage: F(ZZ)
        Traceback (most recent call last):
        ...
        TypeError: x (=Integer Ring) is not in Category of fields

    TESTS::

        sage: R = IdentityFunctor(Rings())
        sage: P, _ = QQ['t'].construction()
        sage: R == P
        False
        sage: P == R
        False
        sage: R == QQ
        False
    """
    def __init__(self, C):
        """
        TESTS::

            sage: from sage.categories.functor import IdentityFunctor_generic
            sage: F = IdentityFunctor_generic(Groups())
            sage: F == IdentityFunctor(Groups())
            True
            sage: F
            The identity functor on Category of groups
        """
        ForgetfulFunctor_generic.__init__(self, C, C)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: F = IdentityFunctor(Groups())
            sage: loads(F.dumps()) == F
            True
        """
        return IdentityFunctor, (self.domain(), )

    def _repr_(self):
        """
        TESTS::

            sage: fields = Fields()
            sage: F = IdentityFunctor(fields)
            sage: F #indirect doctest
            The identity functor on Category of fields
        """
        return "The identity functor on %s" % (self.domain())

    def _apply_functor(self, x):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: fields = Fields()
            sage: F = IdentityFunctor(fields)
            sage: F._apply_functor(QQ)
            Rational Field

        It is not tested here whether the argument belongs to the domain
        (this test is done in the default method ``_coerce_into_domain``)::

            sage: F._apply_functor(ZZ)
            Integer Ring
        """
        return x


def IdentityFunctor(C):
    """
    Construct the identity functor of the given category.

    INPUT:

    - ``C`` -- a category

    OUTPUT: the identity functor in ``C``

    EXAMPLES::

        sage: rings = Rings()
        sage: F = IdentityFunctor(rings)
        sage: F(ZZ['x','y']) is ZZ['x','y']
        True
    """
    return IdentityFunctor_generic(C)


def ForgetfulFunctor(domain, codomain):
    """
    Construct the forgetful function from one category to another.

    INPUT:

    - ``C``, ``D`` -- two categories

    OUTPUT:

    A functor that returns the corresponding object of ``D`` for
    any element of ``C``, by forgetting the extra structure.

    ASSUMPTION:

    The category ``C`` must be a sub-category of ``D``.

    EXAMPLES::

        sage: rings = Rings()
        sage: abgrps = CommutativeAdditiveGroups()
        sage: F = ForgetfulFunctor(rings, abgrps)
        sage: F
        The forgetful functor
         from Category of rings
           to Category of commutative additive groups

    It would be a mistake to call it in opposite order::

        sage: F = ForgetfulFunctor(abgrps, rings)
        Traceback (most recent call last):
        ...
        ValueError: Forgetful functor not supported for domain
        Category of commutative additive groups

    If both categories are equal, the forgetful functor is the same as the
    identity functor::

        sage: ForgetfulFunctor(abgrps, abgrps) == IdentityFunctor(abgrps)
        True
    """
    if domain == codomain:
        return IdentityFunctor(domain)
    if not domain.is_subcategory(codomain):
        raise ValueError("Forgetful functor not supported for domain %s" % domain)
    return ForgetfulFunctor_generic(domain, codomain)
