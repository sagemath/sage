r"""
Anderson motives

Let `\GF{q}[T]` be a polynomial ring with coefficients in a finite
field `\GF{q}` and let `K` be an extension of `\GF{q}` equipped
with a distinguished element `z`.

By definition, an Anderson motive attached to these data is a free
module of finite rank `M` over `K[T]`, equipped with a linear
automorphism

.. MATH::

    \tau_M : \tau^\star M \left[\frac 1{T-z}\right] \to M \left[\frac 1{T-z}\right]

where `\tau^\star M = K \otimes_{K, \text{Frob}} M` and the notation
means that `K` is viewed as an algebra over itself through the
Frobenius `\text{Frob} : x \mapsto x^q`.

.. RUBRIC:: Anderson motives attached to Drinfeld modules

Any Drinfeld module `\phi` over `(A, \gamma)` with `\gamma : A \to K,
T \mapsto z` gives rise to an Anderson motive. By definition, it is
`M(\phi) := K\{\tau\}` (the ring of Ore polynomials with commutation
rule `\tau \lambda = \lambda^q \tau` for `\lambda \in K`) where

- the structure of `\GF{q}[T]`-module is given by right multiplication
  by `\phi_a` (`a \in \GF{q}[T]`),

- the structure of `K`-module is given by left multiplication,

- the automorphism `\tau_{M(\phi)}` is the left multiplication
  by `\tau` in the Ore polynomial ring.

Anderson motives are nevertheless much more general than Drinfeld
modules. Besides, their linear nature allows for importing many
interesting construction of linear and bilinear algebra.

In SageMath, one can create the Anderson motive corresponding to
a Drinfeld module as follows::

    sage: k = GF(5)
    sage: A.<T> = k[]
    sage: K.<z> = k.extension(3)
    sage: phi = DrinfeldModule(A, [z, z^2, z^3, z^4])
    sage: M = phi.anderson_motive()
    sage: M
    Anderson motive of Drinfeld module defined by T |--> (2*z^2 + 2*z)*τ^3 + (2*z + 2)*τ^2 + z^2*τ + z

We see that `M` has rank `3`; it is actually a general fact that
the Anderson motive attached a Drinfeld module has the same rank
than the underlying Drinfeld module.

The canonical basis corresponds to the vectors `\tau^i` for `i`
varying between `0` and `r-1` where `r` is the rank::

    sage: tau = phi.ore_variable()
    sage: M(tau^0)
    (1, 0, 0)
    sage: M(tau^1)
    (0, 1, 0)
    sage: M(tau^2)
    (0, 0, 1)

Higher powers of `\tau` can be rewritten as linear combinations
(over `K[T]`!) of those three ones::

    sage: M(tau^3)
    ((z^2 + 3*z)*T + 2*z^2 + 3*z + 3, 3*z^2 + 2*z + 4, 2*z^2 + 1)
    sage: M(tau^4)
    ((4*z^2 + 4*z + 3)*T + z^2 + 4*z + 2, (z^2 + 4*z)*T + 3, 3*z^2 + 4*z + 4)

The matrix of the operator `\tau_M` can be obtained using the method
:meth:`matrix`::

    sage: M.matrix()
    [                              0                               1                               0]
    [                              0                               0                               1]
    [(z^2 + 3*z)*T + 2*z^2 + 3*z + 3                 3*z^2 + 2*z + 4                       2*z^2 + 1]

.. NOTE::

    Here, as it is conventional in SageMath, we use the row
    representation, meaning that the coordinates of the image
    by `\tau_M(\tau^i)` are written in the `i`-th row.

SageMath provides facilities to pick elements in `M` and perform
basic operations with them::

    sage: u, v, w = M.basis()
    sage: T*u + z*w
    (T, 0, z)
    sage: w.image()  # image by tau_M
    ((z^2 + 3*z)*T + 2*z^2 + 3*z + 3, 3*z^2 + 2*z + 4, 2*z^2 + 1)

It is also possible to give names to the vectors of the canonical
basis and then use when printing::

    sage: psi = DrinfeldModule(A, [z, z+1, z+2])
    sage: N.<e0, e1> = psi.anderson_motive()
    sage: N.random_element()  # random
    ((4*z+4)*T^2+(3*z^2+1)*T+z^2+3*z+3)*e0 + (T^2+(2*z^2+1)*T+3*z^2)*e1

.. RUBRIC:: More Anderson motives

One can also build the dual of the Anderson motive attached to
a Drinfeld simply by setting the attribute ``dual=True``::

    sage: Md = phi.anderson_motive(dual=True)
    sage: Md
    Dual Anderson motive of Drinfeld module defined by T |--> (2*z^2 + 2*z)*τ^3 + (2*z + 2)*τ^2 + z^2*τ + z
    sage: Md.matrix()
    [          z^2/(T + 4*z)                       1                       0]
    [    (2*z + 2)/(T + 4*z)                       0                       1]
    [(2*z^2 + 2*z)/(T + 4*z)                       0                       0]

We observe that some entries of the previous matrix have denominator
`T-z`. This corresponds to the fact that `\tau_M` is only defined
after inverting `T-z` in full generality, and it implies in particular
that ``Md`` does not come itself from a Drinfeld module.

Finally, SageMath also provides a general constructor :func:`AndersonMotive`
which allows in particular to explicitly provide the matrix of `\tau_M`::

    sage: mat = matrix(2, 2, [[T, z], [1, 1]])
    sage: N = AndersonMotive(A, mat)
    sage: N
    Anderson motive of rank 2 over Univariate Polynomial Ring in T over Finite Field in z of size 5^3
    sage: N.matrix()
    [T z]
    [1 1]

.. RUBRIC:: Morphisms between Anderson motives

By definition, a morphism between Anderson motives is a
`A \otimes K`-linear morphism commuting with the action of
`\tau`.

One important class of morphisms of Anderson motives are those
coming from isogenies between Drinfeld modules.
Such morphisms can be built easily as follows::

    sage: u = phi.hom(tau + z)
    sage: u
    Drinfeld Module morphism:
      From: Drinfeld module defined by T |--> (2*z^2 + 2*z)*τ^3 + (2*z + 2)*τ^2 + z^2*τ + z
      To:   Drinfeld module defined by T |--> (4*z^2 + 2*z + 4)*τ^3 + (4*z^2 + 1)*τ^2 + (z^2 + 2)*τ + z
      Defn: τ + z
    sage: Mu = u.anderson_motive()
    sage: Mu
    Morphism:
      From: Anderson motive of Drinfeld module defined by T |--> (4*z^2 + 2*z + 4)*τ^3 + (4*z^2 + 1)*τ^2 + (z^2 + 2)*τ + z
      To:   Anderson motive of Drinfeld module defined by T |--> (2*z^2 + 2*z)*τ^3 + (2*z + 2)*τ^2 + z^2*τ + z
    sage: Mu.matrix()
    [                              z                               1                               0]
    [                              0                 2*z^2 + 4*z + 4                               1]
    [(z^2 + 3*z)*T + 2*z^2 + 3*z + 3                 3*z^2 + 2*z + 4                               2]

Standard methods of linear algebra are available::

    sage: Mu.is_injective()
    True
    sage: Mu.is_surjective()
    False
    sage: Mu.image().basis()
    [(T + 3, 0, 0), (z, 1, 0), (z^2 + 2*z + 1, 0, 1)]

We check below that the characteristic polynomial of the Frobenius of
`\phi` is equal to the characteristic polynomial of the action of the
Frobenius on the motive::

    sage: f = phi.frobenius_endomorphism()
    sage: f
    Endomorphism of Drinfeld module defined by T |--> (2*z^2 + 2*z)*τ^3 + (2*z + 2)*τ^2 + z^2*τ + z
      Defn: τ^3
    sage: Mf = f.anderson_motive()
    sage: Mf.characteristic_polynomial()
    X^3 + (T + 4)*X^2 + 3*T^2*X + 4*T^3 + 2*T + 2

::

    sage: phi.frobenius_charpoly()
    X^3 + (T + 4)*X^2 + 3*T^2*X + 4*T^3 + 2*T + 2

AUTHOR:

- Xavier Caruso, Antoine Leudière (2025-11): initial version
"""

# *****************************************************************************
#        Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.latex import latex

from sage.categories.map import Map
from sage.categories.homset import Homset
from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.anderson_motives import AndersonMotives
from sage.structure.factorization import Factorization

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.ring import CommutativeRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.morphism import RingHomomorphism

from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix

from sage.modules.ore_module import OreModule, OreSubmodule, OreQuotientModule
from sage.modules.ore_module import normalize_names
from sage.modules.ore_module_homspace import OreModule_homspace
from sage.modules.ore_module_morphism import OreModuleMorphism

from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism


# Classes for Anderson motives
##############################

class AndersonMotive_general(OreModule):
    r"""
    General class for Anderson motives.

    TESTS::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = GF(5^3)
        sage: M = AndersonMotive(A, K)
        sage: TestSuite(M).run()
    """
    @staticmethod
    def __classcall_private__(cls, category, tau, twist=0, names=None, normalize=True):
        r"""
        Normalize the input and return an instance of the appropriate class.

        INPUT:

        - ``category`` -- the category of Anderson motives where this
          Anderson motive leaves

        - ``tau`` -- a matrix

        - ``twist`` -- an integer (default: ``0``)

        - ``names`` -- a string or a list of strings (default: ``None``),
          the names of the vector of the canonical basis; if ``None``,
          elements will be represented as row vectors

        - ``normalize`` -- a boolean (default: ``True``)

        The action of `\tau` on the Anderson motive will be given by
        the matrix ``tau * (T - z)**(-twist)`` where `T` is the variable
        of the function ring and `z` is its image in the `A`-field.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: M = AndersonMotive(A, K)
            sage: type(M)
            <class 'sage.rings.function_field.drinfeld_modules.anderson_motive.AndersonMotive_general_with_category'>
        """
        AK = category.base()

        # We normalize the inputs
        twist = ZZ(twist)
        tau = tau.change_ring(AK)
        if normalize:
            divisor = category.divisor()
            exponent = Infinity
            for entry in tau.list():
                if not entry:
                    continue
                e = 0
                while entry.degree() > 0 and e < exponent:
                    entry, R = entry.quo_rem(divisor)
                    if R:
                        break
                    e += 1
                exponent = e
                if exponent == 0:
                    break
            if exponent is not Infinity and exponent > 0:
                denom = divisor ** exponent
                tau = tau.parent()([entry // denom for entry in tau.list()])
                twist -= exponent

        names = normalize_names(names, tau.nrows())
        denominator = Factorization([(category.divisor(), twist)])
        ore = category._ore_polring

        # if (isinstance(K, FractionField_1poly_field)
        #     and category.constant_coefficient() == K.gen()):
        #     from sage.rings.function_field.drinfeld_modules.anderson_motive_rational import AndersonMotive_rational
        #     cls = AndersonMotive_rational

        return cls.__classcall__(cls, tau, ore, denominator, names, category)

    def __init__(self, mat, ore, denominator, names, category) -> None:
        r"""
        Initialize this Anderson motive.
        """
        OreModule.__init__(self, mat, ore, denominator, names, category)
        self._initialize_attributes()

    def _initialize_attributes(self):
        r"""
        Set the main attributes to this Anderson motive.

        .. NOTE::

            Separating this method from `__init__` makes it easier
            to call it in subclasses.
        """
        self._tau = self._pseudohom.matrix()
        if self._denominator:
            [(_, self._twist)] = self._denominator
        else:
            self._twist = 0
        self._general_class = AndersonMotive_general
        self._submodule_class = AndersonSubMotive
        self._quotientModule_class = AndersonQuotientMotive

    def __reduce__(self):
        r"""
        Return the necessary arguments to construct this object,
        as per the pickle protocol.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: M = AndersonMotive(A, K)
            sage: loads(dumps(M)) is M
            True
        """
        return self._general_class, (self._category, self._tau, self._twist, self._names, False)

    @lazy_attribute
    def _dettau(self):
        r"""
        Return the leading coefficient of the determinant of `\tau`
        and its degree.

        Only for internal use.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive()
            sage: M._dettau
            (2*z^2 + 3*z + 3, 1)
        """
        det = self._tau.det()
        return det.leading_coefficient(), det.degree()

    def _repr_(self):
        r"""
        Return a string representation of this Anderson motive.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: M = AndersonMotive(A, K)
            sage: M  # indirect doctest
            Anderson motive of rank 1 over Univariate Polynomial Ring in T over Finite Field in z of size 5^3
        """
        s = "Anderson motive "
        if self._names is None:
            s += "of rank %s " % self.rank()
        else:
            s += "<" + ", ".join(self._names) + "> "
        s += "over %s" % self.base()
        return s

    def _latex_(self):
        r"""
        Return a string representation of this Anderson motive.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: M = AndersonMotive(A, K)
            sage: latex(M)  # indirect doctest
            \texttt{Anderson motive of rank } 1\texttt{ over } \Bold{F}_{5^{3}}[T]

        ::

            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M.<u, v> = phi.anderson_motive()
            sage: latex(M)  # indirect doctest
            \left<u, v\right>_{\Bold{F}_{5^{3}}[T]}
        """
        AK = self.base()
        if self._names is None:
            s = "\\texttt{Anderson motive of rank } %s" % self.rank()
            s += "\\texttt{ over } %s" % latex(AK)
        else:
            s = "\\left<" + ", ".join(self._latex_names) + "\\right>"
            s += "_{%s}" % latex(AK)
        return s

    def _Hom_(self, other, category):
        r"""
        Return the set of morphisms from ``self`` to ``other``.

        INPUT:

        - ``other`` -- the codomain of the homset

        - ``category`` -- the category in which we consider the
          morphisms, usually ``self.category()``

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^4)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3, z^4])
            sage: M = phi.anderson_motive()
            sage: End(M)  # indirect doctest
            Set of Morphisms
            from Anderson motive of Drinfeld module defined by T |--> (z^2 + z + 3)*τ^3 + z^3*τ^2 + z^2*τ + z
            to Anderson motive of Drinfeld module defined by T |--> (z^2 + z + 3)*τ^3 + z^3*τ^2 + z^2*τ + z
            in Category of finite dimensional Ore modules with basis
            over Univariate Polynomial Ring in T over Finite Field in z of size 5^4 twisted by T |--> T, with map of base ring
        """
        if category is None:
            category = self._category
        return AndersonMotive_homspace(self, other, category)

    def hodge_pink_weights(self):
        r"""
        Return the Hodge-Pink weights of this Anderson motive,
        sorted by increasing order.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^4)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3, z^4])
            sage: M = phi.anderson_motive()
            sage: M.hodge_pink_weights()
            [0, 0, 1]

        We check that the Hodge-Pink weights of the dual are the opposite
        of the Hodge-Pink weights of the initial Anderson motive::

            sage: N = phi.anderson_motive(dual=True)
            sage: N.hodge_pink_weights()
            [-1, 0, 0]
        """
        S = self._tau.smith_form(transformation=False)
        return [-self._twist + S[i,i].degree() for i in range(self.rank())]

    def is_effective(self):
        r"""
        Return whether this Anderson module is effective, that is,
        whether the action of `\tau` stabilizes it.
        This is also equivalent to the fact that all Hodge-Pink weights
        are nonnegative.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^4)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive()
            sage: M.is_effective()
            True

        ::

            sage: N = phi.anderson_motive(dual=True)
            sage: N.is_effective()
            False
        """
        return self._twist <= 0


class AndersonMotive_drinfeld(AndersonMotive_general):
    r"""
    A class for Anderson motives coming from Drinfeld modules.

    TESTS::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = GF(5^3)
        sage: phi = DrinfeldModule(A, [z, z^2, z^3])
        sage: M = phi.anderson_motive()
        sage: TestSuite(M).run()
    """
    def __classcall_private__(cls, phi, dual, names):
        r"""
        Normalize the input and construct this Anderson motive.

        INPUT:

        - ``phi`` -- a Drinfeld module

        - ``dual`` -- a boolean

        - ``names`` -- a string or a list of strings (default: ``None``),
          the names of the vector of the canonical basis; if ``None``,
          elements will be represented as row vectors

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive(names='e')
            sage: type(M)
            <class 'sage.rings.function_field.drinfeld_modules.anderson_motive.AndersonMotive_drinfeld_with_category'>

        ::

            sage: N.<e0, e1> = phi.anderson_motive()
            sage: M is N
            True
        """
        category = AndersonMotives(phi.category())
        AK = category.base()
        r = phi.rank()
        tau = matrix(AK, r)
        P = phi.gen()
        if dual:
            divisor = category.divisor()
            for i in range(1, r):
                tau[i-1, i] = divisor
                tau[i-1, 0] = P[i]
            tau[r-1, 0] = P[r]
            denominator = Factorization([(divisor, 1)])
        else:
            tau[r-1, 0] = (AK.gen() - P[0]) / P[r]
            for i in range(1, r):
                tau[i-1, i] = 1
                tau[r-1, i] = -P[i]/P[r]
            denominator = Factorization([])
        names = normalize_names(names, r)
        return cls.__classcall__(cls, tau, category._ore_polring, denominator, names, category, phi, dual)

    def __init__(self, mat, ore, denominator, names, category, phi, dual) -> None:
        r"""
        Initialize this Anderson motive.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: tau = phi.ore_variable()
            sage: M = phi.anderson_motive()
            sage: M(tau)
            (0, 1)
        """
        super().__init__(mat, ore, denominator, names, category)
        if not dual:
            Ktau = phi.ore_polring()
            self.register_coercion(DrinfeldToAnderson(Homset(Ktau, self), phi))
            Ktau.register_conversion(AndersonToDrinfeld(Homset(self, Ktau), phi))
        self._drinfeld_module = phi
        self._dual = dual

    def __reduce__(self):
        r"""
        Return the necessary arguments to construct this object,
        as per the pickle protocol.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive()
            sage: loads(dumps(M)) is M
            True

        ::

            sage: Md = phi.anderson_motive(dual=True)
            sage: loads(dumps(Md)) is Md
            True
        """
        return AndersonMotive_drinfeld, (self._drinfeld_module, self._dual, self._names)

    def _repr_(self):
        r"""
        Return a string representation of this Anderson motive.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M.<u, v> = phi.anderson_motive()
            sage: M  # indirect doctest
            Anderson motive <u, v> of Drinfeld module defined by T |--> (2*z + 2)*τ^2 + z^2*τ + z

        ::

            sage: Md = phi.anderson_motive(dual=True)
            sage: Md  # indirect doctest
            Dual Anderson motive of Drinfeld module defined by T |--> (2*z + 2)*τ^2 + z^2*τ + z
        """
        if self._dual:
            s = "Dual Anderson motive "
        else:
            s = "Anderson motive "
        if self._names is not None:
            s += "<" + ", ".join(self._names) + "> "
        s += "of %s" % self._drinfeld_module
        return s

    def drinfeld_module(self):
        r"""
        Return the Drinfeld module from which this Anderson motive
        was constructed.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^5)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive()
            sage: M.drinfeld_module()
            Drinfeld module defined by T |--> z^3*τ^2 + z^2*τ + z
            sage: M.drinfeld_module() is phi
            True
        """
        return self._drinfeld_module


class AndersonSubMotive(AndersonMotive_general, OreSubmodule):
    r"""
    A class for Anderson motives defined as submodules of an
    other Anderson motive.

    TESTS::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = GF(5^3)
        sage: phi = DrinfeldModule(A, [z, z^2, z^3])
        sage: M.<u, v> = phi.anderson_motive()
        sage: N = M.span(v)
        sage: TestSuite(N).run()
    """
    def __init__(self, ambient, submodule, names):
        r"""
        Initialize this Anderson motive.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M.<u, v> = phi.anderson_motive()
            sage: N = M.span(v)
            sage: N.ambient_module() is M
            True

        ::

            sage: type(N)
            <class 'sage.rings.function_field.drinfeld_modules.anderson_motive.AndersonSubMotive_with_category'>
        """
        OreSubmodule.__init__(self, ambient, submodule, names)
        self._initialize_attributes()

    def __reduce__(self):
        r"""
        Return the necessary arguments to construct this object,
        as per the pickle protocol.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M.<u, v> = phi.anderson_motive()
            sage: N = M.span(v)
            sage: loads(dumps(N)) is N
            True
        """
        return OreSubmodule.__reduce__(self)


class AndersonQuotientMotive(AndersonMotive_general, OreQuotientModule):
    r"""
    A class for Anderson motives defined as quotients of an
    other Anderson motive.

    TESTS::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = GF(5^3)
        sage: phi = DrinfeldModule(A, [z, z^2, z^3])
        sage: M.<u, v> = phi.anderson_motive()
        sage: Q = M.quo(u)
        sage: TestSuite(Q).run()
    """
    def __init__(self, cover, submodule, names):
        r"""
        Initialize this Anderson motive.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: M.<u, v> = AndersonMotive(A, diagonal_matrix([1, T - z]))
            sage: Q = M.quo(u)
            sage: Q.cover() is M
            True
            sage: type(Q)
            <class 'sage.rings.function_field.drinfeld_modules.anderson_motive.AndersonQuotientMotive_with_category'>
        """
        OreQuotientModule.__init__(self, cover, submodule, names)
        self._initialize_attributes()

    def __reduce__(self):
        r"""
        Return the necessary arguments to construct this object,
        as per the pickle protocol.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive()
            sage: loads(dumps(M)) is M
            True
        """
        return OreQuotientModule.__reduce__(self)


# Morphisms
###########

# Morphisms between Anderson modules

class AndersonMotiveMorphism(OreModuleMorphism):
    r"""
    A class for morphisms betweeen Anderson motives.

    TESTS::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = GF(5^3)
        sage: phi = DrinfeldModule(A, [z, z^2, z^3])
        sage: u = phi.scalar_multiplication(T)
        sage: f = u.anderson_motive()
        sage: TestSuite(f).run()
    """
    def _repr_type(self):
        r"""
        Return a string representation of the type of this morphism.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: u = phi.scalar_multiplication(T)
            sage: u.anderson_motive()  # indirect doctest
            Endomorphism of Anderson motive of Drinfeld module defined by T |--> (2*z + 2)*τ^2 + z^2*τ + z
        """
        return None

    def __init__(self, parent, im_gens, check=True):
        r"""
        Initialize this morphism.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: M = AndersonMotive(A, K)
            sage: E = End(M)
            sage: f = E(matrix(1, 1, [T]))
            sage: f
            Endomorphism of Anderson motive of rank 1 over Univariate Polynomial Ring in T over Finite Field in z of size 5^3
        """
        from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive_drinfeld
        if isinstance(im_gens, DrinfeldModuleMorphism):
            domain = parent.domain()
            codomain = parent.codomain()
            if not isinstance(domain, AndersonMotive_drinfeld)\
            or domain.drinfeld_module() is not im_gens.codomain():
                raise ValueError("the domain must be the Anderson module of the codomain of the isogeny")
            if not isinstance(codomain, AndersonMotive_drinfeld)\
            or codomain.drinfeld_module() is not im_gens.domain():
                raise ValueError("the codomain must be the Anderson module of the domain of the isogeny")
            im_gens = im_gens._motive_matrix()
            check = False
        OreModuleMorphism.__init__(self, parent, im_gens, check)

    def characteristic_polynomial(self, var='X'):
        r"""
        Return the characteristic polynomial of this morphism.

        INPUT:

        - ``var`` -- a string (default: ``X``), the name of the variable

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: f = phi.scalar_multiplication(T).anderson_motive()
            sage: chi = f.characteristic_polynomial()
            sage: chi
            X^2 + 3*T*X + T^2
            sage: chi.factor()
            (4*X + T)^2

        We compute the characteristic polynomial of the Frobenius and
        compare the result with the output of the method
        :meth:`sage.rings.function_field.drinfeld_modules.drinfeld_module_finite.frobenius_charpoly`::

            sage: Frob = phi.frobenius_endomorphism().anderson_motive()
            sage: Frob.characteristic_polynomial()
            X^2 + X + 3*T^3 + 4*T + 4
            sage: phi.frobenius_charpoly()
            X^2 + X + 3*T^3 + 4*T + 4
        """
        chi = OreModuleMorphism.characteristic_polynomial(self, var)
        A = self.domain().function_ring()
        return chi.change_ring(A)

    charpoly = characteristic_polynomial


class AndersonMotive_homspace(OreModule_homspace):
    Element = AndersonMotiveMorphism


# Coercion maps

class DrinfeldToAnderson(Map):
    r"""
    The canonical isomorphism `K\{\tau\} \to M(\phi)`
    for a Drinfeld module `\phi : A \to K\{\tau\}`.
    """
    def __init__(self, parent, phi):
        r"""
        Initialize this map.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: Ktau = phi.ore_polring()
            sage: M = phi.anderson_motive()
            sage: f = Ktau.convert_map_from(M)
            sage: type(f)
            <class 'sage.rings.function_field.drinfeld_modules.anderson_motive.AndersonToDrinfeld'>
        """
        Map.__init__(self, parent)
        self._drinfeld_module = phi
        self._motive = parent.codomain()
        self._AK = self._motive.base()

    def _call_(self, f):
        r"""
        Return the image of `f` in the Anderson motive.

        INPUT:

        - ``f`` -- a Ore polynomial

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: M = phi.anderson_motive()
            sage: tau = phi.ore_variable()
            sage: M(tau)  # indirect doctest
            (0, 1)
        """
        phi = self._drinfeld_module
        r = phi.rank()
        phiT = phi.gen()
        coords = []
        for _ in range(r):
            coords.append([])
        while f:
            f, rem = f.right_quo_rem(phiT)
            for i in range(r):
                coords[i].append(rem[i])
        coords = [self._AK(c) for c in coords]
        return self._motive(coords)


class AndersonToDrinfeld(Map):
    r"""
    The canonical isomorphism `M(\phi) \to K\{\tau\}`
    for a Drinfeld module `\phi : A \to K\{\tau\}`.
    """
    def __init__(self, parent, phi):
        r"""
        Initialize this map.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^4])
            sage: Ktau = phi.ore_polring()
            sage: M = phi.anderson_motive()
            sage: f = M.coerce_map_from(Ktau)
            sage: type(f)
            <class 'sage.rings.function_field.drinfeld_modules.anderson_motive.DrinfeldToAnderson'>
        """
        Map.__init__(self, parent)
        self._drinfeld_module = phi
        self._Ktau = parent.codomain()

    def _call_(self, x):
        r"""
        Initialize this map.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: K.<z> = GF(5^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^4])
            sage: tau = phi.ore_variable()
            sage: M.<u, v> = phi.anderson_motive()
            sage: u + tau
            u + v
        """
        phi = self._drinfeld_module
        r = phi.rank()
        phiT = phi.gen()
        S = self._Ktau
        xs = x.list()
        ans = S.zero()
        d = max(xi.degree() for xi in xs)
        for j in range(d, -1, -1):
            ans = ans*phiT + S([xs[i][j] for i in range(r)])
        return ans


# Constructor
#############

def AndersonMotive(arg1, arg2=None, names=None):
    r"""
    Construct an Anderson motive.

    INPUT:

    The two first arguments can be one of the followings:

    - a pair `(A, K)` where `A` is the underlying function
      ring (which currently needs to be of the form `\GF{q}[t]`)
      and `K` is the `A`-field; these parameters correspond to
      the trivial Anderson motive over `A \otimes K`

    - a pair `(A, z)` where `A = \GF{q}[t]` is the function
      base ring and `z` is an element; the `A`-field is then
      then parent `K` of `z` viewed as an algebra over `A`
      through `A \mapsto K, T \mapsto z`.

    - a pair `(A, \tau)` where

      - `A` is either `\GF{q}[t]` or a category (of Drinfeld
        modules or Anderson motives)

      - `\tau` is the matrix defining the Anderson motive

    - a Drinfeld module

    EXAMPLES::

        sage: A.<T> = GF(7)[]
        sage: K.<z> = GF(7^3)

    We first construct the trivial Anderson motive over `K`::

        sage: M = AndersonMotive(A, K)
        sage: M
        Anderson motive of rank 1 over Univariate Polynomial Ring in T over Finite Field in z of size 7^3
        sage: M.matrix()
        [1]

    Here the structure of `A`-field on `K` is given by the map
    that takes `T` to the canonical generator of `K`, namely `z`::

        sage: M.A_field()
        Finite Field in z of size 7^3 over its base
        sage: M.A_field().defining_morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field of size 7
          To:   Finite Field in z of size 7^3 over its base
          Defn: T |--> z

    Specifying another element in `K` leads to a different
    structure of `A`-field::

        sage: N = AndersonMotive(A, z^2)
        sage: N.A_field().defining_morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field of size 7
          To:   Finite Field in z of size 7^3 over its base
          Defn: T |--> z^2

    One can also directly construct the Anderson motive attached
    to a Drinfeld module as follows::

        sage: phi = DrinfeldModule(A, [z, z^2, z^3])
        sage: AndersonMotive(phi)
        Anderson motive of Drinfeld module defined by T |--> (z^2 + 3)*τ^2 + z^2*τ + z

    Finally, another possibility is to give the matrix of `\tau` as an
    argument::

        sage: tau = matrix(2, 2, [[T, z], [z+1, 1]])
        sage: tau
        [    T     z]
        [z + 1     1]
        sage: M = AndersonMotive(A, tau)
        sage: M
        Anderson motive of rank 2 over Univariate Polynomial Ring in T over Finite Field in z of size 7^3
        sage: M.matrix()
        [    T     z]
        [z + 1     1]

    In this case, the structure of `A`-field is automatically inferred::

        sage: M.A_field().defining_morphism()
        Ring morphism:
          From: Univariate Polynomial Ring in T over Finite Field of size 7
          To:   Finite Field in z of size 7^3 over its base
          Defn: T |--> z^2 + z

    TESTS::

        sage: gamma = A.hom([z+1])
        sage: AndersonMotive(A, gamma)
        Anderson motive of rank 1 over Univariate Polynomial Ring in T over Finite Field in z of size 7^3

    ::

        sage: AndersonMotive(phi, tau)
        Traceback (most recent call last):
        ...
        ValueError: too many arguments

    ::

        sage: AndersonMotive(ZZ, K)
        Traceback (most recent call last):
        ...
        ValueError: unable to parse arguments

    ::

        sage: tau = matrix(2, 2, [[T^2, z], [z+1, 1]])
        sage: AndersonMotive(A, tau)
        Traceback (most recent call last):
        ...
        ValueError: the given matrix does not define an Anderson motive in this category

    .. SEEALSO::

        :mod:`sage.rings.function_field.drinfeld_modules.anderson_motive`
    """
    # Options for *args:
    #  . a Drinfeld module
    #  . a category (of Drinfeld modules or AndersonMotives)
    #  . a ring, a matrix
    #  . a ring, a A-field

    # We first try to parse the first argument
    # We have several cases:
    # (1) arg1 is a Drinfeld module
    if isinstance(arg1, DrinfeldModule):
        if arg2 is not None:
            raise ValueError("too many arguments")
        return AndersonMotive_drinfeld(arg1, False, names=names)

    # (2) arg1 is a category
    category = None
    if isinstance(arg1, DrinfeldModules):
        category = AndersonMotives(arg1)
    if isinstance(arg1, AndersonMotives):
        category = arg1
    if category is not None:
        return category.object(arg2, names=names)

    # (3) arg1 is the function ring
    if not isinstance(arg1, PolynomialRing_general):
        raise ValueError("unable to parse arguments")
    A = arg1
    # This case requires to parse the second argument
    # Again we have several cases:
    # (3a) arg2 encodes the base field
    K = None
    if isinstance(arg2, RingHomomorphism) and arg2.domain() is A:
        # (3a-i) arg2 is the morphism of the form A -> K
        K = arg2.codomain()
        gamma = arg2
    elif isinstance(arg2, CommutativeRing):
        # (3a-ii) arg2 is the A-field K
        K = arg2
        if K.has_coerce_map_from(A):
            gamma = K.coerce_map_from(A)
        else:
            gamma = A.hom([K.gen()])
    elif hasattr(arg2, 'parent') and isinstance(arg2.parent(), CommutativeRing):
        # (3a-iii) arg2 is an element in K
        K = arg2.parent()
        gamma = A.hom([arg2])
    if K is not None:
        # The case (3a) was successful:
        # We construct and return the Anderson motive
        category = AndersonMotives(gamma)
        return category.object(names=names)

    # (3b) arg2 is a matrix given the action of tau
    if isinstance(arg2, Matrix):
        tau = arg2
        AK = tau.base_ring()
        if not isinstance(AK, PolynomialRing_general) or AK.variable_name() != A.variable_name():
            raise TypeError("incompatible base rings")
        det = tau.determinant()
        if det == 0:
            raise ValueError("tau does not define an Anderson motive")
        h = det.degree()
        K = AK.base_ring()
        gamma = K.coerce_map_from(A)
        if gamma is None:
            p = A.characteristic()
            if h.gcd(p) == 1:
                theta = -det[h-1] / det[h] / h
            else:
                raise NotImplementedError("cannot determine the structure of A-field")
            gamma = A.hom([theta])
        category = AndersonMotives(gamma)
        return category.object(tau, names=names)

    raise ValueError("unable to parse arguments")
