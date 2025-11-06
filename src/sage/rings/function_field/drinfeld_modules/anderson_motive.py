r"""
Anderson motives

Let `\GF{q}[T]` be a polynomial ring with coefficients in a finite
field `\GF{q}` and let `K` be an extension of `\GF{q}` equipped
with a distinguished element `z`.

By definition, an Anderson motive attached to these data is a free
module of finite rank `M` over `K[T]`, equipped with a linear
automorphism

.. MATH::

    tau_M : \tau^\star M \left[\frac 1[T-z]\right] \to M \left[\frac 1[T-z]\right]

where `\tau^\star M = K \otimes_{K, \text{Frob}} M`.

Any Drinfeld module `phi` over `(A, \gamma)` with `gamma : A \to K,
T \mapsto z` gives rise to an Anderson motive. By definition, it is
`M(\phi) := K\{\tau\}` (the ring of Ore polynomials with commutation
rule `\tau \lambda = \lambda^q \tau` for `\lambda \in K`) where

- the structure of `\GF{q}[T]`-module is given by right multiplication
  by `phi_a` (`a \in \GF{q}[T]`),

- the structure of `K`-module is given by left multiplication,

- the automorphism `\tau_{M(\phi)}` is the left multiplication
  by `tau` in the Ore polynomial ring.

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
    Anderson motive of rank 3 over Univariate Polynomial Ring in T over Finite Field in z of size 5^3

We see that `M` has rank `3`; it is actually a general fact that
the Anderson motive attached a Drinfeld module has the same rank
than the underlying Drinfeld module.

The canonical basis corresponds to the vectors `tau^i` for `i`
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

SageMath provides facilities to pick elements in `M` and perform
basic operations with them::

    sage: u, v, w = M.basis()
    sage: T*u + z*w
    (T, 0, z)
    sage: w.image()  # image by tau_M
    ((z^2 + 3*z)*T + 2*z^2 + 3*z + 3, 3*z^2 + 2*z + 4, 2*z^2 + 1)

Some basic constructions on Anderson modules are also available.
For example, one can form the dual::

    sage: Md = M.dual()
    sage: Md
    Anderson motive of rank 3 over Univariate Polynomial Ring in T over Finite Field in z of size 5^3
    sage: Md.matrix()
    [          z^2/(T + 4*z)                       1                       0]
    [    (2*z + 2)/(T + 4*z)                       0                       1]
    [(2*z^2 + 2*z)/(T + 4*z)                       0                       0]

or Carlitz twists::

    sage: M2 = M.carlitz_twist(2)
    sage: M2
    Anderson motive of rank 3 over Univariate Polynomial Ring in T over Finite Field in z of size 5^3
    sage: M2.matrix()
    [                                    0                 1/(T^2 + 3*z*T + z^2)                                     0]
    [                                    0                                     0                 1/(T^2 + 3*z*T + z^2)]
    [                (z^2 + 3*z)/(T + 4*z) (3*z^2 + 2*z + 4)/(T^2 + 3*z*T + z^2)       (2*z^2 + 1)/(T^2 + 3*z*T + z^2)]

We observe that the entries of the previous matrices have denominators which are
`T-z` or powers of it. This corresponds to the fact that `\tau_M` is only defined
after inverting `T-z` in full generality.

AUTHOR:

- Xavier Caruso (2025-11): initial version
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

import operator

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.latex import latex
from sage.misc.functional import log

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
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField_1poly_field

from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix, block_diagonal_matrix

from sage.modules.ore_module import OreModule, OreSubmodule, OreQuotientModule
from sage.modules.ore_module import OreAction
from sage.modules.ore_module import normalize_names
from sage.modules.ore_module_element import OreModuleElement
from sage.modules.ore_module_homspace import OreModule_homspace
from sage.modules.ore_module_morphism import OreModuleMorphism

from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism


# Classes for Anderson motives
##############################

class AndersonMotiveElement(OreModuleElement):
    def image(self, integral=None):
        if integral is None:
            integral = self.parent().is_effective()
        return super().image(integral=integral)


class AndersonMotive_general(OreModule):
    Element = AndersonMotiveElement

    @staticmethod
    def __classcall_private__(self, category, tau, twist=0, names=None, normalize=True):
        K = category.base()
        AK = category.base_combined()

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

        #if (isinstance(K, FractionField_1poly_field)
        #    and category.constant_coefficient() == K.gen()):
        #    from sage.rings.function_field.drinfeld_modules.anderson_motive_rational import AndersonMotive_rational
        #    cls = AndersonMotive_rational
        #else:
        cls = AndersonMotive_general

        return cls.__classcall__(cls, tau, ore, denominator, names, category)

    def __init__(self, mat, ore, denominator, names, category) -> None:
        OreModule.__init__(self, mat, ore, denominator, names, category)
        self._initialize_attributes()

    def _initialize_attributes(self):
        category = self._category
        self._A = A = category.function_ring()
        self._t_name = A.variable_name()
        self._Fq = Fq = A.base_ring()
        self._q = Fq.cardinality()
        self._deg = ZZ(log(self._q, Fq.characteristic()))
        self._K = self._base = K = category.base()
        self._theta = category.constant_coefficient()
        self._AK = base = category.base_combined()
        self._t = base.gen()
        self._tau = self._pseudohom.matrix()
        if self._denominator:
            self._twist = self._denominator[0][1]
        else:
            self._twist = 0
        self._submodule_class = AndersonSubMotive
        self._quotientModule_class = AndersonQuotientMotive

    @lazy_attribute
    def _dettau(self):
        det = self._tau.det()
        return det.leading_coefficient(), det.degree()

    def _repr_(self):
        s = "Anderson motive "
        if self._names is None:
            s += "of rank %s " % self.rank()
        else:
            s += "<" + ", ".join(self._names) + "> "
        s += "over %s" % self._AK
        return s

    def _latex_(self):
        if self._names is None:
            s = "\\texttt{Anderson motive of rank } %s" % self.rank()
            s += "\\texttt{ over } %s" % latex(self._AK)
        else:
            s = "\\left<" + ", ".join(self._latex_names) + "\\right>"
            s += "_{%s}" % latex(self._AK)
        return s

    def carlitz_twist(self, n, names=None):
        return AndersonMotive_general(self._category, self._tau, self._twist + ZZ(n),
                                      names, normalize=False)

    def dual(self, names=None):
        disc, deg = self._dettau
        scalar = self._K(~disc)
        tau = scalar * self._tau.adjugate().transpose()
        twist = deg - self._twist
        return AndersonMotive_general(self._category, tau, twist, names, normalize=True)

    def _Hom_(self, codomain, category):
        from sage.rings.function_field.drinfeld_modules.anderson_motive_morphism import AndersonMotive_homspace
        return AndersonMotive_homspace(self, codomain)

    def hodge_pink_weights(self):
        S = self._tau.smith_form(transformation=False)
        return [-self._twist + S[i,i].degree() for i in range(self.rank())]

    def is_effective(self):
        return self._twist <= 0

    def ore_variable(self):
        return self._category._ore_polring.gen()

    def ore_polring(self, names=None, action=True):
        if names is None:
            names = self._category._ore_variable_name
        S = self._ore_category.ore_ring(names)
        if action:
            self._unset_coercions_used()
            self.register_action(OreAction(S, self, True, operator.mul))
        return S


class AndersonMotive_drinfeld(AndersonMotive_general):
    def __init__(self, phi, names):
        category = AndersonMotives(phi.category())
        AK = category.base_combined()
        r = phi.rank()
        tau = matrix(AK, r)
        P = phi.gen()
        tau[r-1, 0] = (AK.gen() - P[0]) / P[r]
        for i in range(1, r):
            tau[i-1, i] = 1
            tau[r-1, i] = -P[i]/P[r]
        names = normalize_names(names, r)
        AndersonMotive_general.__init__(self, tau, category._ore_polring, None, names, category)
        Ktau = phi.ore_polring()
        self.register_coercion(DrinfeldToAnderson(Homset(Ktau, self), phi))
        try:
            Ktau.register_conversion(AndersonToDrinfeld(Homset(self, Ktau), phi))
        except AssertionError:
            pass
        self._drinfeld_module = phi

    def drinfeld_module(self):
        return self._drinfeld_module


class AndersonSubMotive(AndersonMotive_general, OreSubmodule):
    def __init__(self, ambient, submodule, names):
        OreSubmodule.__init__(self, ambient, submodule, names)
        self._initialize_attributes()


class AndersonQuotientMotive(AndersonMotive_general, OreQuotientModule):
    def __init__(self, cover, submodule, names):
        OreQuotientModule.__init__(self, cover, submodule, names)
        self._initialize_attributes()


# Morphisms
###########

# Morphisms between Anderson modules

class AndersonMotiveMorphism(OreModuleMorphism):
    def _repr_type(self):
        return "Anderson motive"

    def __init__(self, parent, im_gens, check=True):
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
            u = im_gens._ore_polynomial
            im_gens = {codomain.gen(0): u*domain.gen(0)}
            check = False
        OreModuleMorphism.__init__(self, parent, im_gens, check)

    def characteristic_polynomial(self, var='X'):
        chi = OreModuleMorphism.characteristic_polynomial(self, var)
        A = self.domain().function_ring()
        return chi.change_ring(A)

    charpoly = characteristic_polynomial


class AndersonMotive_homspace(OreModule_homspace):
    Element = AndersonMotiveMorphism


# Coercion maps

class DrinfeldToAnderson(Map):
    def __init__(self, parent, phi):
        Map.__init__(self, parent)
        self._phi = phi
        self._motive = parent.codomain()
        self._AK = self._motive.base_combined()

    def _call_(self, f):
        phi = self._phi
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
    def __init__(self, parent, phi):
        Map.__init__(self, parent)
        self._phi = phi
        self._Ktau = parent.codomain()

    def _call_(self, x):
        phi = self._phi
        r = phi.rank()
        phiT = phi.gen()
        S = self._Ktau
        xs = []
        for i in range(r):
            if x[i].denominator() != 1:
                raise ValueError("not in the Anderson motive")
            xs.append(x[i].numerator())
        ans = S.zero()
        d = max(xi.degree() for xi in xs)
        for j in range(d, -1, -1):
            ans = ans*phiT + S([xs[i][j] for i in range(r)])
        return ans


# Constructor
#############

def AndersonMotive(arg1, tau=None, names=None):
    r"""
    Construct an Anderson motive

    INPUT:

    The input can be one of the followings:

    - a Drinfeld module

    - a pair `(A, \tau)` where

      - `A` is either the underlying function ring (which
        currently needs to be of the form `\FF_q[t]`) or
        a category (of Drinfeld modules or Anderson motives)

      - `\tau` is the matrix defining the Anderson motive

    - a pair '(A, K)` where `A = \FF_q[t]` is the function
      base ring and `K` is the coefficient `A`-field; these
      parameters correspond to the trivial Anderson motive
      over `A \otimes K`

    OUTPUT:

    An anderson motive

    EXAMPLES::


    """
    # Options for *args:
    #  . a Drinfeld module
    #  . a category (of Drinfeld modules or AndersonMotives)
    #  . a ring, a matrix
    #  . a ring, a A-field
    # arg1 is a Drinfeld module
    if isinstance(arg1, DrinfeldModule):
        if tau is not None:
            raise ValueError("")
        category = AndersonMotives(arg1.category())
        A = category.function_ring()
        K = category._base_field
        AK = A.change_ring(K)
        r = arg1.rank()
        tau = matrix(AK, r)
        P = arg1.gen()
        tau[r-1, 0] = (AK.gen() - P[0]) / P[r]
        for i in range(1, r):
            tau[i-1, i] = 1
            tau[r-1, i] = -P[i]/P[r]
        return AndersonMotive_general(category, tau, names=names)

    # arg1 is a category
    category = None
    if isinstance(arg1, DrinfeldModules):
        category = AndersonMotives(arg1)
    if isinstance(arg1, AndersonMotives):
        category = arg1
    if category is not None:
        if tau is None:
            tau = identity_matrix(category.base_combined(), 1)
        det = tau.determinant()
        if det == 0:
            raise ValueError("tau does not define an Anderson motive")
        h = det.degree()
        disc, R = det.quo_rem(category.divisor() ** h)
        if R:
            raise ValueError("tau does not define an Anderson motive")
        M = AndersonMotive_general(category, tau, names=names)
        #M._set_dettau(disc[0], h, 0)
        return M

    # arg1 is the function ring
    if isinstance(arg1, CommutativeRing):
        A = arg1
        if not isinstance(A, PolynomialRing_general):
            raise NotImplementedError("Anderson motives over arbitrary Dedekind domain are not supported")
    else:
        raise ValueError("first argument must be the function ring")

    # tau is the base ring
    K = None
    if isinstance(tau, RingHomomorphism) and tau.domain() is A:
        K = tau.codomain()
        gamma = tau
    elif isinstance(tau, CommutativeRing) and tau.has_coerce_map_from(A):
        K = tau
        gamma = K.coerce_map_from(A)
    elif hasattr(tau, 'parent') and isinstance(tau.parent(), CommutativeRing):
        K = tau.parent()
        gamma = A.hom([tau])
    if K is not None:
        try:
            if K.variable_name() == A.variable_name():
                K = K.base_ring()
        except (AttributeError, ValueError):
            pass
        category = AndersonMotives(gamma)
        AK = category.base_combined()
        tau = identity_matrix(AK, 1)
        return AndersonMotive_general(category, tau, names=names)

    # tau is a matrix
    if isinstance(tau, Matrix):
        AK = tau.base_ring()
        if not isinstance(AK, PolynomialRing_general) or AK.variable_name() != A.variable_name():
            raise ValueError("incompatible base rings")
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
        disc, R = det.quo_rem(category.divisor() ** h)
        if R:
            raise ValueError("tau does not define an Anderson motive")
        M = AndersonMotive_general(category, tau, names=names)
        #M._set_dettau(disc[0], h, 0)
        return M

    raise ValueError("unable to parse arguments")
