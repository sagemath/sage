r"""
Constructor for Anderson motives
"""

# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************


from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.anderson_motives import AndersonMotives

from sage.rings.ring import CommutativeRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.morphism import RingHomomorphism

from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive_general


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
        M._set_dettau(disc[0], h, 0)
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
    elif isinstance(tau, CommutativeRing):
        K = tau
        gamma = A
    if K is not None:
        try:
            if K.variable_name() == A.variable_name():
                K = K.base_ring()
        except (AttributeError, ValueError):
            pass
        category = AndersonMotives(K.over(gamma))
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
        category = AndersonMotives(K.over(gamma))
        disc, R = det.quo_rem(category.divisor() ** h)
        if R:
            raise ValueError("tau does not define an Anderson motive")
        M = AndersonMotive_general(category, tau, names=names)
        M._set_dettau(disc[0], h, 0)
        return M

    raise ValueError("unable to parse arguments")
