# sage.doctest: optional - giac
"""
Wrappers for Giac functions

We provide a python function to compute and convert to sage a Groebner
basis.

AUTHORS:

- Martin Albrecht (2015-07-01): initial version
- Han Frederic (2015-07-01): initial version

EXAMPLES:

Compute and verify a Groebner basis::

    sage: from sage.libs.giac.gb import groebner_basis as gb_giac
    sage: from sage.rings.ideal import Cyclic as CyclicIdeal
    sage: from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    sage: from sage.rings.rational_field import QQ
    sage: P = PolynomialRing(QQ, 6, 'x')
    sage: I = CyclicIdeal(P)
    sage: B = gb_giac(I.gens())  # random
    sage: B
    Polynomial Sequence with 45 Polynomials in 6 Variables
    sage: B.is_groebner()
    True

"""

# *****************************************************************************
#       Copyright (C) 2013 Frederic Han <frederic.han@imj-prg.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.proof.all import polynomial as proof_polynomial
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.libs.giac.context import local_giacsettings
from sage.libs.giac.giac import giacsettings, libgiac


@local_giacsettings
def groebner_basis(gens, proba_epsilon=None, threads=None, prot=False,
                   elim_variables=None, *args, **kwds):
    r"""
    Compute a Groebner Basis of an ideal using ``giacpy_sage``. The result is
    automatically converted to sage.

    Supported term orders of the underlying polynomial ring are ``lex``,
    ``deglex``, ``degrevlex`` and block orders with 2 ``degrevlex`` blocks.

    INPUT:

    - ``gens`` -- an ideal (or a list) of polynomials over a prime field
      of characteristic 0 or `p < 2^31`

    - ``proba_epsilon`` -- (default: ``None``) majoration of the probability
      of a wrong answer when probabilistic algorithms are allowed

        * if ``proba_epsilon`` is None, the value of
          ``sage.structure.proof.all.polynomial()`` is taken. If it is
          false then the global ``giacpy_sage.giacsettings.proba_epsilon`` is
          used.

        * if ``proba_epsilon`` is 0, probabilistic algorithms are
          disabled.

    - ``threads`` -- (default: ``None``) maximal number of threads allowed
      for giac. If ``None``, the global ``giacpy_sage.giacsettings.threads`` is
      considered.

    - ``prot`` -- boolean (default: ``False``); if ``True`` print detailed information

    - ``elim_variables`` -- (default: ``None``) a list of variables to eliminate
      from the ideal

        * if ``elim_variables`` is None, a Groebner basis with respect to the
          term ordering of the parent polynomial ring of the polynomials
          ``gens`` is computed.

        * if ``elim_variables`` is a list of variables, a Groebner basis of the
          elimination ideal with respect to a ``degrevlex`` term order is
          computed, regardless of the term order of the polynomial ring.

    OUTPUT: polynomial sequence of the reduced Groebner basis

    EXAMPLES::

        sage: from sage.arith.misc import previous_prime
        sage: from sage.libs.giac.gb import groebner_basis as gb_giac
        sage: from sage.rings.finite_rings.finite_field_constructor import GF
        sage: from sage.rings.ideal import Cyclic as CyclicIdeal
        sage: from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        sage: P = PolynomialRing(GF(previous_prime(2**31)), 6, 'x')
        sage: I = CyclicIdeal(P)
        sage: B = gb_giac(I.gens())
        sage: B
        Polynomial Sequence with 45 Polynomials in 6 Variables
        sage: B.is_groebner()
        True

    Elimination ideals can be computed by passing ``elim_variables``::

        sage: from sage.libs.giac.gb import groebner_basis as gb_giac
        sage: P = PolynomialRing(GF(previous_prime(2**31)), 5, 'x')
        sage: I = CyclicIdeal(P)
        sage: B = gb_giac(I.gens(), elim_variables=[P.gen(0), P.gen(2)])
        sage: B.is_groebner()
        True
        sage: B.ideal() == I.elimination_ideal([P.gen(0), P.gen(2)])
        True

    Computations over QQ can benefit from a probabilistic lifting::

        sage: from sage.rings.rational_field import QQ
        sage: from sage.rings.ideal import Ideal
        sage: from sage.structure.proof.all import polynomial as proof_polynomial
        sage: P = PolynomialRing(QQ, 5, 'x')
        sage: I = Ideal([P.random_element(3,7) for j in range(5)])
        sage: B1 = gb_giac(I.gens(),1e-16)  # random
        sage: proof_polynomial(True)
        sage: B2 = gb_giac(I.gens())  # random
        sage: B1 == B2
        True
        sage: B1.is_groebner()
        True

    You can get detailed information by setting ``prot=True``, but
    it won't appear in the doctest output because C libraries are
    missed by python's doctest input/output redirection::

        sage: from sage.rings.ideal import Katsura as KatsuraIdeal
        sage: P = PolynomialRing(QQ, 8, 'x')
        sage: I = KatsuraIdeal(P)
        sage: B = gb_giac(I, prot=True)  # random
        sage: B
        Polynomial Sequence with 74 Polynomials in 8 Variables

    TESTS::

        sage: from sage.libs.giac.giac import libgiac
        sage: libgiac("x2:=22; x4:='whywouldyoudothis'")
        22,whywouldyoudothis
        sage: gb_giac(I)
        Traceback (most recent call last):
        ...
        ValueError: Variables names ['x2', 'x4'] conflict in giac. Change them or purge them from in giac with libgiac.purge('x2')
        sage: libgiac.purge('x2'),libgiac.purge('x4')
        (22, whywouldyoudothis)
        sage: gb_giac(I)
        Polynomial Sequence with 74 Polynomials in 8 Variables

        sage: I = Ideal(P(0),P(0))
        sage: I.groebner_basis() == gb_giac(I)
        True

    Test the supported term orderings::

        sage: from sage.rings.ideal import Cyclic as CyclicIdeal
        sage: P = PolynomialRing(QQ, 'x', 4, order='lex')
        sage: B = gb_giac(CyclicIdeal(P))
        ...
        sage: B.is_groebner(), B.ideal() == CyclicIdeal(P)
        (True, True)
        sage: P = P.change_ring(order='deglex')
        sage: B = gb_giac(CyclicIdeal(P))
        ...
        sage: B.is_groebner(), B.ideal() == CyclicIdeal(P)
        (True, True)
        sage: P = P.change_ring(order='degrevlex(2),degrevlex(2)')
        sage: B = gb_giac(CyclicIdeal(P))
        ...
        sage: B.is_groebner(), B.ideal() == CyclicIdeal(P)
        (True, True)
    """
    try:
        iter(gens)
    except TypeError:
        gens = gens.gens()

    # get the ring from gens
    P = next(iter(gens)).parent()
    K = P.base_ring()
    p = K.characteristic()

    # check if the ideal is zero. (giac 1.2.0.19 segfault)
    from sage.rings.ideal import Ideal
    if (Ideal(gens)).is_zero():
        return PolynomialSequence([P(0)], P, immutable=True)

    # check for name confusions
    blackgiacconstants = ['i', 'e'] # NB e^k is expanded to exp(k)
    blacklist = blackgiacconstants + [str(j) for j in libgiac.VARS()]
    problematicnames = sorted(set(P.gens_dict()).intersection(blacklist))

    if problematicnames:
        raise ValueError("Variables names %s conflict in giac. Change them or purge them from in giac with libgiac.purge(\'%s\')"
                         % (problematicnames, problematicnames[0]))

    if K.is_prime_field() and p == 0:
        F = libgiac(gens)
    elif K.is_prime_field() and p < 2**31:
        F = (libgiac(gens) % p)
    else:
        raise NotImplementedError("Only prime fields of cardinal < 2^31 are implemented in Giac for Groebner bases.")

    # proof or probabilistic reconstruction
    if proba_epsilon is None:
        if proof_polynomial():
            giacsettings.proba_epsilon = 0
        else:
            giacsettings.proba_epsilon = 1e-15
    else:
        giacsettings.proba_epsilon = proba_epsilon

    # prot
    if prot:
        libgiac('debug_infolevel(2)')

    # threads
    if threads is not None:
        giacsettings.threads = threads

    if elim_variables is None:
        var_names = P.variable_names()
        order_name = P.term_order().name()
        if order_name == "degrevlex":
            giac_order = "revlex"
        elif order_name == "lex":
            giac_order = "plex"
        elif order_name == "deglex":
            giac_order = "tdeg"
        else:
            blocks = P.term_order().blocks()
            if (len(blocks) == 2 and
                    all(order.name() == "degrevlex" for order in blocks)):
                giac_order = "revlex"
                var_names = var_names[:len(blocks[0])]
            else:
                raise NotImplementedError(
                        "%s is not a supported term order in "
                        "Giac Groebner bases." % P.term_order())

        # compute de groebner basis with giac
        gb_giac = F.gbasis(list(var_names), giac_order)

    else:
        gb_giac = F.eliminate(list(elim_variables), 'gbasis')

    return PolynomialSequence(gb_giac, P, immutable=True)
