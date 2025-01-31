# sage.doctest: needs sage.libs.giac
"""
Wrappers for Giac functions

We provide a python function to compute and convert to sage a Groebner
basis using the ``giacpy_sage`` module.

AUTHORS:

- Martin Albrecht (2015-07-01): initial version
- Han Frederic (2015-07-01): initial version

EXAMPLES::

    sage: from sage.libs.giac import groebner_basis as gb_giac # random
    sage: P = PolynomialRing(QQ, 6, 'x')
    sage: I = sage.rings.ideal.Cyclic(P)
    sage: B = gb_giac(I.gens()) # random
    sage: B
    Polynomial Sequence with 45 Polynomials in 6 Variables
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
from .giac import giacsettings, libgiac

#  Remarks for doctests:
#     1) The first time that the c++ library giac is loaded a message appears.
#        This message is version and arch dependent.
#     2) When proba_epsilon is too bad (>1e-6?) setting it to a better value
#        will give an additional message like the following one:
#       Restoring proba epsilon to 1e-6 from 1e-12
#       (it looks like  in internal giac changes this also to not work with a too bad probability)


class GiacSettingsDefaultContext:
    """
    Context preserve libgiac settings.
    """

    def __enter__(self):
        """
        EXAMPLES::

           sage: from sage.libs.giac import GiacSettingsDefaultContext
           sage: from sage.libs.giac.giac import giacsettings
           sage: giacsettings.proba_epsilon = 1e-16
           sage: with GiacSettingsDefaultContext(): giacsettings.proba_epsilon = 1e-12
           sage: giacsettings.proba_epsilon < 1e-14
           True
        """
        self.proba_epsilon = giacsettings.proba_epsilon
        self.threads = giacsettings.threads
        # Change the debug level at the end to not have messages at each modification
        self.debuginfolevel = libgiac('debug_infolevel()')

    def __exit__(self, typ, value, tb):
        """
        EXAMPLES::

           sage: from sage.libs.giac import GiacSettingsDefaultContext
           sage: from sage.libs.giac.giac import giacsettings
           sage: giacsettings.proba_epsilon = 1e-16
           sage: with GiacSettingsDefaultContext(): giacsettings.proba_epsilon = 1e-30
           sage: giacsettings.proba_epsilon < 1e-20
           False
        """
        # Restore the debug level first to not have messages at each modification
        libgiac('debug_infolevel')(self.debuginfolevel)
        # NB: giacsettings.epsilon has a different meaning that giacsettings.proba_epsilon.
        giacsettings.proba_epsilon = self.proba_epsilon
        giacsettings.threads = self.threads


def local_giacsettings(func):
    """
    Decorator to preserve Giac's proba_epsilon and threads settings.

    EXAMPLES::

        sage: def testf(a, b):
        ....:    giacsettings.proba_epsilon = a/100
        ....:    giacsettings.threads = b+2
        ....:    return (giacsettings.proba_epsilon, giacsettings.threads)

        sage: from sage.libs.giac.giac import giacsettings
        sage: from sage.libs.giac import local_giacsettings
        sage: gporig, gtorig = (giacsettings.proba_epsilon,giacsettings.threads)
        sage: gp, gt = local_giacsettings(testf)(giacsettings.proba_epsilon,giacsettings.threads)
        sage: gporig == giacsettings.proba_epsilon
        True
        sage: gtorig == giacsettings.threads
        True
        sage: gp<gporig, gt-gtorig
        (True, 2)
    """
    from sage.misc.decorators import sage_wraps

    @sage_wraps(func)
    def wrapper(*args, **kwds):
        """
        Execute function in ``GiacSettingsDefaultContext``.
        """
        with GiacSettingsDefaultContext():
            return func(*args, **kwds)

    return wrapper


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
      of characteristic 0 or `p<2^31`

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

        sage: from sage.libs.giac import groebner_basis as gb_giac
        sage: P = PolynomialRing(GF(previous_prime(2**31)), 6, 'x')
        sage: I = sage.rings.ideal.Cyclic(P)
        sage: B = gb_giac(I.gens())
        ...
        sage: B
        Polynomial Sequence with 45 Polynomials in 6 Variables
        sage: B.is_groebner()
        True

    Elimination ideals can be computed by passing ``elim_variables``::

        sage: P = PolynomialRing(GF(previous_prime(2**31)), 5, 'x')
        sage: I = sage.rings.ideal.Cyclic(P)
        sage: B = gb_giac(I.gens(), elim_variables=[P.gen(0), P.gen(2)])
        ...
        sage: B.is_groebner()
        True
        sage: B.ideal() == I.elimination_ideal([P.gen(0), P.gen(2)])
        True

    Computations over QQ can benefit from

    * a probabilistic lifting::

        sage: P = PolynomialRing(QQ,5, 'x')
        sage: I = ideal([P.random_element(3,7) for j in range(5)])
        sage: B1 = gb_giac(I.gens(),1e-16) # long time (1s)
        ...
        sage: sage.structure.proof.all.polynomial(True)
        sage: B2 = gb_giac(I.gens()) # long time (4s)
        ...
        sage: B1 == B2 # long time
        True
        sage: B1.is_groebner() # not tested, too long time (50s)
        True

    * multi threaded operations::

        sage: P = PolynomialRing(QQ, 8, 'x')
        sage: I = sage.rings.ideal.Cyclic(P)
        sage: time B = gb_giac(I.gens(),1e-6,threads=2) # doctest: +SKIP
        ...
        Time: CPU 168.98 s, Wall: 94.13 s

    You can get detailed information by setting ``prot=True``

    ::

        sage: I = sage.rings.ideal.Katsura(P)
        sage: gb_giac(I,prot=True)  # random, long time (3s)
        9381383 begin computing basis modulo 535718473
        9381501 begin new iteration zmod, number of pairs: 8, base size: 8
        ...end, basis size 74 prime number 1
        G=Vector [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,...
        ...creating reconstruction #0
        ...
        ++++++++basis size 74
        checking pairs for i=0, j=
        checking pairs for i=1, j=2,6,12,17,19,24,29,34,39,42,43,48,56,61,64,69,
        ...
        checking pairs for i=72, j=73,
        checking pairs for i=73, j=
        Number of critical pairs to check 373
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++...
        Successful... check of 373 critical pairs
        12380865 end final check
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
        sage: gb_giac(I) # long time (3s)
        ...Polynomial Sequence with 74 Polynomials in 8 Variables

        sage: I = ideal(P(0),P(0))
        sage: I.groebner_basis() == gb_giac(I)
        True

    Test the supported term orderings::

        sage: from sage.rings.ideal import Cyclic
        sage: P = PolynomialRing(QQ, 'x', 4, order='lex')
        sage: B = gb_giac(Cyclic(P))
        ...
        sage: B.is_groebner(), B.ideal() == Cyclic(P)
        (True, True)
        sage: P = P.change_ring(order='deglex')
        sage: B = gb_giac(Cyclic(P))
        ...
        sage: B.is_groebner(), B.ideal() == Cyclic(P)
        (True, True)
        sage: P = P.change_ring(order='degrevlex(2),degrevlex(2)')
        sage: B = gb_giac(Cyclic(P))
        ...
        sage: B.is_groebner(), B.ideal() == Cyclic(P)
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
    blackgiacconstants = ['i', 'e']  # NB e^k is expanded to exp(k)
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
