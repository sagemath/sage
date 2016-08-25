r"""
Common Asymptotic Expansions

.. WARNING::

    As this code is experimental, a warning is thrown when an
    asymptotic ring (or an associated structure) is created for the
    first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: AsymptoticRing(growth_group='z^ZZ * log(z)^QQ', coefficient_ring=ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Asymptotic Ring <z^ZZ * log(z)^QQ> over Integer Ring


Asymptotic expansions in SageMath can be built through the
``asymptotic_expansions`` object. It contains generators for common
asymptotic expressions. For example,
::

    sage: asymptotic_expansions.Stirling('n', precision=5)
    sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
    1/12*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
    1/288*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
    O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))

generates the first 5 summands of Stirling's approximation formula for
factorials.

To construct an asymptotic expression manually, you can use the class
:class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticRing`. See
:doc:`asymptotic ring <asymptotic_ring>` for more details and a lot of
examples.


**Asymptotic Expansions**

.. list-table::
   :class: contentstable
   :widths: 4 12
   :header-rows: 0

   * - :meth:`~AsymptoticExpansionGenerators.HarmonicNumber`
     - harmonic numbers

   * - :meth:`~AsymptoticExpansionGenerators.Stirling`
     - Stirling's approximation formula for factorials

   * - :meth:`~AsymptoticExpansionGenerators.log_Stirling`
     - the logarithm of Stirling's approximation formula for factorials

   * - :meth:`~AsymptoticExpansionGenerators.Binomial_kn_over_n`
     - an asymptotic expansion of the binomial coefficient

   * - :meth:`~AsymptoticExpansionGenerators.SingularityAnalysis`
     - an asymptotic expansion obtained by singularity analysis


AUTHORS:

- Daniel Krenn (2015)
- Clemens Heuberger (2016)
- Benjamin Hackl (2016)


ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.


Classes and Methods
===================
"""

#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
# Copyright (C) 2016 Clemens Heuberger <clemens.heuberger@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from sage.structure.sage_object import SageObject


class AsymptoticExpansionGenerators(SageObject):
    r"""
    A collection of constructors for several common asymptotic expansions.

    A list of all asymptotic expansions in this database is available via tab
    completion. Type "``asymptotic_expansions.``" and then hit tab to see which
    expansions are available.

    The asymptotic expansions currently in this class include:

    - :meth:`~HarmonicNumber`
    - :meth:`~Stirling`
    - :meth:`~log_Stirling`
    - :meth:`~Binomial_kn_over_n`
    - :meth:`~SingularityAnalysis`
    """

    @staticmethod
    def Stirling(var, precision=None, skip_constant_factor=False):
        r"""
        Return Stirling's approximation formula for factorials.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- (default: ``None``) an integer `\ge 3`. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_factor`` -- (default: ``False``) a
          boolean. If set, then the constant factor `\sqrt{2\pi}` is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.Stirling('n', precision=5)
            sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
            1/12*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
            1/288*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
            O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: _.parent()
            Asymptotic Ring <(e^(n*log(n)))^QQ * (e^n)^QQ * n^QQ * log(n)^QQ>
            over Symbolic Constants Subring

        .. SEEALSO::

            :meth:`log_Stirling`,
            :meth:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticExpansion.factorial`.

        TESTS::

            sage: expansion = asymptotic_expansions.Stirling('n', precision=5)
            sage: n = expansion.parent().gen()
            sage: expansion.compare_with_values(n, lambda x: x.factorial(), [5, 10, 20])  # rel tol 1e-6
            [(5, 0.00675841118?), (10, 0.0067589306?), (20, 0.006744925?)]
            sage: asymptotic_expansions.Stirling('n', precision=5,
            ....:                                skip_constant_factor=True)
            e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
            1/12*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
            1/288*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
            O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: _.parent()
            Asymptotic Ring <(e^(n*log(n)))^QQ * (e^n)^QQ * n^QQ * log(n)^QQ>
            over Rational Field
            sage: asymptotic_expansions.Stirling('m', precision=4)
            sqrt(2)*sqrt(pi)*e^(m*log(m))*(e^m)^(-1)*m^(1/2) +
            O(e^(m*log(m))*(e^m)^(-1)*m^(-1/2))
            sage: asymptotic_expansions.Stirling('m', precision=3)
            O(e^(m*log(m))*(e^m)^(-1)*m^(1/2))
            sage: asymptotic_expansions.Stirling('m', precision=2)
            Traceback (most recent call last):
            ...
            ValueError: precision must be at least 3
        """
        if precision < 3:
            raise ValueError("precision must be at least 3")
        log_Stirling = AsymptoticExpansionGenerators.log_Stirling(
            var, precision=precision, skip_constant_summand=True)

        P = log_Stirling.parent().change_parameter(
            growth_group='(e^({n}*log({n})))^QQ * (e^{n})^QQ * {n}^QQ * log({n})^QQ'.format(n=var))
        from sage.functions.log import exp
        result = exp(P(log_Stirling))

        if not skip_constant_factor:
            from sage.symbolic.ring import SR
            SCR = SR.subring(no_variables=True)
            result *= (2*SCR('pi')).sqrt()

        return result


    @staticmethod
    def log_Stirling(var, precision=None, skip_constant_summand=False):
        r"""
        Return the logarithm of Stirling's approximation formula
        for factorials.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_summand`` -- (default: ``False``) a
          boolean. If set, then the constant summand `\log(2\pi)/2` is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.log_Stirling('n', precision=7)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            - 1/360*n^(-3) + 1/1260*n^(-5) + O(n^(-7))

        .. SEEALSO::

            :meth:`Stirling`,
            :meth:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticExpansion.factorial`.

        TESTS::

            sage: expansion = asymptotic_expansions.log_Stirling('n', precision=7)
            sage: n = expansion.parent().gen()
            sage: expansion.compare_with_values(n, lambda x: x.factorial().log(), [5, 10, 20])  # rel tol 1e-6
            [(5, 0.000564287?), (10, 0.0005870?), (20, 0.0006?)]
            sage: asymptotic_expansions.log_Stirling('n')
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            - 1/360*n^(-3) + 1/1260*n^(-5) - 1/1680*n^(-7) + 1/1188*n^(-9)
            - 691/360360*n^(-11) + 1/156*n^(-13) - 3617/122400*n^(-15)
            + 43867/244188*n^(-17) - 174611/125400*n^(-19) + 77683/5796*n^(-21)
            - 236364091/1506960*n^(-23) + 657931/300*n^(-25)
            - 3392780147/93960*n^(-27) + 1723168255201/2492028*n^(-29)
            - 7709321041217/505920*n^(-31) + O(n^(-33))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Symbolic Constants Subring

        ::

            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=7, skip_constant_summand=True)
            n*log(n) - n + 1/2*log(n) + 1/12*n^(-1) - 1/360*n^(-3) +
            1/1260*n^(-5) + O(n^(-7))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Rational Field
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=0)
            O(n*log(n))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=1)
            n*log(n) + O(n)
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=2)
            n*log(n) - n + O(log(n))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=3)
            n*log(n) - n + 1/2*log(n) + O(1)
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=4)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + O(n^(-1))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=5)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            + O(n^(-3))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'm', precision=7, skip_constant_summand=True)
            m*log(m) - m + 1/2*log(m) + 1/12*m^(-1) - 1/360*m^(-3) +
            1/1260*m^(-5) + O(m^(-7))
        """
        if not skip_constant_summand:
            from sage.symbolic.ring import SR
            coefficient_ring = SR.subring(no_variables=True)
        else:
            from sage.rings.rational_field import QQ
            coefficient_ring = QQ

        from .asymptotic_ring import AsymptoticRing
        A = AsymptoticRing(growth_group='{n}^ZZ * log({n})^ZZ'.format(n=var),
                           coefficient_ring=coefficient_ring)
        n = A.gen()

        if precision is None:
            precision = AsymptoticRing.__default_prec__

        from sage.functions.log import log
        result = A.zero()
        if precision >= 1:
            result += n * log(n)
        if precision >= 2:
            result += -n
        if precision >= 3:
            result += log(n) / 2
        if precision >= 4 and not skip_constant_summand:
            result += log(2*coefficient_ring('pi')) / 2

        result += AsymptoticExpansionGenerators._log_StirlingNegativePowers_(
            var, precision - 4)

        if precision < 1:
            result += (n * log(n)).O()
        elif precision == 1:
            result += n.O()
        elif precision == 2:
            result += log(n).O()
        elif precision == 3:
            result += A(1).O()

        return result


    @staticmethod
    def _log_StirlingNegativePowers_(var, precision):
        r"""
        Helper function to calculate the logarithm of Stirling's approximation
        formula from the negative powers of ``var`` on, i.e., it skips the
        summands `n \log n - n + (\log n)/2 + \log(2\pi)/2`.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- an integer specifying the number of exact summands.
          If this is negative, then the result is `0`.

        OUTPUT:

        An asymptotic expansion.

        TESTS::

            sage: asymptotic_expansions._log_StirlingNegativePowers_(
            ....:     'm', precision=-1)
            0
            sage: asymptotic_expansions._log_StirlingNegativePowers_(
            ....:     'm', precision=0)
            O(m^(-1))
            sage: asymptotic_expansions._log_StirlingNegativePowers_(
            ....:     'm', precision=3)
            1/12*m^(-1) - 1/360*m^(-3) + 1/1260*m^(-5) + O(m^(-7))
            sage: _.parent()
            Asymptotic Ring <m^ZZ> over Rational Field
        """
        from .asymptotic_ring import AsymptoticRing
        from sage.rings.rational_field import QQ

        A = AsymptoticRing(growth_group='{n}^ZZ'.format(n=var),
                           coefficient_ring=QQ)
        if precision < 0:
            return A.zero()
        n = A.gen()

        from sage.arith.all import bernoulli
        from sage.arith.srange import srange

        result = sum((bernoulli(k) / k / (k-1) / n**(k-1)
                      for k in srange(2, 2*precision + 2, 2)),
                     A.zero())
        return result + (1 / n**(2*precision + 1)).O()


    @staticmethod
    def HarmonicNumber(var, precision=None, skip_constant_summand=False):
        r"""
        Return the asymptotic expansion of a harmonic number.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_summand`` -- (default: ``False``) a
          boolean. If set, then the constant summand ``euler_gamma`` is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.HarmonicNumber('n', precision=5)
            log(n) + euler_gamma + 1/2*n^(-1) - 1/12*n^(-2) + 1/120*n^(-4) + O(n^(-6))

        TESTS::

            sage: ex = asymptotic_expansions.HarmonicNumber('n', precision=5)
            sage: n = ex.parent().gen()
            sage: ex.compare_with_values(n,                      # rel tol 1e-6
            ....:      lambda x: sum(1/k for k in srange(1, x+1)), [5, 10, 20])
            [(5, 0.0038125360?), (10, 0.00392733?), (20, 0.0039579?)]
            sage: asymptotic_expansions.HarmonicNumber('n')
            log(n) + euler_gamma + 1/2*n^(-1) - 1/12*n^(-2) + 1/120*n^(-4)
            - 1/252*n^(-6) + 1/240*n^(-8) - 1/132*n^(-10)
            + 691/32760*n^(-12) - 1/12*n^(-14) + 3617/8160*n^(-16)
            - 43867/14364*n^(-18) + 174611/6600*n^(-20) - 77683/276*n^(-22)
            + 236364091/65520*n^(-24) - 657931/12*n^(-26)
            + 3392780147/3480*n^(-28) - 1723168255201/85932*n^(-30)
            + 7709321041217/16320*n^(-32)
            - 151628697551/12*n^(-34) + O(n^(-36))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Symbolic Constants Subring

        ::

            sage: asymptotic_expansions.HarmonicNumber(
            ....:     'n', precision=5, skip_constant_summand=True)
            log(n) + 1/2*n^(-1) - 1/12*n^(-2) + 1/120*n^(-4) + O(n^(-6))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Rational Field
            sage: for p in range(5):
            ....:     print(asymptotic_expansions.HarmonicNumber(
            ....:         'n', precision=p))
            O(log(n))
            log(n) + O(1)
            log(n) + euler_gamma + O(n^(-1))
            log(n) + euler_gamma + 1/2*n^(-1) + O(n^(-2))
            log(n) + euler_gamma + 1/2*n^(-1) - 1/12*n^(-2) + O(n^(-4))
            sage: asymptotic_expansions.HarmonicNumber('m', precision=5)
            log(m) + euler_gamma + 1/2*m^(-1) - 1/12*m^(-2) + 1/120*m^(-4) + O(m^(-6))
        """
        if not skip_constant_summand:
            from sage.symbolic.ring import SR
            coefficient_ring = SR.subring(no_variables=True)
        else:
            from sage.rings.rational_field import QQ
            coefficient_ring = QQ

        from .asymptotic_ring import AsymptoticRing
        A = AsymptoticRing(growth_group='{n}^ZZ * log({n})^ZZ'.format(n=var),
                           coefficient_ring=coefficient_ring)
        n = A.gen()

        if precision is None:
            precision = A.default_prec

        from sage.functions.log import log
        result = A.zero()
        if precision >= 1:
            result += log(n)
        if precision >= 2 and not skip_constant_summand:
            from sage.symbolic.constants import euler_gamma
            result += coefficient_ring(euler_gamma)
        if precision >= 3:
            result += 1 / (2 * n)

        from sage.arith.srange import srange
        from sage.arith.all import bernoulli
        for k in srange(2, 2*precision - 4, 2):
            result += -bernoulli(k) / k / n**k

        if precision < 1:
            result += (log(n)).O()
        elif precision == 1:
            result += A(1).O()
        elif precision == 2:
            result += (1 / n).O()
        else:
            result += (1 / n**(2*precision - 4)).O()

        return result


    @staticmethod
    def Binomial_kn_over_n(var, k, precision=None, skip_constant_factor=False):
        r"""
        Return the asymptotic expansion of the binomial coefficient
        `kn` choose `n`.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``k`` -- a number or symbolic constant.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_factor`` -- (default: ``False``) a
          boolean. If set, then the constant factor `\sqrt{k/(2\pi(k-1))}`
          is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=2, precision=3)
            1/sqrt(pi)*4^n*n^(-1/2)
            - 1/8/sqrt(pi)*4^n*n^(-3/2)
            + 1/128/sqrt(pi)*4^n*n^(-5/2)
            + O(4^n*n^(-7/2))
            sage: _.parent()
            Asymptotic Ring <QQ^n * n^QQ> over Symbolic Constants Subring

        ::

            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=3, precision=3)
            1/2*sqrt(3)/sqrt(pi)*(27/4)^n*n^(-1/2)
            - 7/144*sqrt(3)/sqrt(pi)*(27/4)^n*n^(-3/2)
            + 49/20736*sqrt(3)/sqrt(pi)*(27/4)^n*n^(-5/2)
            + O((27/4)^n*n^(-7/2))

        ::

            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=7/5, precision=3)
            1/2*sqrt(7)/sqrt(pi)*(7/10*7^(2/5)*2^(3/5))^n*n^(-1/2)
            - 13/112*sqrt(7)/sqrt(pi)*(7/10*7^(2/5)*2^(3/5))^n*n^(-3/2)
            + 169/12544*sqrt(7)/sqrt(pi)*(7/10*7^(2/5)*2^(3/5))^n*n^(-5/2)
            + O((7/10*7^(2/5)*2^(3/5))^n*n^(-7/2))
            sage: _.parent()
            Asymptotic Ring <(Symbolic Constants Subring)^n * n^QQ>
            over Symbolic Constants Subring

        TESTS::

            sage: expansion = asymptotic_expansions.Binomial_kn_over_n('n', k=7/5, precision=3)
            sage: n = expansion.parent().gen()
            sage: expansion.compare_with_values(n, lambda x: binomial(7/5*x, x), [5, 10, 20])  # rel tol 1e-6
            [(5, -0.0287383845047?), (10, -0.030845971026?), (20, -0.03162833549?)]
            sage: asymptotic_expansions.Binomial_kn_over_n(
            ....:     'n', k=5, precision=3, skip_constant_factor=True)
            (3125/256)^n*n^(-1/2)
            - 7/80*(3125/256)^n*n^(-3/2)
            + 49/12800*(3125/256)^n*n^(-5/2)
            + O((3125/256)^n*n^(-7/2))
            sage: _.parent()
            Asymptotic Ring <QQ^n * n^QQ> over Rational Field
            sage: asymptotic_expansions.Binomial_kn_over_n(
            ....:     'n', k=4, precision=1, skip_constant_factor=True)
            (256/27)^n*n^(-1/2) + O((256/27)^n*n^(-3/2))

        ::

            sage: S = asymptotic_expansions.Stirling('n', precision=5)
            sage: n = S.parent().gen()
            sage: all(  # long time
            ....:     SR(asymptotic_expansions.Binomial_kn_over_n(
            ....:         'n', k=k, precision=3)).canonicalize_radical() ==
            ....:     SR(S.subs(n=k*n) / (S.subs(n=(k-1)*n) * S)).canonicalize_radical()
            ....:     for k in [2, 3, 4])
            True
        """
        from sage.symbolic.ring import SR
        SCR = SR.subring(no_variables=True)
        try:
            SCR.coerce(k)
        except TypeError as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                TypeError('Cannot use k={}.'.format(k)), e)

        S = AsymptoticExpansionGenerators._log_StirlingNegativePowers_(
                var, precision=max(precision - 2,0))
        n = S.parent().gen()
        result = (S.subs(n=k*n) - S.subs(n=(k-1)*n) - S).exp()

        from sage.rings.rational_field import QQ

        P = S.parent().change_parameter(
                growth_group='QQ^{n} * {n}^QQ'.format(n=var),
                coefficient_ring=QQ)
        n = P.gen()

        b = k**k / (k-1)**(k-1)
        if b.parent() is SR:
            b = SCR(b).canonicalize_radical()
        result *= n.rpow(b)
        result *= n**(-QQ(1)/QQ(2))
        if not skip_constant_factor:
            result *= (k/((k-1)*2*SCR('pi'))).sqrt()

        return result


    @staticmethod
    def SingularityAnalysis(var, zeta=1, alpha=0, beta=0, delta=0,
                            precision=None, normalized=True):
        r"""
        Return the asymptotic expansion of the coefficients of
        an power series with specified pole and logarithmic singularity.

        More precisely, this extracts the `n`-th coefficient

        .. MATH::

            [z^n] \left(\frac{1}{1-z/\zeta}\right)^\alpha
            \left(\frac{1}{z/\zeta} \log \frac{1}{1-z/\zeta}\right)^\beta
            \left(\frac{1}{z/\zeta} \log
            \left(\frac{1}{z/\zeta} \log \frac{1}{1-z/\zeta}\right)\right)^\delta

        (if ``normalized=True``, the default) or

        .. MATH::

            [z^n] \left(\frac{1}{1-z/\zeta}\right)^\alpha
            \left(\log \frac{1}{1-z/\zeta}\right)^\beta
            \left(\log
            \left(\frac{1}{z/\zeta} \log \frac{1}{1-z/\zeta}\right)\right)^\delta

        (if ``normalized=False``).

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``zeta`` -- (default: `1`) the location of the singularity.

        - ``alpha`` -- (default: `0`) the pole order of the singularty.

        - ``beta`` -- (default: `0`) the order of the logarithmic singularity.

        - ``delta`` -- (default: `0`) the order of the log-log singularity.
          Not yet implemented for ``delta != 0``.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``normalized`` -- (default: ``True``) a boolean, see above.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.SingularityAnalysis('n', alpha=1)
            1
            sage: asymptotic_expansions.SingularityAnalysis('n', alpha=2)
            n + 1
            sage: asymptotic_expansions.SingularityAnalysis('n', alpha=3)
            1/2*n^2 + 3/2*n + 1
            sage: _.parent()
            Asymptotic Ring <n^ZZ> over Rational Field

        ::

            sage: asymptotic_expansions.SingularityAnalysis('n', alpha=-3/2,
            ....:     precision=3)
            3/4/sqrt(pi)*n^(-5/2)
            + 45/32/sqrt(pi)*n^(-7/2)
            + 1155/512/sqrt(pi)*n^(-9/2)
            + O(n^(-11/2))
            sage: asymptotic_expansions.SingularityAnalysis('n', alpha=-1/2,
            ....:     precision=3)
            -1/2/sqrt(pi)*n^(-3/2)
            - 3/16/sqrt(pi)*n^(-5/2)
            - 25/256/sqrt(pi)*n^(-7/2)
            + O(n^(-9/2))
            sage: asymptotic_expansions.SingularityAnalysis('n', alpha=1/2,
            ....:     precision=4)
            1/sqrt(pi)*n^(-1/2)
            - 1/8/sqrt(pi)*n^(-3/2)
            + 1/128/sqrt(pi)*n^(-5/2)
            + 5/1024/sqrt(pi)*n^(-7/2)
            + O(n^(-9/2))
            sage: _.parent()
            Asymptotic Ring <n^QQ> over Symbolic Constants Subring

        ::

            sage: S = SR.subring(rejecting_variables=('n',))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=S.var('a'),
            ....:     precision=4).map_coefficients(lambda c: c.factor())
            1/gamma(a)*n^(a - 1)
            + (1/2*(a - 1)*a/gamma(a))*n^(a - 2)
            + (1/24*(3*a - 1)*(a - 1)*(a - 2)*a/gamma(a))*n^(a - 3)
            + (1/48*(a - 1)^2*(a - 2)*(a - 3)*a^2/gamma(a))*n^(a - 4)
            + O(n^(a - 5))
            sage: _.parent()
            Asymptotic Ring <n^(Symbolic Subring rejecting the variable n)>
            over Symbolic Subring rejecting the variable n

        ::

            sage: ae = asymptotic_expansions.SingularityAnalysis('n',
            ....:          alpha=1/2, beta=1, precision=4); ae
            1/sqrt(pi)*n^(-1/2)*log(n) + ((euler_gamma + 2*log(2))/sqrt(pi))*n^(-1/2)
            - 5/8/sqrt(pi)*n^(-3/2)*log(n) + (1/8*(3*euler_gamma + 6*log(2) - 8)/sqrt(pi)
            - (euler_gamma + 2*log(2) - 2)/sqrt(pi))*n^(-3/2) + O(n^(-5/2)*log(n))
            sage: n = ae.parent().gen()
            sage: ae.subs(n=n-1).map_coefficients(lambda x: x.canonicalize_radical())
            1/sqrt(pi)*n^(-1/2)*log(n)
            + ((euler_gamma + 2*log(2))/sqrt(pi))*n^(-1/2)
            - 1/8/sqrt(pi)*n^(-3/2)*log(n)
            + (-1/8*(euler_gamma + 2*log(2))/sqrt(pi))*n^(-3/2)
            + O(n^(-5/2)*log(n))

        ::

            sage: asymptotic_expansions.SingularityAnalysis('n',
            ....:     alpha=1, beta=1/2, precision=4)
            log(n)^(1/2) + 1/2*euler_gamma*log(n)^(-1/2)
            + (-1/8*euler_gamma^2 + 1/48*pi^2)*log(n)^(-3/2)
            + (1/16*euler_gamma^3 - 1/32*euler_gamma*pi^2 + 1/8*zeta(3))*log(n)^(-5/2)
            + O(log(n)^(-7/2))

        ::

            sage: ae = asymptotic_expansions.SingularityAnalysis('n',
            ....:     alpha=0, beta=2, precision=14)
            sage: n = ae.parent().gen()
            sage: ae.subs(n=n-2)
            2*n^(-1)*log(n) + 2*euler_gamma*n^(-1) - n^(-2) - 1/6*n^(-3) + O(n^(-5))

        ::

            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=-1/2, beta=1, precision=2, normalized=False)
            -1/2/sqrt(pi)*n^(-3/2)*log(n)
            + (-1/2*(euler_gamma + 2*log(2) - 2)/sqrt(pi))*n^(-3/2)
            + O(n^(-5/2)*log(n))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1/2, alpha=0, beta=1, precision=3, normalized=False)
            2^n*n^(-1) + O(2^n*n^(-2))


        ALGORITHM:

        See [FS2009]_ together with the
        `errata list <http://algo.inria.fr/flajolet/Publications/AnaCombi/errata.pdf>`_.

        REFERENCES:

        .. [FS2009] Philippe Flajolet and Robert Sedgewick,
           `Analytic combinatorics <http://algo.inria.fr/flajolet/Publications/AnaCombi/book.pdf>`_.
           Cambridge University Press, Cambridge, 2009.

        TESTS::

            sage: ex = asymptotic_expansions.SingularityAnalysis('n', alpha=-1/2,
            ....:     precision=4)
            sage: n = ex.parent().gen()
            sage: coefficients = ((1-x)^(1/2)).series(
            ....:     x, 21).truncate().coefficients(x, sparse=False)
            sage: ex.compare_with_values(n,    # rel tol 1e-6
            ....:     lambda k: coefficients[k], [5, 10, 20])
            [(5, 0.015778873294?), (10, 0.01498952777?), (20, 0.0146264622?)]
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=3, precision=2)
            1/2*n^2 + 3/2*n + O(1)
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=3, precision=3)
            1/2*n^2 + 3/2*n + 1
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=3, precision=4)
            1/2*n^2 + 3/2*n + 1

        ::

            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=0)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: The error term in the result is O(0)
            which means 0 for sufficiently large n.
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=-1)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: The error term in the result is O(0)
            which means 0 for sufficiently large n.

        ::

            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'm', alpha=-1/2, precision=3)
            -1/2/sqrt(pi)*m^(-3/2)
            - 3/16/sqrt(pi)*m^(-5/2)
            - 25/256/sqrt(pi)*m^(-7/2)
            + O(m^(-9/2))
            sage: _.parent()
            Asymptotic Ring <m^QQ> over Symbolic Constants Subring

        Location of the singularity::

            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=1, zeta=2, precision=3)
            (1/2)^n
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=1, zeta=1/2, precision=3)
            2^n
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=1, zeta=CyclotomicField(3).gen(),
            ....:      precision=3)
            (-zeta3 - 1)^n
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=4, zeta=2, precision=3)
            1/6*(1/2)^n*n^3 + (1/2)^n*n^2 + 11/6*(1/2)^n*n + O((1/2)^n)
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=-1, zeta=2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: The error term in the result is O(0)
            which means 0 for sufficiently large n.
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=1/2, zeta=2, precision=3)
            1/sqrt(pi)*(1/2)^n*n^(-1/2) - 1/8/sqrt(pi)*(1/2)^n*n^(-3/2)
            + 1/128/sqrt(pi)*(1/2)^n*n^(-5/2) + O((1/2)^n*n^(-7/2))

        The following tests correspond to Table VI.5 in [FS2009]_. ::

            sage: A.<n> = AsymptoticRing('n^QQ * log(n)^QQ', QQ)
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=-1/2, beta=1, precision=2,
            ....:     normalized=False) * (- sqrt(pi*n^3))
            1/2*log(n) + 1/2*euler_gamma + log(2) - 1 + O(n^(-1)*log(n))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=0, beta=1, precision=3,
            ....:     normalized=False)
            n^(-1) + O(n^(-2))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=0, beta=2,  precision=14,
            ....:     normalized=False) * n
            2*log(n) + 2*euler_gamma - n^(-1) - 1/6*n^(-2) +  O(n^(-4))
            sage: (asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=1/2, beta=1, precision=4,
            ....:     normalized=False) * sqrt(pi*n)).\
            ....:     map_coefficients(lambda x: x.expand())
            log(n) + euler_gamma + 2*log(2) - 1/8*n^(-1)*log(n) +
            (-1/8*euler_gamma - 1/4*log(2))*n^(-1) + O(n^(-2)*log(n))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=1, beta=1, precision=13,
            ....:     normalized=False)
            log(n) + euler_gamma + 1/2*n^(-1) - 1/12*n^(-2) + 1/120*n^(-4)
            + O(n^(-6))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=1, beta=2, precision=4,
            ....:     normalized=False)
            log(n)^2 + 2*euler_gamma*log(n) + euler_gamma^2 - 1/6*pi^2
            + O(n^(-1)*log(n))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=3/2, beta=1, precision=3,
            ....:     normalized=False) * sqrt(pi/n)
            2*log(n) + 2*euler_gamma + 4*log(2) - 4 + 3/4*n^(-1)*log(n)
            + O(n^(-1))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=2, beta=1, precision=5,
            ....:     normalized=False)
            n*log(n) + (euler_gamma - 1)*n + log(n) + euler_gamma + 1/2
            + O(n^(-1))
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, alpha=2, beta=2, precision=4,
            ....:     normalized=False) / n
            log(n)^2 + (2*euler_gamma - 2)*log(n)
            - 2*euler_gamma + euler_gamma^2 - 1/6*pi^2 + 2
            + n^(-1)*log(n)^2 + O(n^(-1)*log(n))

        Be aware that the last result does *not* coincide with [FS2009]_,
        they do have a different error term.

        Checking parameters::

            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, 1, 1/2, precision=0, normalized=False)
            Traceback (most recent call last):
            ...
            ValueError: beta and delta must be integers
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', 1, 1, 1, 1/2, normalized=False)
            Traceback (most recent call last):
            ...
            ValueError: beta and delta must be integers

        ::

            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=0, beta=0, delta=1, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented for delta!=0
        """
        from itertools import islice, count
        from .asymptotic_ring import AsymptoticRing
        from .growth_group import ExponentialGrowthGroup, \
                MonomialGrowthGroup
        from sage.arith.all import falling_factorial
        from sage.categories.cartesian_product import cartesian_product
        from sage.functions.other import binomial, gamma
        from sage.calculus.calculus import limit
        from sage.misc.cachefunc import cached_function
        from sage.arith.srange import srange
        from sage.rings.rational_field import QQ
        from sage.rings.integer_ring import ZZ
        from sage.symbolic.ring import SR

        SCR = SR.subring(no_variables=True)
        s = SR('s')
        iga = 1/gamma(alpha)
        if iga.parent() is SR:
            try:
                iga = SCR(iga)
            except TypeError:
                pass

        coefficient_ring = iga.parent()
        if beta != 0:
            coefficient_ring = SCR

        @cached_function
        def inverse_gamma_derivative(shift, r):
            """
            Return value of `r`-th derivative of 1/Gamma
            at alpha-shift.
            """
            if r == 0:
                result = iga*falling_factorial(alpha-1, shift)
            else:
                result = limit((1/gamma(s)).diff(s, r), s=alpha-shift)

            try:
                return coefficient_ring(result)
            except TypeError:
                return result

        if isinstance(alpha, int):
            alpha = ZZ(alpha)
        if isinstance(beta, int):
            beta = ZZ(beta)
        if isinstance(delta, int):
            delta = ZZ(delta)

        if precision is None:
            precision = AsymptoticRing.__default_prec__


        if not normalized and not (beta in ZZ and delta in ZZ):
            raise ValueError("beta and delta must be integers")
        if delta != 0:
            raise NotImplementedError("not implemented for delta!=0")

        groups = []
        if zeta != 1:
            groups.append(ExponentialGrowthGroup((1/zeta).parent(), var))

        groups.append(MonomialGrowthGroup(alpha.parent(), var))
        if beta != 0:
            groups.append(MonomialGrowthGroup(beta.parent(), 'log({})'.format(var)))

        group = cartesian_product(groups)
        A = AsymptoticRing(growth_group=group, coefficient_ring=coefficient_ring,
                           default_prec=precision)
        n = A.gen()

        if zeta == 1:
            exponential_factor = 1
        else:
            exponential_factor = n.rpow(1/zeta)

        if beta in ZZ and beta >= 0:
            it = ((k, r)
                  for k in count()
                  for r in srange(beta+1))
            k_max = precision
        else:
            it = ((0, r)
                  for r in count())
            k_max = 0

        if beta != 0:
            log_n = n.log()
        else:
            # avoid construction of log(n)
            # because it does not exist in growth group.
            log_n = 1

        it = reversed(list(islice(it, precision+1)))
        if normalized:
            beta_denominator = beta
        else:
            beta_denominator = 0
        L = _sa_coefficients_lambda_(max(1, k_max), beta=beta_denominator)
        (k, r) = next(it)
        result = (n**(-k) * log_n**(-r)).O()

        if alpha in ZZ and beta == 0:
            if alpha > 0 and alpha <= precision:
                result = A(0)
            elif alpha <= 0 and precision > 0:
                from .misc import NotImplementedOZero
                raise NotImplementedOZero(A)

        for (k, r) in it:
            result += binomial(beta, r) * \
                sum(L[(k, ell)] * (-1)**ell *
                    inverse_gamma_derivative(ell, r)
                    for ell in srange(k, 2*k+1)
                    if (k, ell) in L) * \
                n**(-k) * log_n**(-r)

        result *= exponential_factor * n**(alpha-1) * log_n**beta

        return result


def _sa_coefficients_lambda_(K, beta=0):
    r"""
    Return the coefficients `\lambda_{k, \ell}(\beta)` used in singularity analysis.

    INPUT:

    - ``K`` -- an integer.

    - ``beta`` -- (default: `0`) the order of the logarithmic
      singularity.

    OUTPUT:

    A dictionary mapping pairs of indices to rationals.

    .. SEEALSO::

        :meth:`~AsymptoticExpansionGenerators.SingularityAnalysis`

    TESTS::

        sage: from sage.rings.asymptotic.asymptotic_expansion_generators \
        ....:     import _sa_coefficients_lambda_
        sage: _sa_coefficients_lambda_(3)
        {(0, 0): 1,
         (1, 1): -1,
         (1, 2): 1/2,
         (2, 2): 1,
         (2, 3): -5/6,
         (2, 4): 1/8,
         (3, 3): -1,
         (3, 4): 13/12,
         (4, 4): 1}
        sage: _sa_coefficients_lambda_(3, beta=1)
        {(0, 0): 1,
         (1, 1): -2,
         (1, 2): 1/2,
         (2, 2): 3,
         (2, 3): -4/3,
         (2, 4): 1/8,
         (3, 3): -4,
         (3, 4): 29/12,
         (4, 4): 5}
    """
    from sage.rings.laurent_series_ring import LaurentSeriesRing
    from sage.rings.power_series_ring import PowerSeriesRing
    from sage.rings.rational_field import QQ

    V = LaurentSeriesRing(QQ, names='v', default_prec=K)
    v = V.gen()
    T = PowerSeriesRing(V, names='t', default_prec=2*K-1)
    t = T.gen()

    S = (t - (1+1/v+beta) * (1+v*t).log()).exp()
    return dict(((k + L.valuation(), ell), c)
                for ell, L in enumerate(S.list())
                for k, c in enumerate(L.list()))


# Easy access to the asymptotic expansions generators from the command line:
asymptotic_expansions = AsymptoticExpansionGenerators()
r"""
A collection of several common asymptotic expansions.

This is an instance of :class:`AsymptoticExpansionGenerators` whose documentation
provides more details.
"""
