r"""
Common Asymptotic Expansions

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

   * - :meth:`~AsymptoticExpansionGenerators.ImplicitExpansion`
     - the singular expansion of a function `y(z)` satisfying `y(z) = z \Phi(y(z))`

   * - :meth:`~AsymptoticExpansionGenerators.ImplicitExpansionPeriodicPart`
     - the singular expansion of the periodic part of a function `y(z)` satisfying `y(z) = z\Phi(y(z))`

   * - :meth:`~AsymptoticExpansionGenerators.InverseFunctionAnalysis`
     - coefficient growth of a function `y(z)` defined implicitly by `y(z) = z \Phi(y(z))`


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

from sage.misc.superseded import experimental
from sage.structure.sage_object import SageObject
from sage.misc.defaults import series_precision


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
    - :meth:`~ImplicitExpansion`
    - :meth:`~ImplicitExpansionPeriodicPart`
    - :meth:`~InverseFunctionAnalysis`
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

        Check that :trac:`20066` is resolved::

            sage: set_series_precision(5)
            sage: asymptotic_expansions.Stirling('n')
            sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(1/2) + 
            ... + O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: set_series_precision(20)  # restore series precision default
        """
        if precision is None:
            precision = series_precision()

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
            precision = series_precision()

        log = A.locals()['log']
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
            precision = series_precision()

        log = A.locals()['log']
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

        Check that :trac:`20066` is resolved::

            sage: set_series_precision(3)
            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=2)
            1/sqrt(pi)*4^n*n^(-1/2) - 1/8/sqrt(pi)*4^n*n^(-3/2) + ... + O(4^n*n^(-7/2))
            sage: set_series_precision(20)  # restore series precision default
        """
        from sage.symbolic.ring import SR
        SCR = SR.subring(no_variables=True)
        try:
            SCR.coerce(k)
        except TypeError as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                TypeError('Cannot use k={}.'.format(k)), e)

        if precision is None:
            precision = series_precision()

        S = AsymptoticExpansionGenerators._log_StirlingNegativePowers_(
                var, precision=max(precision - 2,0))
        n = S.parent().gen()
        result = (S.subs(n=k*n) - S.subs(n=(k-1)*n) - S).exp()

        from sage.rings.rational_field import QQ

        P = S.parent().change_parameter(
                growth_group='(QQ_+)^{n} * {n}^QQ'.format(n=var),
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

        - ``alpha`` -- (default: `0`) the pole order of the singularity.

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

        See [FS2009]_.


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
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large n.
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=-1)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large n.

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
            (e^(I*arg(-zeta3 - 1)))^n
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=4, zeta=2, precision=3)
            1/6*(1/2)^n*n^3 + (1/2)^n*n^2 + 11/6*(1/2)^n*n + O((1/2)^n)
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', alpha=-1, zeta=2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large n.
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

        ::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: Signs = SignGroup()
            sage: asymptotic_expansions.SingularityAnalysis(
            ....:     'n', Signs(-1), alpha=2, beta=1, precision=5,
            ....:     normalized=False)
            n*log(n)*(-1)^n + (euler_gamma - 1)*n*(-1)^n + log(n)*(-1)^n
            + (euler_gamma + 1/2)*(-1)^n + O(n^(-1)*(-1)^n)
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ * Signs^n> over Symbolic Constants Subring
        """
        from itertools import islice, count
        from .asymptotic_ring import AsymptoticRing
        from .growth_group import ExponentialGrowthGroup, \
                MonomialGrowthGroup, GenericNonGrowthGroup
        from sage.arith.all import falling_factorial
        from sage.categories.cartesian_product import cartesian_product
        from sage.functions.other import binomial
        from sage.functions.gamma import gamma
        from sage.calculus.calculus import limit
        from sage.misc.cachefunc import cached_function
        from sage.arith.srange import srange
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
            precision = series_precision()


        if not normalized and not (beta in ZZ and delta in ZZ):
            raise ValueError("beta and delta must be integers")
        if delta != 0:
            raise NotImplementedError("not implemented for delta!=0")

        groups = []
        non_growth_groups = []
        if zeta != 1:
            E = ExponentialGrowthGroup.factory((~zeta).parent(), var,
                                               return_factors=True)
            for factor in E:
                if isinstance(factor, GenericNonGrowthGroup):
                    non_growth_groups.append(factor)
                else:
                    groups.append(factor)
        groups.append(MonomialGrowthGroup(alpha.parent(), var))
        if beta != 0:
            groups.append(MonomialGrowthGroup(beta.parent(), 'log({})'.format(var)))
        groups.extend(non_growth_groups)
        group = cartesian_product(groups)
        A = AsymptoticRing(growth_group=group, coefficient_ring=coefficient_ring,
                           default_prec=precision)
        n = A.gen()

        if zeta == 1:
            exponential_factor = 1
        else:
            exponential_factor = A(n.rpow(~zeta))

        polynomial_factor = A(n**(alpha-1))

        if beta != 0:
            log_n = n.log()
            logarithmic_factor = log_n**beta
        else:
            # avoid construction of log(n)
            # because it does not exist in growth group.
            log_n = 1
            logarithmic_factor = 1

        if beta in ZZ and beta >= 0:
            it = ((k, r)
                  for k in count()
                  for r in srange(beta+1))
            k_max = precision
        else:
            it = ((0, r)
                  for r in count())
            k_max = 0

        it = reversed(list(islice(it, int(precision) + 1)))
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
                raise NotImplementedOZero(A, exact_part=A.zero())

        for (k, r) in it:
            result += binomial(beta, r) * \
                sum(L[(k, ell)] * (-1)**ell *
                    inverse_gamma_derivative(ell, r)
                    for ell in srange(k, 2*k+1)
                    if (k, ell) in L) * \
                n**(-k) * log_n**(-r)

        result *= exponential_factor * polynomial_factor * logarithmic_factor

        return result


    @staticmethod
    @experimental(20050)
    def ImplicitExpansion(var, phi, tau=None, precision=None):
        r"""
        Return the singular expansion for a function `y(z)` defined
        implicitly by `y(z) = z \Phi(y(z))`.

        The function `\Phi` is assumed to be analytic around `0`. Furthermore,
        `\Phi` is not allowed to be an affine-linear function and we require
        `\Phi(0) \neq 0`.

        Furthermore, it is assumed that there is a unique positive solution `\tau`
        of `\Phi(\tau) - \tau\Phi'(\tau) = 0`.

        All details are given in Chapter VI.7 of [FS2009]_.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``phi`` -- the function `\Phi`. See the extended description for
          assumptions on `\Phi`.

        - ``tau`` -- (default: ``None``) the fundamental constant described
          in the extended description. If ``None``, then `\tau` is determined
          automatically if possible.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.


        OUTPUT:

        An asymptotic expansion.


        .. NOTE::

            In the given case, the radius of convergence of the function of
            interest is known to be `\rho = \tau/\Phi(\tau)`.  Until :trac:`20050`
            is implemented, the variable in the returned asymptotic expansion
            represents a singular element of the form `(1 - z/\rho)^{-1}`,
            for the variable `z\to\rho`.


        EXAMPLES:

        We can, for example, determine the singular expansion of the well-known
        tree function `T` (which satisfies `T(z) = z \exp(T(z))`)::

            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=exp, precision=8)
            doctest:warning
            ...
            FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See http://trac.sagemath.org/20050 for details.
            1 - sqrt(2)*Z^(-1/2) + 2/3*Z^(-1) - 11/36*sqrt(2)*Z^(-3/2) +
            43/135*Z^(-2) - 769/4320*sqrt(2)*Z^(-5/2) + 1768/8505*Z^(-3) + O(Z^(-7/2))

        Another classical example in this context is the generating function `B(z)`
        enumerating binary trees with respect to the number of inner nodes. The
        function satisfies `B(z) = z (1 + 2B(z) + B(z)^2)`, which can also be
        solved explicitly, yielding `B(z) = \frac{1 - \sqrt{1 - 4z}}{2z} - 1`. We
        compare the expansions from both approaches::

            sage: def B(z):
            ....:     return (1 - sqrt(1 - 4*z))/(2*z) - 1
            sage: A.<Z> = AsymptoticRing('Z^QQ', QQ, default_prec=3)
            sage: B((1-1/Z)/4)
            1 - 2*Z^(-1/2) + 2*Z^(-1) - 2*Z^(-3/2) + 2*Z^(-2)
            - 2*Z^(-5/2) + O(Z^(-3))
            sage: asymptotic_expansions.ImplicitExpansion(Z, phi=lambda u: 1 + 2*u + u^2, precision=7)
            1 - 2*Z^(-1/2) + 2*Z^(-1) - 2*Z^(-3/2) + 2*Z^(-2)
            - 2*Z^(-5/2) + O(Z^(-3))

        Neither `\tau` nor `\Phi` have to be known explicitly, they can
        also be passed symbolically::

            sage: tau = var('tau')
            sage: phi = function('phi')
            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=phi, tau=tau, precision=3)  # long time
            tau + (-sqrt(2)*sqrt(-tau*phi(tau)^2/(2*tau*diff(phi(tau), tau)^2
            - tau*phi(tau)*diff(phi(tau), tau, tau)
            - 2*phi(tau)*diff(phi(tau), tau))))*Z^(-1/2) + O(Z^(-1))

        Note that we do not check whether a passed `\tau` actually
        satisfies the requirements. Only the first of the following
        expansions is correct::

            sage: asymptotic_expansions.ImplicitExpansion('Z',
            ....:     phi=lambda u: 1 + 2*u + u^2, precision=5) # correct expansion
            1 - 2*Z^(-1/2) + 2*Z^(-1) - 2*Z^(-3/2) + O(Z^(-2))
            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=lambda u: 1 + 2*u + u^2, tau=2, precision=5)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: symbolic division by zero
            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=lambda u: 1 + 2*u + u^2, tau=3, precision=5)
            3 - 4*I*sqrt(3)*Z^(-1/2) + 6*I*sqrt(3)*Z^(-3/2) + O(Z^(-2))

        .. SEEALSO::

            :meth:`~AsymptoticExpansionGenerators.ImplicitExpansionPeriodicPart`,
            :meth:`~AsymptoticExpansionGenerators.InverseFunctionAnalysis`.

        TESTS::

            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=lambda u: 1 + 42*u, precision=5)
            Traceback (most recent call last):
            ...
            ValueError: The function phi does not satisfy the requirements
            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=lambda u: 42*u + u^2, precision=5)
            Traceback (most recent call last):
            ...
            ValueError: The function phi does not satisfy the requirements
            sage: asymptotic_expansions.ImplicitExpansion('Z', phi=lambda u: 1 + u^2 + u^42, precision=5)
            Traceback (most recent call last):
            ...
            ValueError: Fundamental constant tau could not be determined

        """
        from sage.symbolic.ring import SR
        from sage.rings.rational_field import QQ
        from sage.rings.integer_ring import ZZ
        from sage.rings.asymptotic.asymptotic_ring import AsymptoticRing
        from sage.arith.srange import srange
        y, u = SR('y'), SR('u')
        one_half = QQ(1)/2

        if phi(QQ(0)).is_zero() or phi(u) == phi(0) + u*phi(u).diff(u)(u=0):
            raise ValueError('The function phi does not satisfy the requirements')

        if tau is None:
            tau = _fundamental_constant_implicit_function_(phi=phi)

        def H(y):
            return tau/phi(tau) - y/phi(y)

        A = AsymptoticRing(growth_group='{Z}^QQ'.format(Z=var),
                           coefficient_ring=SR,
                           default_prec=precision)
        if precision is None:
            precision = ZZ(A.default_prec)
        Z = A.gen()

        def ansatz(prec=precision):
            if prec < 1:
                return A(1).O()
            if prec == 1:
                return ((1/Z)**one_half).O()
            return (-(2*tau/phi(tau)/H(y).diff(y, 2)(y=tau)).sqrt() * (1/Z)**one_half
                    + sum(SR("d{}".format(j)) * (1/Z)**(j * one_half) for j in srange(2, prec))
                    + ((1/Z)**(prec * one_half)).O())

        # we compare coefficients between a "single" Z and the
        # following expansion, this allows us to compute the constants d_j
        z = SR('z')
        z_expansion = sum(H(z).diff(z, k)(z=tau)/k.factorial() *
                          ansatz(prec=precision+2-k)**k
                          for k in srange(2, precision)) + ((1/Z)**(precision * one_half)).O()

        solution_dict = dict()
        for k in srange(2, precision-1):
            coef = z_expansion.monomial_coefficient((1/Z)**((k+1) * one_half))
            current_var = SR('d{k}'.format(k=k))
            solution_dict[current_var] = coef.subs(solution_dict).simplify_rational().solve(current_var)[0].rhs()

        return A(tau) + ansatz(prec=precision-1).map_coefficients(lambda term: term.subs(solution_dict).simplify_rational())


    @staticmethod
    @experimental(20050)
    def ImplicitExpansionPeriodicPart(var, phi, period, tau=None, precision=None):
        r"""
        Return the singular expansion for the periodic part of a function `y(z)`
        defined implicitly by `y(z) = z \Phi(y(z))`.

        The function `\Phi` is assumed to be analytic around `0`. Furthermore,
        `\Phi` is not allowed to be an affine-linear function and we require
        `\Phi(0) \neq 0`. For an integer `p`, `\Phi` is called `p`-periodic
        if we have `\Psi(u^p) = \Phi(u)` for a power series `\Psi`
        where `p` is maximal.

        Furthermore, it is assumed that there is a unique positive solution `\tau`
        of `\Phi(\tau) - \tau\Phi'(\tau) = 0`.

        If `\Phi` is `p`-periodic, then we have `y(z) = z g(z^p)`. This method
        returns the singular expansion of `g(z)`.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``phi`` -- the function `\Phi`. See the extended description for
          assumptions on `\Phi`.

        - ``period`` -- the period of the function `\Phi`. See the
          extended description for details.

        - ``tau`` -- (default: ``None``) the fundamental constant described
          in the extended description. If ``None``, then `\tau` is determined
          automatically if possible.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.


        OUTPUT:

        An asymptotic expansion.


        .. NOTE::

            In the given case, the radius of convergence of the function of
            interest is known to be `\rho = \tau/\Phi(\tau)`. Until :trac:`20050`
            is implemented, the variable in the returned asymptotic expansion
            represents a singular element of the form `(1 - z/\rho)^{-1}`,
            for the variable `z\to\rho`.

        .. SEEALSO::

            :meth:`~AsymptoticExpansionGenerators.ImplicitExpansion`,
            :meth:`~AsymptoticExpansionGenerators.InverseFunctionAnalysis`.

        EXAMPLES:

        The generating function enumerating binary trees with respect to
        tree size satisfies `B(z) = z (1 + B(z)^2)`. This function can be
        written as `B(z) = z g(z^2)`, and as `B(z)` can be determined
        explicitly we have `g(z) = \frac{1 - \sqrt{1 - 4z}}{2z}`. We
        compare the corresponding expansions::

            sage: asymptotic_expansions.ImplicitExpansionPeriodicPart('Z', phi=lambda u: 1 + u^2,
            ....:                                                     period=2, precision=7)
            doctest:warning
            ...
            FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See http://trac.sagemath.org/20050 for details.
            2 - 2*Z^(-1/2) + 2*Z^(-1) - 2*Z^(-3/2) + 2*Z^(-2) - 2*Z^(-5/2) + O(Z^(-3))
            sage: def g(z):
            ....:     return (1 - sqrt(1 - 4*z))/(2*z)
            sage: A.<Z> = AsymptoticRing('Z^QQ', QQ, default_prec=3)
            sage: g((1 - 1/Z)/4)
            2 - 2*Z^(-1/2) + 2*Z^(-1) - 2*Z^(-3/2) + 2*Z^(-2) - 2*Z^(-5/2) + O(Z^(-3))

        """
        if tau is None:
            tau = _fundamental_constant_implicit_function_(phi=phi)

        tau_p = tau**period
        aperiodic_expansion = asymptotic_expansions.ImplicitExpansion(var,
                                            phi=lambda u: phi(u**(1/period))**period,
                                            tau=tau_p, precision=precision)

        rho = tau/phi(tau)
        Z = aperiodic_expansion.parent().gen()
        return 1/rho * (aperiodic_expansion/(1 - 1/Z))**(1/period)


    @staticmethod
    def InverseFunctionAnalysis(var, phi, tau=None, period=1, precision=None):
        r"""
        Return the coefficient growth of a function `y(z)` defined implicitly
        by `y(z) = z \Phi(y(z))`.

        The function `\Phi` is assumed to be analytic around `0`. Furthermore,
        `\Phi` is not allowed to be an affine-linear function and we require
        `\Phi(0) \neq 0`. For an integer `p`, `\Phi` is called `p`-periodic
        if we have `\Psi(u^p) = \Phi(u)` for a power series `\Psi`
        where `p` is maximal.

        Furthermore, it is assumed that there is a unique positive solution `\tau`
        of `\Phi(\tau) - \tau\Phi'(\tau) = 0`.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``phi`` -- the function `\Phi`. See the extended description for
          assumptions on `\Phi`.

        - ``tau`` -- (default: ``None``) the fundamental constant described
          in the extended description. If ``None``, then `\tau` is determined
          automatically if possible.

        - ``period`` -- (default: `1`) the period of the function `\Phi`. See
          the extended description for details.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.


        OUTPUT:

        An asymptotic expansion.


        .. NOTE::

            It is not checked that the passed period actually fits to
            the passed function `\Phi`.

            The resulting asymptotic expansion is only valid
            for `n \equiv 1 \mod p`, where `p` is the period. All other
            coefficients are `0`.


        EXAMPLES:

        There are `C_n` (the `n`-th Catalan number) different binary trees
        of size `2n+1`, and there are no binary trees with an even number of
        nodes. The corresponding generating function satisfies
        `B(z) = z (1 + B(z)^2)`, which allows us to compare the asymptotic
        expansions for the number of binary trees of size `n` obtained via
        `C_n` and obtained via the analysis of `B(z)`::

            sage: assume(SR.an_element() > 0)
            sage: A.<n> = AsymptoticRing('QQ^n * n^QQ', SR)
            sage: binomial_expansion = asymptotic_expansions.Binomial_kn_over_n(n, k=2, precision=3)
            sage: catalan_expansion = binomial_expansion / (n+1)
            sage: catalan_expansion.subs(n=(n-1)/2)
            2*sqrt(1/2)/sqrt(pi)*2^n*n^(-3/2) - 3/2*sqrt(1/2)/sqrt(pi)*2^n*n^(-5/2)
            + 25/16*sqrt(1/2)/sqrt(pi)*2^n*n^(-7/2) + O(2^n*n^(-9/2))
            sage: asymptotic_expansions.InverseFunctionAnalysis(n, phi=lambda u: 1 + u^2, period=2,
            ....:                                               tau=1, precision=8)
            2*sqrt(1/2)/sqrt(pi)*2^n*n^(-3/2) - 3/2*sqrt(1/2)/sqrt(pi)*2^n*n^(-5/2)
            + 25/16*sqrt(1/2)/sqrt(pi)*2^n*n^(-7/2) + O(2^n*n^(-9/2))

        The code in the aperiodic case is more efficient, however. Therefore,
        it is recommended to use combinatorial identities to reduce to the
        aperiodic case. In the example above, this is well-known: we now count
        binary trees with `n` internal nodes. The corresponding generating function
        satisfies `B(z) = z (1 + 2B(z) + B(z)^2)`::

            sage: catalan_expansion
            1/sqrt(pi)*4^n*n^(-3/2) - 9/8/sqrt(pi)*4^n*n^(-5/2)
            + 145/128/sqrt(pi)*4^n*n^(-7/2) + O(4^n*n^(-9/2))
            sage: asymptotic_expansions.InverseFunctionAnalysis(n, phi=lambda u: 1 + 2*u + u^2,
            ....:                                               tau=1, precision=8)
            1/sqrt(pi)*4^n*n^(-3/2) - 9/8/sqrt(pi)*4^n*n^(-5/2)
            + 145/128/sqrt(pi)*4^n*n^(-7/2) + O(4^n*n^(-9/2))

            sage: forget()

        .. SEEALSO::

            :meth:`~AsymptoticExpansionGenerators.ImplicitExpansion`,
            :meth:`~AsymptoticExpansionGenerators.ImplicitExpansionPeriodicPart`.


        TESTS:

        Omitting the precision parameter does not lead to an error (per default,
        the default series precision is a python integer, which led to an error
        in an earlier version of the code)::

            sage: set_series_precision(int(5))
            sage: asymptotic_expansions.InverseFunctionAnalysis('n', phi=lambda u: 1 + 2*u + u^2,
            ....:                                               tau=1)
            1/sqrt(pi)*4^n*n^(-3/2) - 9/8/sqrt(pi)*4^n*n^(-5/2) + O(4^n*n^(-3))
        """
        if tau is None:
            tau = _fundamental_constant_implicit_function_(phi=phi)

        rho = tau/phi(tau)

        if period == 1:
            expansion = asymptotic_expansions.ImplicitExpansion(var=var, phi=phi,
                                                                tau=tau, precision=precision)
            return expansion._singularity_analysis_(var, zeta=rho, precision=precision)
        expansion = asymptotic_expansions.ImplicitExpansionPeriodicPart(var=var, phi=phi,
                                                     period=period, tau=tau, precision=precision)
        growth = expansion._singularity_analysis_(var, zeta=rho**period, precision=precision)
        n = growth.parent().gen()
        return growth.subs({n: (n-1)/period})

def _fundamental_constant_implicit_function_(phi):
    r"""
    Return the fundamental constant `\tau` occurring in the analysis of
    implicitly defined functions.

    For a function `y(z)` satisfying `y(z) = z \Phi(y(z))`, the fundamental
    constant `\tau` is the unique positive solution of the equation
    `\Phi(\tau) - \tau \Phi'(\tau) = 0`.

    INPUT:

    - ``phi`` -- the function `\Phi`.

    .. SEEALSO::

        :meth:`~AsymptoticExpansionGenerators.ImplicitExpansion`,
        :meth:`~AsymptoticExpansionGenerators.ImplicitExpansionPeriodicPart`.

    TESTS::

        sage: from sage.rings.asymptotic.asymptotic_expansion_generators \
        ....:     import _fundamental_constant_implicit_function_
        sage: _fundamental_constant_implicit_function_(phi=exp)
        1
        sage: _fundamental_constant_implicit_function_(phi=lambda u: 1 + u^2)
        1
        sage: _fundamental_constant_implicit_function_(phi=lambda u: 1 + 2*u + 2*u^2)
        1/2*sqrt(2)

    """
    from sage.symbolic.ring import SR
    u = SR('u')
    positive_solution = [s for s in (phi(u) - u*phi(u).diff(u)).solve(u)
                         if s.rhs() > 0]
    if len(positive_solution) == 1:
        return positive_solution[0].rhs()
    raise ValueError('Fundamental constant tau could not be determined')

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
