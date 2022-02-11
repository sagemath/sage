r"""
Jacobi Elliptic Functions

This module implements the 12 Jacobi elliptic functions, along with their
inverses and the Jacobi amplitude function.

Jacobi elliptic functions can be thought of as generalizations
of both ordinary and hyperbolic trig functions. There are twelve
Jacobian elliptic functions. Each of the twelve corresponds to an
arrow drawn from one corner of a rectangle to another.

::

                n ------------------- d
                |                     |
                |                     |
                |                     |
                s ------------------- c

Each of the corners of the rectangle are labeled, by convention, ``s``,
``c``, ``d``, and ``n``. The rectangle is understood to be lying on the complex
plane, so that ``s`` is at the origin, ``c`` is on the real axis, and ``n`` is
on the imaginary axis. The twelve Jacobian elliptic functions are
then `\operatorname{pq}(x)`, where ``p`` and ``q`` are one of the letters
``s``, ``c``, ``d``, ``n``.

The Jacobian elliptic functions are then the unique
doubly-periodic, meromorphic functions satisfying the following
three properties:

#. There is a simple zero at the corner ``p``, and a simple pole at the
   corner ``q``.
#. The step from ``p`` to ``q`` is equal to half the period of the function
   `\operatorname{pq}(x)`; that is, the function `\operatorname{pq}(x)` is
   periodic in the direction ``pq``, with the period being twice the distance
   from ``p`` to ``q``. `\operatorname{pq}(x)` is periodic in the other two
   directions as well, with a period such that the distance from ``p`` to one
   of the other corners is a quarter period.
#. If the function `\operatorname{pq}(x)` is expanded in terms of `x` at one of
   the corners, the leading term in the expansion has a coefficient of 1.
   In other words, the leading term of the expansion of `\operatorname{pq}(x)`
   at the corner ``p`` is `x`; the leading term of the expansion at the corner
   ``q`` is `1/x`, and the leading term of an expansion at the other two
   corners is 1.

We can write

.. MATH::

    \operatorname{pq}(x) = \frac{\operatorname{pr}(x)}{\operatorname{qr}(x)}

where ``p``, ``q``, and ``r`` are any of the
letters ``s``, ``c``, ``d``, ``n``, with
the understanding that `\mathrm{ss} = \mathrm{cc} = \mathrm{dd}
= \mathrm{nn} = 1`.

Let

.. MATH::

    u = \int_0^{\phi} \frac{d\theta} {\sqrt {1-m \sin^2 \theta}},

then the *Jacobi elliptic function* `\operatorname{sn}(u)` is given by

.. MATH::

    \operatorname{sn}{u} = \sin{\phi}

and `\operatorname{cn}(u)` is given by

.. MATH::

    \operatorname{cn}{u} = \cos{\phi}

and

.. MATH::

    \operatorname{dn}{u} = \sqrt{1 - m\sin^2 \phi}.

To emphasize the dependence on `m`, one can write
`\operatorname{sn}(u|m)` for example (and similarly for `\mathrm{cn}` and
`\mathrm{dn}`). This is the notation used below.

For a given `k` with `0 < k < 1` they therefore are
solutions to the following nonlinear ordinary differential
equations:

- `\operatorname{sn}\,(x;k)` solves the differential equations

  .. MATH::

      \frac{d^2 y}{dx^2} + (1+k^2) y - 2 k^2 y^3 = 0
      \quad \text{ and } \quad
      \left(\frac{dy}{dx}\right)^2 = (1-y^2) (1-k^2 y^2).

- `\operatorname{cn}(x;k)` solves the differential equations

  .. MATH::

      \frac{d^2 y}{dx^2} + (1-2k^2) y + 2 k^2 y^3 = 0
      \quad \text{ and } \quad
      \left(\frac{dy}{dx}\right)^2 = (1-y^2)(1-k^2 + k^2 y^2).

- `\operatorname{dn}(x;k)` solves the differential equations

  .. MATH::

      \frac{d^2 y}{dx^2} - (2 - k^2) y + 2 y^3 = 0
      \quad \text{ and } \quad
      \left(\frac{dy}{dx}\right)^2 = y^2 (1 - k^2 - y^2).

  If `K(m)` denotes the complete elliptic integral of the
  first kind (named ``elliptic_kc`` in Sage), the elliptic functions
  `\operatorname{sn}(x|m)` and `\operatorname{cn}(x|m)` have real periods
  `4K(m)`, whereas `\operatorname{dn}(x|m)` has a period
  `2K(m)`. The limit `m \rightarrow 0` gives
  `K(0) = \pi/2` and trigonometric functions:
  `\operatorname{sn}(x|0) = \sin{x}`, `\operatorname{cn}(x|0) = \cos{x}`,
  `\operatorname{dn}(x|0) = 1`. The limit `m \rightarrow 1` gives
  `K(1) \rightarrow \infty` and hyperbolic functions:
  `\operatorname{sn}(x|1) = \tanh{x}`,
  `\operatorname{cn}(x|1) = \operatorname{sech}{x}`,
  `\operatorname{dn}(x|1) = \operatorname{sech}{x}`.

REFERENCES:

- :wikipedia:`Jacobi%27s_elliptic_functions`

- [KS2002]_

AUTHORS:

- David Joyner (2006): initial version

- Eviatar Bach (2013): complete rewrite, new numerical evaluation, and
  addition of the Jacobi amplitude function
"""
# ****************************************************************************
#       Copyright (C) 2006 David Joyner <wdj@usna.edu>
#       Copyright (C) 2013 Eviatar Bach <eviatarbach@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.symbolic.function import BuiltinFunction
from sage.functions.trig import (arctan, arcsin, arccos, arccot, arcsec,
                                 arccsc, csc, sec, sin, cos, tan, cot)
from sage.functions.hyperbolic import (arctanh, arccosh, arcsinh, arcsech,
                                       arccsch, arccoth, cosh, coth, sech,
                                       csch, tanh, sinh)
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.functions.special import elliptic_e, elliptic_kc
from sage.libs.mpmath import utils
from sage.misc.latex import latex

HALF = QQ((1, 2))


class Jacobi(BuiltinFunction):
    """
    Base class for the Jacobi elliptic functions.
    """
    def __init__(self, kind):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.functions.jacobi import Jacobi
            sage: Jacobi('sn')
            jacobi_sn

        TESTS::

            sage: N(jacobi("sn", I, 1/2))   # abs tol 1e-12
            -8.59454886300046e-73 + 1.34737147138542*I
        """
        if kind not in ['nd', 'ns', 'nc', 'dn', 'ds', 'dc', 'sn', 'sd',
                        'sc', 'cn', 'cd', 'cs']:
            raise ValueError("kind must be one of 'nd', 'ns', 'nc', 'dn', "
                             "'ds', 'dc', 'sn', 'sd', 'sc', 'cn', 'cd', 'cs'.")
        self.kind = kind
        BuiltinFunction.__init__(self,
                                 name='jacobi_{}'.format(kind),
                                 nargs=2, evalf_params_first=False,
                                 conversions=dict(maple=
                                                  ('Jacobi{}'
                                                   .format(kind.upper())),
                                                  mathematica=
                                                  ('Jacobi{}'
                                                   .format(kind.upper())),
                                                  maxima=
                                                  ('jacobi_{}'
                                                   .format(kind))))

    def _eval_(self, x, m):
        r"""
        TESTS:

        Check that the simplifications are correct::

            sage: from mpmath import almosteq
            sage: almosteq(n(jacobi_nd(8, 0, hold=True)), n(jacobi_nd(8, 0)))
            True
            sage: almosteq(n(jacobi_nd(1, 1, hold=True)), n(jacobi_nd(1, 1)))
            True
            sage: almosteq(n(jacobi_nd(0, -5, hold=True)), n(jacobi_nd(0, -5)))
            True
            sage: almosteq(n(jacobi_ns(-4, 0, hold=True)), n(jacobi_ns(-4, 0)))
            True
            sage: almosteq(n(jacobi_ns(-2, 1, hold=True)), n(jacobi_ns(-2, 1)))
            True
            sage: almosteq(n(jacobi_nc(2, 0, hold=True)), n(jacobi_nc(2, 0)))
            True
            sage: almosteq(n(jacobi_nc(1, 1, hold=True)), n(jacobi_nc(1, 1)))
            True
            sage: almosteq(n(jacobi_nc(0, 0, hold=True)), n(jacobi_nc(0, 0)))
            True
            sage: almosteq(n(jacobi_dn(-10, 0, hold=True)), n(jacobi_dn(-10, 0)))
            True
            sage: almosteq(n(jacobi_dn(-1, 1, hold=True)), n(jacobi_dn(-1, 1)))
            True
            sage: almosteq(n(jacobi_dn(0, 3, hold=True)), n(jacobi_dn(0, 3)))
            True
            sage: almosteq(n(jacobi_ds(2, 0, hold=True)), n(jacobi_ds(2, 0)))
            True
            sage: almosteq(n(jacobi_dc(-1, 0, hold=True)), n(jacobi_dc(-1, 0)))
            True
            sage: almosteq(n(jacobi_dc(-8, 1, hold=True)), n(jacobi_dc(-8, 1)))
            True
            sage: almosteq(n(jacobi_dc(0, -10, hold=True)), n(jacobi_dc(0, -10)))
            True
            sage: almosteq(n(jacobi_sn(-7, 0, hold=True)), n(jacobi_sn(-7, 0)))
            True
            sage: almosteq(n(jacobi_sn(-3, 1, hold=True)), n(jacobi_sn(-3, 1)))
            True
            sage: almosteq(n(jacobi_sn(0, -6, hold=True)), n(jacobi_sn(0, -6)))
            True
            sage: almosteq(n(jacobi_sd(4, 0, hold=True)), n(jacobi_sd(4, 0)))
            True
            sage: almosteq(n(jacobi_sd(0, 1, hold=True)), n(jacobi_sd(0, 1)))
            True
            sage: almosteq(n(jacobi_sd(0, 3, hold=True)), n(jacobi_sd(0, 3)))
            True
            sage: almosteq(n(jacobi_sc(-9, 0, hold=True)), n(jacobi_sc(-9, 0)))
            True
            sage: almosteq(n(jacobi_sc(0, 1, hold=True)), n(jacobi_sc(0, 1)))
            True
            sage: almosteq(n(jacobi_sc(0, -10, hold=True)), n(jacobi_sc(0, -10)))
            True
            sage: almosteq(n(jacobi_cn(-2, 0, hold=True)), n(jacobi_cn(-2, 0)))
            True
            sage: almosteq(n(jacobi_cn(6, 1, hold=True)), n(jacobi_cn(6, 1)))
            True
            sage: almosteq(n(jacobi_cn(0, -10, hold=True)), n(jacobi_cn(0, -10)))
            True
            sage: almosteq(n(jacobi_cd(9, 0, hold=True)), n(jacobi_cd(9, 0)))
            True
            sage: almosteq(n(jacobi_cd(-8, 1, hold=True)), n(jacobi_cd(-8, 1)))
            True
            sage: almosteq(n(jacobi_cd(0, 1, hold=True)), n(jacobi_cd(0, 1)))
            True
            sage: almosteq(n(jacobi_cs(-9, 0, hold=True)), n(jacobi_cs(-9, 0)))
            True
            sage: almosteq(n(jacobi_cs(-6, 1, hold=True)), n(jacobi_cs(-6, 1)))
            True
        """
        if self.kind == 'nd':
            if m == 0:
                return Integer(1)
            elif m == 1:
                return cosh(x)
            elif x == 0:
                return Integer(1)
        elif self.kind == 'ns':
            if m == 0:
                return csc(x)
            elif m == 1:
                return coth(x)
        elif self.kind == 'nc':
            if m == 0:
                return sec(x)
            elif m == 1:
                return cosh(x)
            elif x == 0:
                return Integer(1)
        elif self.kind == 'dn':
            if m == 0:
                return Integer(1)
            elif m == 1:
                return sech(x)
            elif x == 0:
                return Integer(1)
        elif self.kind == 'ds':
            if m == 0:
                return csc(x)
        elif self.kind == 'dc':
            if m == 0:
                return sec(x)
            elif m == 1:
                return Integer(1)
            elif x == 0:
                return Integer(1)
        elif self.kind == 'sn':
            if m == 0:
                return sin(x)
            elif m == 1:
                return tanh(x)
            elif x == 0:
                return Integer(0)
        elif self.kind == 'sd':
            if m == 0:
                return sin(x)
            elif m == 1:
                return sinh(x)
            elif x == 0:
                return Integer(0)
        elif self.kind == 'sc':
            if m == 0:
                return tan(x)
            elif m == 1:
                return sinh(x)
            elif x == 0:
                return Integer(0)
        elif self.kind == 'cn':
            if m == 0:
                return cos(x)
            elif m == 1:
                return sech(x)
            elif x == 0:
                return Integer(1)
        elif self.kind == 'cd':
            if m == 0:
                return cos(x)
            elif m == 1:
                return Integer(1)
            elif x == 0:
                return Integer(1)
        elif self.kind == 'cs':
            if m == 0:
                return cot(x)
            elif m == 1:
                return csch(x)
        return

    def _evalf_(self, x, m, parent, algorithm=None):
        r"""
        TESTS::

            sage: jacobi_sn(3, 4).n(100)
            -0.33260000892770027112809652714 + 1.7077912301715219199143891076e-33*I
            sage: jacobi_dn(I, I).n()
            0.874189950651018 + 0.667346865048825*I
        """
        from mpmath import ellipfun
        return utils.call(ellipfun, self.kind, x, m, parent=parent)

    def _derivative_(self, x, m, diff_param):
        r"""
        TESTS:

        sn, cn, and dn are analytic for all real ``x``, so we can check
        that the derivatives are correct by computing the series::

            sage: from mpmath import almosteq
            sage: a = 0.9327542442482303
            sage: b = 0.7402326293643771
            sage: almosteq(jacobi_sn(x, b).series(x, 10).subs(x=a),
            ....:          jacobi_sn(a, b), abs_eps=0.01)
            True
            sage: almosteq(jacobi_cn(x, b).series(x, 10).subs(x=a),
            ....:          jacobi_cn(a, b), abs_eps=0.01)
            True
            sage: almosteq(jacobi_dn(x, b).series(x, 10).subs(x=a),
            ....:          jacobi_dn(a, b), abs_eps=0.01)
            True
        """
        if diff_param == 0:
            # From Wolfram Functions Site
            if self.kind == 'cd':
                return (m - Integer(1)) * jacobi_nd(x, m) * jacobi_sd(x, m)
            elif self.kind == 'cn':
                return -jacobi_sn(x, m) * jacobi_dn(x, m)
            elif self.kind == 'cs':
                return -jacobi_ds(x, m) * jacobi_ns(x, m)
            elif self.kind == 'dc':
                return (Integer(1) - m) * jacobi_nc(x, m) * jacobi_sc(x, m)
            elif self.kind == 'dn':
                return -m * jacobi_sn(x, m) * jacobi_cn(x, m)
            elif self.kind == 'ds':
                return -jacobi_cs(x, m) * jacobi_ns(x, m)
            elif self.kind == 'nc':
                return jacobi_dc(x, m) * jacobi_sc(x, m)
            elif self.kind == 'nd':
                return m * jacobi_cd(x, m) * jacobi_sd(x, m)
            elif self.kind == 'ns':
                return -jacobi_cs(x, m) * jacobi_ds(x, m)
            elif self.kind == 'sc':
                return jacobi_dc(x, m) * jacobi_nc(x, m)
            elif self.kind == 'sd':
                return jacobi_cd(x, m) * jacobi_nd(x, m)
            elif self.kind == 'sn':
                return jacobi_cn(x, m) * jacobi_dn(x, m)
        elif diff_param == 1:
            # From Maxima
            if self.kind == 'nd':
                return (HALF*((x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                              (m - Integer(1)))*jacobi_sn(x, m)*jacobi_cn(x, m) -
                              jacobi_dn(x, m)*jacobi_sn(x, m)**Integer(2)/(m - Integer(1)))/
                        jacobi_dn(x, m)**Integer(2))
            elif self.kind == 'ns':
                return (HALF*(jacobi_sn(x, m)*jacobi_cn(x, m)**Integer(2)/(m - Integer(1)) -
                              (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                               (m - Integer(1)))*jacobi_dn(x, m)*jacobi_cn(x, m)/m)/
                        jacobi_sn(x, m)**Integer(2))
            elif self.kind == 'nc':
                return (-HALF*(jacobi_sn(x, m)**Integer(2)*jacobi_cn(x, m)/(m - Integer(1)) -
                               (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                                (m - Integer(1)))*jacobi_dn(x, m)*
                               jacobi_sn(x, m)/m)/jacobi_cn(x, m)**Integer(2))
            elif self.kind == 'dn':
                return (-HALF*(x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                               (m - Integer(1)))*jacobi_sn(x, m)*jacobi_cn(x, m) +
                        HALF*jacobi_dn(x, m)*jacobi_sn(x, m)**Integer(2)/(m - Integer(1)))
            elif self.kind == 'ds':
                return (HALF*(jacobi_sn(x, m)*jacobi_cn(x, m)**Integer(2)/(m - Integer(1)) -
                        (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                         (m - Integer(1)))*jacobi_dn(x, m)*jacobi_cn(x, m)/m)*
                        jacobi_dn(x, m)/jacobi_sn(x, m)**Integer(2) -
                        HALF*((x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                              (m - Integer(1)))*jacobi_sn(x, m)*jacobi_cn(x, m) -
                              jacobi_dn(x, m)*jacobi_sn(x, m)**Integer(2)/(m - Integer(1)))/
                        jacobi_sn(x, m))
            elif self.kind == 'dc':
                return (-HALF*(jacobi_sn(x, m)**Integer(2)*jacobi_cn(x, m)/(m - Integer(1)) -
                               (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                                (m - Integer(1)))*jacobi_dn(x, m)*
                               jacobi_sn(x, m)/m)*jacobi_dn(x, m)/
                        jacobi_cn(x, m)**Integer(2) -
                        HALF*((x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                               (m - Integer(1)))*jacobi_sn(x, m)*jacobi_cn(x, m) -
                              jacobi_dn(x, m)*jacobi_sn(x, m)**Integer(2)/(m - Integer(1)))/
                        jacobi_cn(x, m))
            elif self.kind == 'sn':
                return (-HALF*jacobi_sn(x, m)*jacobi_cn(x, m)**Integer(2)/(m - Integer(1)) +
                        HALF*(x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                              (m - Integer(1)))*jacobi_dn(x, m)*jacobi_cn(x, m)/m)
            elif self.kind == 'sd':
                return (-HALF*(jacobi_sn(x, m)*jacobi_cn(x, m)**Integer(2)/(m - Integer(1)) -
                        (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                         (m - Integer(1)))*jacobi_dn(x, m)*jacobi_cn(x, m)/m)/
                        jacobi_dn(x, m) + HALF*
                        ((x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                          (m - Integer(1)))*jacobi_sn(x, m)*jacobi_cn(x, m) -
                         jacobi_dn(x, m)*jacobi_sn(x, m)**Integer(2)/(m - Integer(1)))*
                        jacobi_sn(x, m)/jacobi_dn(x, m)**Integer(2))
            elif self.kind == 'sc':
                return (-HALF*(jacobi_sn(x, m)*jacobi_cn(x, m)**Integer(2)/(m - Integer(1)) -
                               (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                                (m - Integer(1)))*jacobi_dn(x, m)*
                               jacobi_cn(x, m)/m)/jacobi_cn(x, m) -
                        HALF*(jacobi_sn(x, m)**Integer(2)*jacobi_cn(x, m)/(m - Integer(1)) -
                              (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                               (m - Integer(1)))*jacobi_dn(x, m)*jacobi_sn(x, m)/m)*
                        jacobi_sn(x, m)/jacobi_cn(x, m)**Integer(2))
            elif self.kind == 'cn':
                return (HALF*jacobi_sn(x, m)**Integer(2)*jacobi_cn(x, m)/(m - Integer(1)) -
                        HALF*(x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                              (m - Integer(1)))*jacobi_dn(x, m)*jacobi_sn(x, m)/m)
            elif self.kind == 'cd':
                return (HALF*(jacobi_sn(x, m)**Integer(2)*jacobi_cn(x, m)/(m - Integer(1)) -
                        (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                         (m - Integer(1)))*jacobi_dn(x, m)*jacobi_sn(x, m)/m)/
                        jacobi_dn(x, m) +
                        HALF*((x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                               (m - Integer(1)))*jacobi_sn(x, m)*jacobi_cn(x, m) -
                              jacobi_dn(x, m)*jacobi_sn(x, m)**Integer(2)/(m - Integer(1)))*
                        jacobi_cn(x, m)/jacobi_dn(x, m)**Integer(2))
            elif self.kind == 'cs':
                return (HALF*(jacobi_sn(x, m)*jacobi_cn(x, m)**Integer(2)/(m - Integer(1)) -
                        (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                         (m - Integer(1)))*jacobi_dn(x, m)*jacobi_cn(x, m)/m)*
                        jacobi_cn(x, m)/jacobi_sn(x, m)**Integer(2) +
                        HALF*(jacobi_sn(x, m)**Integer(2)*jacobi_cn(x, m)/(m - Integer(1)) -
                              (x + elliptic_e(arcsin(jacobi_sn(x, m)), m)/
                               (m - Integer(1)))*jacobi_dn(x, m)*jacobi_sn(x, m)/m)/
                        jacobi_sn(x, m))

    def _latex_(self):
        r"""
        TESTS::

            sage: latex(jacobi_sn)
            \operatorname{sn}
        """
        return r"\operatorname{{{}}}".format(self.kind)

    def _print_latex_(self, x, m):
        r"""
        TESTS::

            sage: latex(jacobi_sn(x, 3))
            \operatorname{sn}\left(x\middle|3\right)
        """
        return r"\operatorname{{{}}}\left({}\middle|{}\right)".format(self.kind,
                                                                      latex(x),
                                                                      latex(m))

jacobi_nd = Jacobi('nd')
jacobi_ns = Jacobi('ns')
jacobi_nc = Jacobi('nc')
jacobi_dn = Jacobi('dn')
jacobi_ds = Jacobi('ds')
jacobi_dc = Jacobi('dc')
jacobi_sn = Jacobi('sn')
jacobi_sd = Jacobi('sd')
jacobi_sc = Jacobi('sc')
jacobi_cn = Jacobi('cn')
jacobi_cd = Jacobi('cd')
jacobi_cs = Jacobi('cs')


class InverseJacobi(BuiltinFunction):
    r"""
    Base class for the inverse Jacobi elliptic functions.
    """
    def __init__(self, kind):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.functions.jacobi import InverseJacobi
            sage: InverseJacobi('sn')
            inverse_jacobi_sn
        """
        if kind not in ['nd', 'ns', 'nc', 'dn', 'ds', 'dc', 'sn', 'sd',
                        'sc', 'cn', 'cd', 'cs']:
            raise ValueError("kind must be one of 'nd', 'ns', 'nc', 'dn', "
                             "'ds', 'dc', 'sn', 'sd', 'sc', 'cn', 'cd', 'cs'.")
        self.kind = kind
        BuiltinFunction.__init__(self,
                                 name='inverse_jacobi_{}'.format(kind),
                                 nargs=2, evalf_params_first=False,
                                 conversions=dict(maple=
                                                  ('InverseJacobi{}'
                                                   .format(kind.upper())),
                                                  mathematica=
                                                  ('InverseJacobi{}'
                                                   .format(kind.upper())),
                                                  maxima=
                                                  ('inverse_jacobi_{}'
                                                   .format(kind))))

    def _eval_(self, x, m):
        r"""
        TESTS:

        Check that the simplifications are correct::

            sage: from mpmath import almosteq
            sage: almosteq(n(inverse_jacobi_cd(1, -8, hold=True)),
            ....:          n(inverse_jacobi_cd(1, -8)))
            True
            sage: almosteq(n(inverse_jacobi_cn(0, -5, hold=True)),
            ....:          n(inverse_jacobi_cn(0, -5)))
            True
            sage: almosteq(n(inverse_jacobi_cn(1, -8, hold=True)),
            ....:          n(inverse_jacobi_cn(1, -8)))
            True
            sage: almosteq(n(inverse_jacobi_cs(7, 1, hold=True)),
            ....:          n(inverse_jacobi_cs(7, 1)))
            True
            sage: almosteq(n(inverse_jacobi_dc(3, 0, hold=True)),
            ....:          n(inverse_jacobi_dc(3, 0)))
            True
            sage: almosteq(n(inverse_jacobi_dc(1, 7, hold=True)),
            ....:          n(inverse_jacobi_dc(1, 7)))
            True
            sage: almosteq(n(inverse_jacobi_dn(1, -1, hold=True)),
            ....:          n(inverse_jacobi_dn(1, -1)))
            True
            sage: almosteq(n(inverse_jacobi_ds(7, 0, hold=True)),
            ....:          n(inverse_jacobi_ds(7, 0)))
            True
            sage: almosteq(n(inverse_jacobi_ds(5, 1, hold=True)),
            ....:          n(inverse_jacobi_ds(5, 1)))
            True
            sage: almosteq(n(inverse_jacobi_nc(-2, 0, hold=True)),
            ....:          n(inverse_jacobi_nc(-2, 0)))
            True
            sage: almosteq(n(inverse_jacobi_nc(-1, 1, hold=True)),
            ....:          n(inverse_jacobi_nc(-1, 1)))
            True
            sage: almosteq(n(inverse_jacobi_nc(1, 4, hold=True)),
            ....:          n(inverse_jacobi_nc(1, 4)))
            True
            sage: almosteq(n(inverse_jacobi_nd(9, 1, hold=True)),
            ....:          n(inverse_jacobi_nd(9, 1)))
            True
            sage: almosteq(n(inverse_jacobi_nd(1, -9, hold=True)),
            ....:          n(inverse_jacobi_nd(1, -9)))
            True
            sage: almosteq(n(inverse_jacobi_ns(-6, 0, hold=True)),
            ....:          n(inverse_jacobi_ns(-6, 0)))
            True
            sage: almosteq(n(inverse_jacobi_ns(6, 1, hold=True)),
            ....:          n(inverse_jacobi_ns(6, 1)))
            True
            sage: almosteq(n(inverse_jacobi_sc(9, 0, hold=True)),
            ....:          n(inverse_jacobi_sc(9, 0)))
            True
            sage: almosteq(n(inverse_jacobi_sc(8, 1, hold=True)),
            ....:          n(inverse_jacobi_sc(8, 1)))
            True
            sage: almosteq(n(inverse_jacobi_sc(0, -8, hold=True)),
            ....:          n(inverse_jacobi_sc(0, -8)))
            True
            sage: almosteq(n(inverse_jacobi_sd(-1, 0, hold=True)),
            ....:          n(inverse_jacobi_sd(-1, 0)))
            True
            sage: almosteq(n(inverse_jacobi_sd(-2, 1, hold=True)),
            ....:          n(inverse_jacobi_sd(-2, 1)))
            True
            sage: almosteq(n(inverse_jacobi_sd(0, -2, hold=True)),
            ....:          n(inverse_jacobi_sd(0, -2)))
            True
            sage: almosteq(n(inverse_jacobi_sn(0, 0, hold=True)),
            ....:          n(inverse_jacobi_sn(0, 0)))
            True
            sage: almosteq(n(inverse_jacobi_sn(0, 6, hold=True)),
            ....:          n(inverse_jacobi_sn(0, 6)))
            True
        """
        if self.kind == 'cd':
            if m == 0:
                return arccos(x)
            elif x == 1:
                return Integer(0)
        elif self.kind == 'cn':
            if m == 0:
                return arccos(x)
            elif m == 1:
                return arcsech(x)
            elif x == 0:
                return elliptic_kc(m)
            elif x == 1:
                return Integer(0)
        elif self.kind == 'cs':
            if m == 0:
                return arccot(x)
            elif m == 1:
                return arccsch(x)
        elif self.kind == 'dc':
            if m == 0:
                return arcsec(x)
            elif x == 1:
                return Integer(0)
        elif self.kind == 'dn':
            if m == 1:
                return arcsech(x)
            elif x == 1:
                return Integer(0)
        elif self.kind == 'ds':
            if m == 0:
                return arccsc(x)
            elif m == 1:
                return arccsch(x)
        elif self.kind == 'nc':
            if m == 0:
                return arcsec(x)
            elif m == 1:
                return arccosh(x)
            elif x == 1:
                return Integer(0)
        elif self.kind == 'nd':
            if m == 1:
                return arccosh(x)
            elif x == 1:
                return Integer(0)
        elif self.kind == 'ns':
            if m == 0:
                return arccsc(x)
            elif m == 1:
                return arccoth(x)
        elif self.kind == 'sc':
            if m == 0:
                return arctan(x)
            elif m == 1:
                return arcsinh(x)
            elif x == 0:
                return Integer(0)
        elif self.kind == 'sd':
            if m == 0:
                return arcsin(x)
            elif m == 1:
                return arcsinh(x)
            elif x == 0:
                return Integer(0)
        elif self.kind == 'sn':
            if m == 0:
                return arcsin(x)
            elif m == 1:
                return arctanh(x)
            elif x == 0:
                return Integer(0)
        return

    def _evalf_(self, x, m, parent, algorithm=None):
        r"""
        TESTS::

            sage: inverse_jacobi_cn(2, 3).n()
            0.859663746362987*I
            sage: inverse_jacobi_cd(3, 4).n(100)
            -0.67214752201235862490069823239 + 2.1565156474996432354386749988*I
        """
        return utils.call(inverse_jacobi_f, self.kind, x, m, parent=parent)

    def _derivative_(self, x, m, diff_param):
        r"""
        TESTS:

        Check that ``dy/dx * dx/dy == 1``, where ``y = jacobi_pq(x, m)`` and
        ``x = inverse_jacobi_pq(y, m)``::

            sage: from mpmath import almosteq
            sage: a = 0.130103220857094
            sage: b = 0.437176765041986
            sage: m = var('m')
            sage: almosteq(abs((diff(jacobi_cd(x, m), x) *
            ....:               diff(inverse_jacobi_cd(x, m), x).subs(x=jacobi_cd(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_cn(x, m), x) *
            ....:               diff(inverse_jacobi_cn(x, m), x).subs(x=jacobi_cn(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_cs(x, m), x) *
            ....:               diff(inverse_jacobi_cs(x, m), x).subs(x=jacobi_cs(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_dc(x, m), x) *
            ....:               diff(inverse_jacobi_dc(x, m), x).subs(x=jacobi_dc(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_dn(x, m), x) *
            ....:               diff(inverse_jacobi_dn(x, m), x).subs(x=jacobi_dn(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_ds(x, m), x) *
            ....:               diff(inverse_jacobi_ds(x, m), x).subs(x=jacobi_ds(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_nc(x, m), x) *
            ....:               diff(inverse_jacobi_nc(x, m), x).subs(x=jacobi_nc(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_nd(x, m), x) *
            ....:               diff(inverse_jacobi_nd(x, m), x).subs(x=jacobi_nd(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_ns(x, m), x) *
            ....:               diff(inverse_jacobi_ns(x, m), x).subs(x=jacobi_ns(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_sc(x, m), x) *
            ....:               diff(inverse_jacobi_sc(x, m), x).subs(x=jacobi_sc(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_sd(x, m), x) *
            ....:               diff(inverse_jacobi_sd(x, m), x).subs(x=jacobi_sd(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
            sage: almosteq(abs((diff(jacobi_sn(x, m), x) *
            ....:               diff(inverse_jacobi_sn(x, m), x).subs(x=jacobi_sn(x, m))).subs(x=a, m=b)),
            ....:          1, abs_eps=1e-14)
            True
        """
        # From Wolfram Functions Site
        if diff_param == 0:
            if self.kind == 'cd':
                return (jacobi_sn(inverse_jacobi_cd(x, m), m) /
                        (x ** Integer(2) - Integer(1)))
            elif self.kind == 'cn':
                return (jacobi_ds(inverse_jacobi_cn(x, m), m) /
                        (m * x ** Integer(2) - m + Integer(1)))
            elif self.kind == 'cs':
                return (jacobi_nd(inverse_jacobi_cs(x, m), m) /
                        (x ** Integer(2) + Integer(1)))
            elif self.kind == 'dc':
                return (jacobi_sn(inverse_jacobi_dc(x, m), m) /
                        (x ** Integer(2) - Integer(1)))
            elif self.kind == 'dn':
                return -(jacobi_cs(inverse_jacobi_dn(x, m), m) /
                         (x ** Integer(2) + m - Integer(1)))
            elif self.kind == 'ds':
                return (jacobi_nc(inverse_jacobi_ds(x, m), m) /
                        (x ** Integer(2) + m))
            elif self.kind == 'nc':
                return (jacobi_ds(inverse_jacobi_nc(x, m), m) /
                        (-m * x ** Integer(2) + x ** Integer(2) + m))
            elif self.kind == 'nd':
                return (jacobi_sc(inverse_jacobi_nd(x, m), m) /
                        (x ** Integer(2) - Integer(1)))
            elif self.kind == 'ns':
                return Integer(1) / (jacobi_cs(inverse_jacobi_ns(x, m), m) *
                                     jacobi_ds(inverse_jacobi_ns(x, m), m))
            elif self.kind == 'sc':
                return (jacobi_nd(inverse_jacobi_sc(x, m), m) /
                        (x ** Integer(2) + Integer(1)))
            elif self.kind == 'sd':
                return (jacobi_cn(inverse_jacobi_sd(x, m), m) /
                        ((m - Integer(1)) * x ** Integer(2) + Integer(1)))
            elif self.kind == 'sn':
                return (jacobi_cd(inverse_jacobi_sn(x, m), m) /
                        (Integer(1) - x ** Integer(2)))
        elif diff_param == 1:
            if self.kind == 'cd':
                return ((Integer(1) / (Integer(2) * (Integer(1) - m) * m)) *
                        ((m - Integer(1)) * inverse_jacobi_cd(x, m) +
                         elliptic_e(jacobi_am(inverse_jacobi_cd(x, m), m),
                                    m)))
            elif self.kind == 'cn':
                return ((-(Integer(1) / (Integer(2) * (-Integer(1) + m) * m))) *
                        (elliptic_e(jacobi_am(inverse_jacobi_cn(x, m), m),
                                    m) + (-Integer(1) + m) *
                         inverse_jacobi_cn(x, m) - m * x *
                         jacobi_sd(inverse_jacobi_cn(x, m), m)))
            elif self.kind == 'cs':
                return ((-(Integer(1) / (Integer(2) * (-Integer(1) + m) * m * (Integer(1) + x ** Integer(2))))) *
                        ((Integer(1) + x ** Integer(2)) *
                         elliptic_e(jacobi_am(inverse_jacobi_cs(x, m), m),
                                    m) + (-Integer(1) + m) * (Integer(1) + x ** Integer(2)) *
                         inverse_jacobi_cs(x, m) - m * x *
                         jacobi_nd(inverse_jacobi_cs(x, m), m)))
            elif self.kind == 'dc':
                return ((Integer(1) / (Integer(2) * (Integer(1) - m) * m)) *
                        (elliptic_e(jacobi_am(inverse_jacobi_dc(x, m), m),
                                    m) - (Integer(1) - m) *
                         inverse_jacobi_dc(x, m)))
            elif self.kind == 'dn':
                return ((Integer(1) / (Integer(2) * (Integer(1) - m) * m)) * ((m - Integer(1)) *
                        inverse_jacobi_dn(x, m) +
                        elliptic_e(jacobi_am(inverse_jacobi_dn(x, m), m), m) -
                        x * jacobi_sc(inverse_jacobi_dn(x, m), m)))
            elif self.kind == 'ds':
                return ((-(Integer(1) / (Integer(2) * (-Integer(1) + m) * m))) *
                        (elliptic_e(jacobi_am(inverse_jacobi_ds(x, m), m), m) +
                         (-Integer(1) + m) * inverse_jacobi_ds(x, m) -
                         (m * x * jacobi_nc(inverse_jacobi_ds(x, m), m)) /
                         (m + x ** Integer(2))))
            elif self.kind == 'nc':
                return ((Integer(1) / (Integer(2) * (-Integer(1) + m) * m * x)) * ((-x) *
                        (elliptic_e(jacobi_am(inverse_jacobi_nc(x, m), m), m) +
                         (-Integer(1) + m) * inverse_jacobi_nc(x, m)) + m *
                        jacobi_sd(inverse_jacobi_nc(x, m), m)))
            elif self.kind == 'nd':
                return ((Integer(1) / (Integer(2) * (m - Integer(1)) * m)) *
                        ((Integer(1) - m) * inverse_jacobi_nd(x, m) -
                         elliptic_e(jacobi_am(inverse_jacobi_nd(x, m), m), m) +
                         (Integer(1) / x) * jacobi_sc(inverse_jacobi_nd(x, m), m)))
            elif self.kind == 'ns':
                return ((Integer(1)/(Integer(2) * (m - Integer(1)) * m)) *
                        ((Integer(1) - m) * inverse_jacobi_ns(x, m) -
                         elliptic_e(jacobi_am(inverse_jacobi_ns(x, m), m), m) +
                         (m / x) * jacobi_cd(inverse_jacobi_ns(x, m), m)))
            elif self.kind == 'sc':
                return ((-(Integer(1) / (Integer(2) * (-Integer(1) + m) * m * (Integer(1) + x ** Integer(2))))) *
                        ((Integer(1) + x ** Integer(2)) *
                         elliptic_e(jacobi_am(inverse_jacobi_sc(x, m), m), m) +
                         (-Integer(1) + m) * (Integer(1) + x ** Integer(2)) * inverse_jacobi_sc(x, m) -
                         m * x * jacobi_nd(inverse_jacobi_sc(x, m), m)))
            elif self.kind == 'sd':
                return ((-(Integer(1) / (Integer(2) * (-Integer(1) + m) * m))) *
                        (elliptic_e(jacobi_am(inverse_jacobi_sd(x, m), m), m) +
                         (-Integer(1) + m) * inverse_jacobi_sd(x, m) -
                         (m * x * jacobi_nc(inverse_jacobi_sd(x, m), m)) /
                         (Integer(1) + m * x ** Integer(2))))
            elif self.kind == 'sn':
                return ((Integer(1) / (Integer(2) * (Integer(1) - m) * m)) *
                        (elliptic_e(jacobi_am(inverse_jacobi_sn(x, m), m), m) +
                         (-Integer(1) + m) * inverse_jacobi_sn(x, m) - m * x *
                         jacobi_cd(inverse_jacobi_sn(x, m), m)))

    def _latex_(self):
        r"""
        TESTS::

            sage: latex(inverse_jacobi_dn)
            \operatorname{arcdn}
        """
        return r"\operatorname{{arc{}}}".format(self.kind)

    def _print_latex_(self, x, m):
        r"""
        TESTS::

            sage: latex(inverse_jacobi_dn(x, 3))
            \operatorname{arcdn}\left(x\middle|3\right)
        """
        return r"\operatorname{{arc{}}}\left({}\middle|{}\right)".format(self.kind,
                                                                         latex(x),
                                                                         latex(m))

inverse_jacobi_nd = InverseJacobi('nd')
inverse_jacobi_ns = InverseJacobi('ns')
inverse_jacobi_nc = InverseJacobi('nc')
inverse_jacobi_dn = InverseJacobi('dn')
inverse_jacobi_ds = InverseJacobi('ds')
inverse_jacobi_dc = InverseJacobi('dc')
inverse_jacobi_sn = InverseJacobi('sn')
inverse_jacobi_sd = InverseJacobi('sd')
inverse_jacobi_sc = InverseJacobi('sc')
inverse_jacobi_cn = InverseJacobi('cn')
inverse_jacobi_cd = InverseJacobi('cd')
inverse_jacobi_cs = InverseJacobi('cs')


def jacobi(kind, z, m, **kwargs):
    r"""
    The 12 Jacobi elliptic functions.

    INPUT:

    - ``kind`` -- a string of the form ``'pq'``, where ``p``, ``q`` are in
      ``c``, ``d``, ``n``, ``s``
    - ``z`` -- a complex number
    - ``m`` -- a complex number; note that `m = k^2`, where `k` is
      the elliptic modulus

    EXAMPLES::

        sage: jacobi('sn', 1, 1)
        tanh(1)
        sage: jacobi('cd', 1, 1/2)
        jacobi_cd(1, 1/2)
        sage: RDF(jacobi('cd', 1, 1/2))
        0.7240097216593705
        sage: (RDF(jacobi('cn', 1, 1/2)), RDF(jacobi('dn', 1, 1/2)),
        ....:  RDF(jacobi('cn', 1, 1/2) / jacobi('dn', 1, 1/2)))
        (0.5959765676721407, 0.8231610016315962, 0.7240097216593705)
        sage: jsn = jacobi('sn', x, 1)
        sage: P = plot(jsn, 0, 1)
    """
    if kind == 'nd':
        return jacobi_nd(z, m, **kwargs)
    elif kind == 'ns':
        return jacobi_ns(z, m, **kwargs)
    elif kind == 'nc':
        return jacobi_nc(z, m, **kwargs)
    elif kind == 'dn':
        return jacobi_dn(z, m, **kwargs)
    elif kind == 'ds':
        return jacobi_ds(z, m, **kwargs)
    elif kind == 'dc':
        return jacobi_dc(z, m, **kwargs)
    elif kind == 'sn':
        return jacobi_sn(z, m, **kwargs)
    elif kind == 'sd':
        return jacobi_sd(z, m, **kwargs)
    elif kind == 'sc':
        return jacobi_sc(z, m, **kwargs)
    elif kind == 'cn':
        return jacobi_cn(z, m, **kwargs)
    elif kind == 'cd':
        return jacobi_cd(z, m, **kwargs)
    elif kind == 'cs':
        return jacobi_cs(z, m, **kwargs)
    else:
        raise ValueError("kind must be one of 'nd', 'ns', 'nc', 'dn', "
                         "'ds', 'dc', 'sn', 'sd', 'sc', 'cn', 'cd', 'cs'.")


def inverse_jacobi(kind, x, m, **kwargs):
    r"""
    The inverses of the 12 Jacobi elliptic functions. They have the property
    that

    .. MATH::

        \operatorname{pq}(\operatorname{arcpq}(x|m)|m) =
        \operatorname{pq}(\operatorname{pq}^{-1}(x|m)|m) = x.

    INPUT:

    - ``kind`` -- a string of the form ``'pq'``, where ``p``, ``q`` are in
      ``c``, ``d``, ``n``, ``s``
    - ``x`` -- a real number
    - ``m`` -- a real number; note that `m = k^2`, where `k` is the elliptic
      modulus

    EXAMPLES::

        sage: jacobi('dn', inverse_jacobi('dn', 3, 0.4), 0.4)
        3.00000000000000
        sage: inverse_jacobi('dn', 10, 1/10).n(digits=50)
        2.4777736267904273296523691232988240759001423661683*I
        sage: inverse_jacobi_dn(x, 1)
        arcsech(x)
        sage: inverse_jacobi_dn(1, 3)
        0
        sage: m = var('m')
        sage: z = inverse_jacobi_dn(x, m).series(x, 4).subs(x=0.1, m=0.7)
        sage: jacobi_dn(z, 0.7)
        0.0999892750039819...
        sage: inverse_jacobi_nd(x, 1)
        arccosh(x)
        sage: inverse_jacobi_nd(1, 2)
        0
        sage: inverse_jacobi_ns(10^-5, 3).n()
        5.77350269202456e-6 + 1.17142008414677*I
        sage: jacobi('sn', 1/2, 1/2)
        jacobi_sn(1/2, 1/2)
        sage: jacobi('sn', 1/2, 1/2).n()
        0.470750473655657
        sage: inverse_jacobi('sn', 0.47, 1/2)
        0.499098231322220
        sage: inverse_jacobi('sn', 0.4707504, 0.5)
        0.499999911466555
        sage: P = plot(inverse_jacobi('sn', x, 0.5), 0, 1)
    """
    if kind == 'nd':
        return inverse_jacobi_nd(x, m, **kwargs)
    elif kind == 'ns':
        return inverse_jacobi_ns(x, m, **kwargs)
    elif kind == 'nc':
        return inverse_jacobi_nc(x, m, **kwargs)
    elif kind == 'dn':
        return inverse_jacobi_dn(x, m, **kwargs)
    elif kind == 'ds':
        return inverse_jacobi_ds(x, m, **kwargs)
    elif kind == 'dc':
        return inverse_jacobi_dc(x, m, **kwargs)
    elif kind == 'sn':
        return inverse_jacobi_sn(x, m, **kwargs)
    elif kind == 'sd':
        return inverse_jacobi_sd(x, m, **kwargs)
    elif kind == 'sc':
        return inverse_jacobi_sc(x, m, **kwargs)
    elif kind == 'cn':
        return inverse_jacobi_cn(x, m, **kwargs)
    elif kind == 'cd':
        return inverse_jacobi_cd(x, m, **kwargs)
    elif kind == 'cs':
        return inverse_jacobi_cs(x, m, **kwargs)
    else:
        raise ValueError("kind must be one of 'nd', 'ns', 'nc', 'dn', "
                         "'ds', 'dc', 'sn', 'sd', 'sc', 'cn', 'cd', 'cs'.")

class JacobiAmplitude(BuiltinFunction):
    r"""
    The Jacobi amplitude function
    `\operatorname{am}(x|m) = \int_0^x \operatorname{dn}(t|m) dt` for
    `-K(m) \leq x \leq K(m)`, `F(\operatorname{am}(x|m)|m) = x`.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.functions.jacobi import JacobiAmplitude
            sage: JacobiAmplitude()
            jacobi_am
        """
        BuiltinFunction.__init__(self, name='jacobi_am', nargs=2,
                                 conversions=dict(maple='JacobiAM',
                                                  mathematica=
                                                  'JacobiAmplitude'),
                                 evalf_params_first=False)

    def _eval_(self, x, m):
        r"""
        TESTS::

            sage: jacobi_am(x, 0)
            x
            sage: jacobi_am(0, x)
            0
            sage: jacobi_am(3, 4.)
            -0.339059208303591
        """
        if m == 0:
            return x
        elif x == 0:
            return Integer(0)
        return

    def _evalf_(self, x, m, parent, algorithm=None):
        r"""
        TESTS::

            sage: jacobi_am(1, 2).n(100)
            0.73704379494724574105101929735
        """
        return utils.call(jacobi_am_f, x, m, parent=parent)

    def _derivative_(self, x, m, diff_param):
        r"""
        TESTS::

            sage: diff(jacobi_am(x, 3), x)
            jacobi_dn(x, 3)
            sage: diff(jacobi_am(3, x), x)
            -1/2*(x*jacobi_cn(3, x)*jacobi_sn(3, x) -...
            (3*x + elliptic_e(jacobi_am(3, x), x) - 3)*jacobi_dn(3, x))/((x - 1)*x)
        """
        if diff_param == 0:
            return jacobi_dn(x, m)
        elif diff_param == 1:
            return (((Integer(-1) + m) * x + elliptic_e(jacobi_am(x, m), m)) *
                    jacobi('dn', x, m) - m * jacobi('cn', x, m) *
                    jacobi('sn', x, m)) / (Integer(2) * (Integer(-1) + m) * m)

    def _latex_(self):
        r"""
        TESTS::

            sage: latex(jacobi_am)
            \operatorname{am}
        """
        return r"\operatorname{am}"

    def _print_latex_(self, x, m):
        r"""
        TESTS::

            sage: latex(jacobi_am(3,x))
            \operatorname{am}\left(3\middle|x\right)
        """
        return r"\operatorname{{am}}\left({}\middle|{}\right)".format(latex(x),
                                                                      latex(m))


jacobi_am = JacobiAmplitude()


def inverse_jacobi_f(kind, x, m):
    r"""
    Internal function for numerical evaluation of a continuous complex branch
    of each inverse Jacobi function, as described in [Tee1997]_. Only accepts
    real arguments.

    TESTS::

        sage: from mpmath import ellipfun, chop
        sage: from sage.functions.jacobi import inverse_jacobi_f

        sage: chop(ellipfun('sn', inverse_jacobi_f('sn', 0.6, 0), 0))
        mpf('0.59999999999999998')
        sage: chop(ellipfun('sn', inverse_jacobi_f('sn', 0.6, 1), 1))
        mpf('0.59999999999999998')
        sage: chop(ellipfun('sn', inverse_jacobi_f('sn', 0, -3), -3))
        mpf('0.0')
        sage: chop(ellipfun('sn', inverse_jacobi_f('sn', -1, 4), 4))
        mpf('-1.0')
        sage: chop(ellipfun('sn', inverse_jacobi_f('sn', 0.3, 4), 4))
        mpf('0.29999999999999999')
        sage: chop(ellipfun('sn', inverse_jacobi_f('sn', 0.8, 4), 4))
        mpf('0.80000000000000004')

        sage: chop(ellipfun('ns', inverse_jacobi_f('ns', 0.8, 0), 0))
        mpf('0.80000000000000004')
        sage: chop(ellipfun('ns', inverse_jacobi_f('ns', -0.7, 1), 1))
        mpf('-0.69999999999999996')
        sage: chop(ellipfun('ns', inverse_jacobi_f('ns', 0.01, 2), 2))
        mpf('0.01')
        sage: chop(ellipfun('ns', inverse_jacobi_f('ns', 0, 2), 2))
        mpf('0.0')
        sage: chop(ellipfun('ns', inverse_jacobi_f('ns', -10, 6), 6))
        mpf('-10.0')

        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', -10, 0), 0))
        mpf('-9.9999999999999982')
        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', 50, 1), 1))
        mpf('50.000000000000071')
        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', 1, 5), 5))
        mpf('1.0')
        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', 0.5, -5), -5))
        mpf('0.5')
        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', -0.75, -15), -15))
        mpf('-0.75000000000000022')
        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', 10, 0.8), 0.8))
        mpf('9.9999999999999982')
        sage: chop(ellipfun('cn', inverse_jacobi_f('cn', -2, 0.9), 0.9))
        mpf('-2.0')

        sage: chop(ellipfun('nc', inverse_jacobi_f('nc', -4, 0), 0))
        mpf('-3.9999999999999987')
        sage: chop(ellipfun('nc', inverse_jacobi_f('nc', 7, 1), 1))
        mpf('7.0000000000000009')
        sage: chop(ellipfun('nc', inverse_jacobi_f('nc', 7, 3), 3))
        mpf('7.0')
        sage: chop(ellipfun('nc', inverse_jacobi_f('nc', 0, 2), 2))
        mpf('0.0')
        sage: chop(ellipfun('nc', inverse_jacobi_f('nc', -18, -4), -4))
        mpf('-17.999999999999925')

        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', -0.3, 1), 1))
        mpf('-0.29999999999999999')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', 1, -1), -1))
        mpf('1.0')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', 0.8, 0.5), 0.5))
        mpf('0.80000000000000004')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', 5, -4), -4))
        mpf('5.0')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', 0.4, 0.5), 0.5))
        mpf('0.40000000000000002')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', -0.4, 0.5), 0.5))
        mpf('-0.40000000000000002')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', -0.9, 0.5), 0.5))
        mpf('-0.90000000000000002')
        sage: chop(ellipfun('dn', inverse_jacobi_f('dn', -1.9, 0.2), 0.2))
        mpf('-1.8999999999999999')

        sage: chop(ellipfun('nd', inverse_jacobi_f('nd', -1.9, 1), 1))
        mpf('-1.8999999999999999')
        sage: chop(ellipfun('nd', inverse_jacobi_f('nd', 1, -1), -1))
        mpf('1.0')
        sage: chop(ellipfun('nd', inverse_jacobi_f('nd', 11, -6), -6))
        mpf('11.0')
        sage: chop(ellipfun('nd', inverse_jacobi_f('nd', 0, 8), 8))
        mpf('0.0')
        sage: chop(ellipfun('nd', inverse_jacobi_f('nd', -3, 0.8), 0.8))
        mpf('-2.9999999999999996')

        sage: chop(ellipfun('sc', inverse_jacobi_f('sc', -3, 0), 0))
        mpf('-3.0')
        sage: chop(ellipfun('sc', inverse_jacobi_f('sc', 2, 1), 1))
        mpf('2.0')
        sage: chop(ellipfun('sc', inverse_jacobi_f('sc', 0, 9), 9))
        mpf('0.0')
        sage: chop(ellipfun('sc', inverse_jacobi_f('sc', -7, 3), 3))
        mpf('-7.0')

        sage: chop(ellipfun('cs', inverse_jacobi_f('cs', -7, 0), 0))
        mpf('-6.9999999999999991')
        sage: chop(ellipfun('cs', inverse_jacobi_f('cs', 8, 1), 1))
        mpf('8.0')
        sage: chop(ellipfun('cs', inverse_jacobi_f('cs', 2, 6), 6))
        mpf('2.0')
        sage: chop(ellipfun('cs', inverse_jacobi_f('cs', 0, 4), 4))
        mpf('0.0')
        sage: chop(ellipfun('cs', inverse_jacobi_f('cs', -6, 8), 8))
        mpf('-6.0000000000000018')

        sage: chop(ellipfun('cd', inverse_jacobi_f('cd', -6, 0), 0))
        mpf('-6.0000000000000009')
        sage: chop(ellipfun('cd', inverse_jacobi_f('cd', 1, 3), 3))
        mpf('1.0')
        sage: chop(ellipfun('cd', inverse_jacobi_f('cd', 6, 8), 8))
        mpf('6.0000000000000027')

        sage: chop(ellipfun('dc', inverse_jacobi_f('dc', 5, 0), 0))
        mpf('5.0000000000000018')
        sage: chop(ellipfun('dc', inverse_jacobi_f('dc', -4, 2), 2))
        mpf('-4.0000000000000018')

        sage: chop(ellipfun('sd', inverse_jacobi_f('sd', -4, 0), 0))
        mpf('-3.9999999999999991')
        sage: chop(ellipfun('sd', inverse_jacobi_f('sd', 7, 1), 1))
        mpf('7.0')
        sage: chop(ellipfun('sd', inverse_jacobi_f('sd', 0, 9), 9))
        mpf('0.0')
        sage: chop(ellipfun('sd', inverse_jacobi_f('sd', 8, 0.8), 0.8))
        mpf('7.9999999999999991')

        sage: chop(ellipfun('ds', inverse_jacobi_f('ds', 4, 0.25), 0.25))
        mpf('4.0')
    """
    from mpmath import mp

    ctx = mp
    prec = ctx.prec
    try:
        x = ctx.convert(x)
        m = ctx.convert(m)
        if not isinstance(x, ctx.mpf) or not isinstance(x, ctx.mpf):
            raise ValueError('arguments must be real')
        if kind == 'sn':
            if m == 0:
                return ctx.asin(x)
            elif m == 1:
                return ctx.atanh(x)
            elif x == 0:
                return ctx.zero
            sign = ctx.sign(x)  # sn is odd in x, so operate with abs(x) and
            x = abs(x)          # include the sign at the end
            if x <= 1:
                ctx.prec += 10
                phi = ctx.asin(x)
                return sign * ctx.ellipf(phi, m)
            elif x <= 1 / ctx.sqrt(m):
                K = ctx.ellipk(m)
                ctx.prec += 10
                xpn2 = x ** (-2)
                m1 = 1 - m
                ctx.prec += 10
                omxpn2 = 1 - xpn2
                ctx.prec += 10
                omxpn2dm1 = omxpn2 / m1
                ctx.prec += 10
                phi = ctx.asin(omxpn2dm1.sqrt())
                return sign * ctx.mpc(K, ctx.ellipf(phi, m1))
            else:
                ctx.prec += 10
                m1 = 1 - m
                K_prime = ctx.ellipk(m1)
                sqrtm = ctx.sqrt(m)
                ctx.prec += 10
                xsqrtm = x * sqrtm
                ctx.prec += 10
                phi = ctx.asin(1 / xsqrtm)
                ctx.prec += 10
                return sign * ctx.mpc(ctx.ellipf(phi, m), K_prime)
        if kind == 'ns':
            if m == 0:
                return ctx.acsc(x)
            elif m == 1:
                return ctx.acoth(x)
            elif x > 0:
                ctx.prec += 10
                return inverse_jacobi_f('sn', 1 / x, m)
            elif x == 0:
                ctx.prec += 10
                return ctx.j * ctx.ellipk(1 - m)
            else:
                ctx.prec += 10
                K_prime = ctx.ellipk(1 - m)
                odx = 1 / x
                ctx.prec += 10
                arcsnodx = inverse_jacobi_f('sn', odx, m)
                itK_prime = ctx.j * 2 * K_prime
                ctx.prec += 10
                return arcsnodx + itK_prime
        if kind == 'cn':
            if m == 0:
                return ctx.acos(x)
            elif m == 1:
                return ctx.asech(x)
            elif x == 1:
                return ctx.zero
            elif 0 <= x < 1:
                ctx.prec += 10
                x2 = x ** 2
                ctx.prec += 10
                osx2 = 1 - x2
                ctx.prec += 10
                return ctx.ellipf(ctx.asin(ctx.sqrt(osx2)), m)
            elif -1 <= x < 0:
                K = ctx.ellipk(m)
                ctx.prec += 10
                x2 = x ** 2
                ctx.prec += 10
                osx2 = 1 - x2
                ctx.prec += 10
                return (2 * K) - ctx.ellipf(ctx.asin(ctx.sqrt(osx2)), m)
            elif x > 1:
                ctx.prec += 10
                m1 = 1 - m
                xn2 = x ** (-2)
                ctx.prec += 10
                osx2 = 1 - xn2
                ctx.prec += 10
                return ctx.j * ctx.ellipf(ctx.asin(ctx.sqrt(osx2)), m1)
            elif x < -1:
                K = ctx.ellipk(m)
                ctx.prec += 10
                m1 = 1 - m
                xn2 = x ** (-2)
                tK = 2 * K
                ctx.prec += 10
                osx2 = 1 - xn2
                ctx.prec += 10
                phi = ctx.asin(ctx.sqrt(osx2))
                ctx.prec += 10
                return tK - ctx.j * ctx.ellipf(phi, m1)
        if kind == 'nc':
            if m == 0:
                return ctx.asec(x)
            elif m == 1:
                return ctx.acosh(x)
            elif x == 1:
                return ctx.zero
            elif x > 0:
                ctx.prec += 10
                return inverse_jacobi_f('cn', 1 / x, m)
            elif x == 0:
                ctx.prec += 10
                return ctx.j * ctx.ellipk(1 - m)
            else:
                K = ctx.ellipk(m)
                ctx.prec += 10
                K_prime = ctx.ellipk(1 - m)
                odx = 1 / x
                ctx.prec += 10
                arccnodx = inverse_jacobi_f('cn', odx, m)
                tK = 2 * K
                ctx.prec += 10
                return arccnodx - tK + ctx.j * 2 * K_prime
        if kind == 'dn':
            if x == 1:
                return ctx.zero
            if not m <= 1:
                raise ValueError('m must be <= 1')
            if m == 1:
                return ctx.asech(x)
            ctx.prec += 10
            m1 = 1 - m
            sqrtm1 = ctx.sqrt(m1)
            if sqrtm1 <= x < 1:
                ctx.prec += 10
                x2 = x ** 2
                ctx.prec += 10
                osx2 = 1 - x2
                ctx.prec += 10
                osx2dm = osx2 / m
                ctx.prec += 10
                return ctx.ellipf(ctx.asin(ctx.sqrt(osx2dm)), m)
            elif x > 1:
                ctx.prec += 10
                xn2 = x ** (-2)
                ctx.prec += 10
                osxn2 = 1 - xn2
                m1xn2 = m1 * xn2
                ctx.prec += 10
                osm1xn2 = 1 - m1xn2
                ctx.prec += 10
                sqrtosxn2dosm1xn2 = ctx.sqrt(osxn2 / osm1xn2)
                ctx.prec += 10
                return ctx.j * ctx.ellipf(ctx.asin(sqrtosxn2dosm1xn2), m1)
            elif 0 <= x < sqrtm1:
                K = ctx.ellipk(m)
                ctx.prec += 10
                x2 = x ** 2
                ctx.prec += 10
                x2dm1 = x2 / m1
                osx2 = 1 - x2
                ctx.prec += 10
                osx2dm1 = 1 - x2dm1
                ctx.prec += 10
                osx2dm1dosx2 = osx2dm1 / osx2
                ctx.prec += 10
                sqrtall = ctx.sqrt(osx2dm1dosx2)
                ctx.prec += 10
                phi = ctx.asin(sqrtall)
                ctx.prec += 10
                return K + ctx.j * ctx.ellipf(phi, m1)
            elif -sqrtm1 <= x < 0:
                K = ctx.ellipk(m)
                K_prime = ctx.ellipk(m1)
                ctx.prec += 10
                tK_prime = 2 * K_prime
                x2 = x ** 2
                ctx.prec += 10
                x2dm1 = x2 / m1
                osx2 = 1 - x2
                ctx.prec += 10
                osx2dm1 = 1 - x2dm1
                ctx.prec += 10
                osx2dm1dosx2 = osx2dm1 / osx2
                ctx.prec += 10
                sqrtall = ctx.sqrt(osx2dm1dosx2)
                ctx.prec += 10
                phi = ctx.asin(sqrtall)
                ctx.prec += 10
                return K + ctx.j * (tK_prime - ctx.ellipf(phi, m1))
            elif -1 <= x < -sqrtm1:
                K = ctx.ellipk(m)
                K_prime = ctx.ellipk(m1)
                ctx.prec += 10
                x2 = x ** 2
                tK = 2 * K
                # Note that the factor of 2 is missing in the reference
                # (formula (81)), probably mistakenly so
                tK_prime = 2 * K_prime
                ctx.prec += 10
                osx2 = 1 - x2
                ctx.prec += 10
                osx2dm = osx2 / m
                sqrtall = ctx.sqrt(osx2dm)
                ctx.prec += 10
                phi = ctx.asin(sqrtall)
                ctx.prec += 10
                return (tK - ctx.ellipf(phi, m)) + (ctx.j * tK_prime)
            elif x < -1:
                K = ctx.ellipk(m)
                K_prime = ctx.ellipk(m1)
                ctx.prec += 10
                tK = 2 * K
                tK_prime = 2 * K_prime
                xn2 = x ** (-2)
                ctx.prec += 10
                osxn2 = 1 - xn2
                m1xn2 = m1 * xn2
                ctx.prec += 10
                osm1xn2 = 1 - m1xn2
                ctx.prec += 10
                sqrtosxn2dosm1xn2 = ctx.sqrt(osxn2 / osm1xn2)
                ctx.prec += 10
                phi = ctx.asin(sqrtosxn2dosm1xn2)
                ctx.prec += 10
                return tK + ctx.j * (tK_prime - ctx.ellipf(phi, m1))
        if kind == 'nd':
            if m == 1:
                return ctx.acosh(x)
            elif x == 1:
                return ctx.zero
            elif x > 0:
                ctx.prec += 10
                return inverse_jacobi_f('dn', 1 / x, m)
            elif x == 0:
                ctx.prec += 10
                return ctx.j * ctx.ellipk(1 - m)
            else:
                K = ctx.ellipk(m)
                ctx.prec += 10
                tK = 2 * K
                ctx.prec += 10
                return inverse_jacobi_f('dn', 1 / x, m) - tK
        if kind == 'sc':
            if m == 0:
                return ctx.atan(x)
            elif m == 1:
                return ctx.asinh(x)
            elif x == 0:
                return ctx.zero
            else:
                ctx.prec += 10
                atanx = ctx.atan(x)
                return ctx.ellipf(atanx, m)
        if kind == 'cs':
            if m == 0:
                return ctx.acot(x)
            elif m == 1:
                return ctx.acsch(x)
            elif x > 0:
                ctx.prec += 10
                odx = 1 / x
                ctx.prec += 10
                return ctx.ellipf(ctx.atan(odx), m)
            elif x == 0:
                return ctx.ellipk(m)
            else:
                K = ctx.ellipk(m)
                ctx.prec += 10
                odx = 1 / x
                ctx.prec += 10
                phi = ctx.atan(odx)
                ctx.prec += 10
                return ctx.ellipf(phi, m) + (2 * K)
        if kind == 'cd':
            if m == 0:
                return ctx.acos(x)
            elif x == 1:
                return ctx.zero
            else:
                K = ctx.ellipk(m)
                ctx.prec += 10
                return inverse_jacobi_f('sn', x, m) - K
        if kind == 'dc':
            if m == 0:
                return ctx.asec(x)
            K = ctx.ellipk(m)
            ctx.prec += 10
            return inverse_jacobi_f('ns', x, m) - K
        if kind == 'sd':
            if m == 0:
                return ctx.asin(x)
            elif m == 1:
                return ctx.asinh(x)
            elif x == 0:
                return ctx.zero
            else:
                if m > 1:
                    raise ValueError('m must be <= 1')
                K = ctx.ellipk(m)
                ctx.prec += 10
                m1 = 1 - m
                ctx.prec += 10
                sqrtm1 = ctx.sqrt(m1)
                ctx.prec += 10
                xsqrtm1 = x * sqrtm1
                ctx.prec += 10
                return inverse_jacobi_f('cn', xsqrtm1, m) + K
        if kind == 'ds':
            if m == 0:
                return ctx.acsc(x)
            elif m == 1:
                return ctx.acsch(x)
            else:
                if m > 1:
                    raise ValueError('m must be <= 1')
                K = ctx.ellipk(m)
                ctx.prec += 10
                m1 = 1 - m
                ctx.prec += 10
                sqrtm1 = ctx.sqrt(m1)
                ctx.prec += 10
                xdsqrtm1 = x / sqrtm1
                ctx.prec += 10
                return inverse_jacobi_f('nc', xdsqrtm1, m) + K
    finally:
        ctx.prec = prec


def jacobi_am_f(x, m):
    r"""
    Internal function for numeric evaluation of the Jacobi amplitude function
    for real arguments. Procedure described in [Eh2013]_.

    TESTS::

        sage: from mpmath import ellipf
        sage: from sage.functions.jacobi import jacobi_am_f
        sage: ellipf(jacobi_am_f(0.5, 1), 1)
        mpf('0.5')
        sage: ellipf(jacobi_am(3, 0.3), 0.3)
        mpf('3.0')
        sage: ellipf(jacobi_am_f(2, -0.5), -0.5)
        mpf('2.0')
        sage: jacobi_am_f(2, -0.5)
        mpf('2.2680930777934176')
        sage: jacobi_am_f(-2, -0.5)
        mpf('-2.2680930777934176')
        sage: jacobi_am_f(-3, 2)
        mpf('0.36067407399586108')
    """
    from mpmath import mp

    ctx = mp
    prec = ctx.prec
    try:
        x = ctx.convert(x)
        m = ctx.convert(m)
        if not isinstance(x, ctx.mpf) or not isinstance(m, ctx.mpf):
            raise ValueError('arguments must be real')
        if abs(m) == 1:
            # gd(x)
            ctx.prec += 10
            tanhx = ctx.tanh(x)
            ctx.prec += 10
            return ctx.asin(tanhx)
        elif abs(m) > 1:
            ctx.prec += 10
            # Real values needed for atan2; as per "Handbook of Elliptic
            # Integrals for Engineers and Scientists" 121.02, sn is real for
            # real x. The imaginary components can thus be safely discarded.
            snx = ctx.ellipfun('sn', x, m).real
            cnx = ctx.ellipfun('cn', x, m).real
            ctx.prec += 10
            return ctx.atan2(snx, cnx)
        else:
            ctx.prec += 10
            K = ctx.ellipk(m)
            if abs(x) <= K:
                snx = ctx.ellipfun('sn', x, m).real
                cnx = ctx.ellipfun('cn', x, m).real
                ctx.prec += 10
                return ctx.atan2(snx, cnx)
            else:
                # Do argument reduction on x to end up with z = x - 2nK, with
                # abs(z) <= K
                ctx.prec += 10
                tK = 2 * K
                ctx.prec += 10
                n = ctx.floor(x / tK)
                ctx.prec += 10
                tnK = n * tK
                npi = n * ctx.pi()
                ctx.prec += 10
                z = x - tnK
                ctx.prec += 10
                # z (and therefore sn(z, m) and cn(z, m)) is real because K(m)
                # is real for abs(m) <= 1.
                snz = ctx.ellipfun('sn', z, m).real
                cnz = ctx.ellipfun('cn', z, m).real
                ctx.prec += 10
                return ctx.atan2(snz, cnz) + npi
    finally:
        ctx.prec = prec
