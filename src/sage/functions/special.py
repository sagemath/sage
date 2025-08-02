r"""
Miscellaneous special functions

This module provides easy access to many of Maxima and PARI's
special functions.

Maxima's special functions package (which includes spherical
harmonic functions, spherical Bessel functions (of the 1st and 2nd
kind), and spherical Hankel functions (of the 1st and 2nd kind))
was written by Barton Willis of the University of Nebraska at
Kearney. It is released under the terms of the General Public
License (GPL).

Support for elliptic functions and integrals was written by Raymond
Toy. It is placed under the terms of the General Public License
(GPL) that governs the distribution of Maxima.

Next, we summarize some of the properties of the functions
implemented here.

- **Spherical harmonics**: Laplace's equation in spherical coordinates is:

    .. MATH::

        \frac{1}{r^2} \frac{\partial}{\partial r}
        \left( r^2 \frac{\partial f}{\partial r} \right) +
        \frac{1}{r^2\sin\theta} \frac{\partial}{\partial \theta}
        \left( \sin\theta \frac{\partial f}{\partial \theta} \right) +
        \frac{1}{r^2\sin^2\theta} \frac{\partial^2 f}{\partial \varphi^2} = 0.

  Note that the spherical coordinates `\theta` and `\varphi` are defined here
  as follows: `\theta` is the colatitude or polar angle, ranging from
  `0\leq\theta\leq\pi` and `\varphi` the azimuth or longitude, ranging from
  `0\leq\varphi<2\pi`.

  The general solution which remains finite towards infinity is a linear
  combination of functions of the form

    .. MATH::

        r^{-1-\ell} \cos (m \varphi) P_\ell^m (\cos{\theta} )

   and

    .. MATH::

        r^{-1-\ell} \sin (m \varphi) P_\ell^m (\cos{\theta} )

  where `P_\ell^m` are the associated Legendre polynomials
  (cf. :class:`~sage.functions.orthogonal_polys.Func_assoc_legendre_P`),
  and with integer parameters `\ell \ge 0` and `m` from `0` to `\ell`. Put in
  another way, the solutions with integer parameters `\ell \ge 0` and
  `- \ell\leq m\leq \ell`, can be written as linear combinations of:

    .. MATH::

        U_{\ell,m}(r,\theta , \varphi ) =
        r^{-1-\ell} Y_\ell^m( \theta , \varphi )

  where the functions `Y` are the spherical harmonic functions with
  parameters `\ell`, `m`, which can be written as:

    .. MATH::

        Y_\ell^m( \theta , \varphi ) =
        \sqrt{ \frac{(2\ell+1)}{4\pi} \frac{(\ell-m)!}{(\ell+m)!} }
        \, e^{i m \varphi } \, P_\ell^m ( \cos{\theta} ) .

  The spherical harmonics obey the normalisation condition

    .. MATH::

        \int_{\theta=0}^\pi\int_{\varphi=0}^{2\pi}
        Y_\ell^mY_{\ell'}^{m'*}\,d\Omega =
        \delta_{\ell\ell'}\delta_{mm'}\quad\quad d\Omega =
        \sin\theta\,d\varphi\,d\theta .

- The **incomplete elliptic integrals** (of the first kind, etc.) are:

    .. MATH::

        \begin{array}{c}
        \displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2}}\, dx,\\
        \displaystyle\int_0^\phi \sqrt{1 - m\sin(x)^2}\, dx,\\
        \displaystyle\int_0^\phi \frac{\sqrt{1-mt^2}}{\sqrt(1 - t^2)}\, dx,\\
        \displaystyle\int_0^\phi
        \frac{1}{\sqrt{1 - m\sin(x)^2\sqrt{1 - n\sin(x)^2}}}\, dx,
        \end{array}

  and the complete ones are obtained by taking `\phi =\pi/2`.

.. WARNING::

    SciPy's versions are poorly documented and seem less accurate than the
    Maxima and PARI versions. Typically they are limited by hardware floats
    precision.

REFERENCES:

- Abramowitz and Stegun: *Handbook of Mathematical Functions* [AS1964]_

- :wikipedia:`Spherical_harmonics`

- :wikipedia:`Helmholtz_equation`

- `Online Encyclopedia of Special Functions
  <http://algo.inria.fr/esf/index.html>`_

AUTHORS:

- David Joyner (2006-13-06): initial version

- David Joyner (2006-30-10): bug fixes to pari wrappers of Bessel
  functions, hypergeometric_U

- William Stein (2008-02): Impose some sanity checks.

- David Joyner (2008-02-16): optional calls to scipy and replace all ``#random`` by ``...``

- David Joyner (2008-04-23): addition of elliptic integrals

- Eviatar Bach (2013): making elliptic integrals symbolic

- Eric Gourgoulhon (2022): add Condon-Shortley phase to spherical harmonics
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc

from sage.misc.functional import sqrt
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.element import parent as s_parent
from sage.symbolic.function import BuiltinFunction

lazy_import('sage.functions.jacobi', 'jacobi_am_f')
lazy_import('sage.functions.log', ['exp'])
lazy_import('sage.functions.trig', ['sin', 'cot'])

lazy_import('sage.misc.latex', 'latex')

lazy_import('sage.symbolic.constants', ['I', 'pi'])

lazy_import('sage.libs.mpmath.utils', 'call', as_='_mpmath_utils_call')
lazy_import('mpmath',
            ['spherharm', 'ellipe', 'ellipf', 'ellipk', 'ellippi'],
            as_=['_mpmath_spherharm', '_mpmath_ellipe', '_mpmath_ellipf',
                 '_mpmath_ellipk', '_mpmath_ellippi'])


class SphericalHarmonic(BuiltinFunction):
    r"""
    Return the spherical harmonic function `Y_n^m(\theta, \varphi)`.

    For integers `n > -1`, `|m| \leq n`, simplification is done automatically.
    Numeric evaluation is supported for complex `n` and `m`.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: x, y = var('x, y')
        sage: spherical_harmonic(3, 2, x, y)
        1/8*sqrt(30)*sqrt(7)*cos(x)*e^(2*I*y)*sin(x)^2/sqrt(pi)
        sage: spherical_harmonic(3, 2, 1, 2)
        1/8*sqrt(30)*sqrt(7)*cos(1)*e^(4*I)*sin(1)^2/sqrt(pi)
        sage: spherical_harmonic(3 + I, 2., 1, 2)
        -0.351154337307488 - 0.415562233975369*I
        sage: latex(spherical_harmonic(3, 2, x, y, hold=True))
        Y_{3}^{2}\left(x, y\right)
        sage: spherical_harmonic(1, 2, x, y)
        0

    The degree `n` and the order `m` can be symbolic::

        sage: # needs sage.symbolic
        sage: n, m = var('n m')
        sage: spherical_harmonic(n, m, x, y)
        spherical_harmonic(n, m, x, y)
        sage: latex(spherical_harmonic(n, m, x, y))
        Y_{n}^{m}\left(x, y\right)
        sage: diff(spherical_harmonic(n, m, x, y), x)
        m*cot(x)*spherical_harmonic(n, m, x, y)
         + sqrt(-(m + n + 1)*(m - n))*e^(-I*y)*spherical_harmonic(n, m + 1, x, y)
        sage: diff(spherical_harmonic(n, m, x, y), y)
        I*m*spherical_harmonic(n, m, x, y)

    The convention regarding the Condon-Shortley phase `(-1)^m` is the same
    as for SymPy's spherical harmonics and :wikipedia:`Spherical_harmonics`::

        sage: # needs sage.symbolic
        sage: spherical_harmonic(1, 1, x, y)
        -1/4*sqrt(3)*sqrt(2)*e^(I*y)*sin(x)/sqrt(pi)
        sage: from sympy import Ynm                                                     # needs sympy
        sage: Ynm(1, 1, x, y).expand(func=True)                                         # needs sympy
        -sqrt(6)*exp(I*y)*sin(x)/(4*sqrt(pi))
        sage: spherical_harmonic(1, 1, x, y) - Ynm(1, 1, x, y)                          # needs sympy
        0

    It also agrees with SciPy's spherical harmonics::

        sage: spherical_harmonic(1, 1, pi/2, pi).n()  # abs tol 1e-14                   # needs sage.symbolic
        0.345494149471335
        sage: import numpy as np                                                        # needs scipy
        sage: if int(np.version.short_version[0]) > 1:                                  # needs scipy
        ....:     _ = np.set_printoptions(legacy="1.25")                                # needs scipy
        sage: import scipy.version
        sage: if scipy.version.version < '1.15.0':
        ....:     from scipy.special import sph_harm # NB: arguments x and y are swapped   # needs scipy
        ....:     sph_harm(1, 1, pi.n(), (pi/2).n())  # abs tol 1e-14                   # needs scipy sage.symbolic
        ....: else:
        ....:     from scipy.special import sph_harm_y                                  # needs scipy
        ....:     sph_harm_y(1, 1, (pi/2).n(), pi.n()).item()  # abs tol 1e-9           # needs scipy sage.symbolic
        (0.3454941494713355-4.231083042742082e-17j)

    Note that this convention differs from the one in Maxima, as revealed by
    the sign difference for odd values of `m`::

        sage: maxima.spherical_harmonic(1, 1, x, y).sage()                              # needs sage.symbolic
        1/2*sqrt(3/2)*e^(I*y)*sin(x)/sqrt(pi)

    It follows that, contrary to Maxima, SageMath uses the same sign convention
    for spherical harmonics as SymPy, SciPy, Mathematica and
    :wikipedia:`Table_of_spherical_harmonics`.

    REFERENCES:

    - :wikipedia:`Spherical_harmonics`
    """
    def __init__(self):
        r"""
        TESTS::

            sage: n, m, theta, phi = var('n m theta phi')                               # needs sage.symbolic
            sage: spherical_harmonic(n, m, theta, phi)._sympy_()                        # needs sympy sage.symbolic
            Ynm(n, m, theta, phi)
        """
        BuiltinFunction.__init__(self, 'spherical_harmonic', nargs=4,
                                 conversions=dict(
                                     maple='SphericalY',
                                     mathematica='SphericalHarmonicY',
                                     maxima='spherical_harmonic',
                                     sympy='Ynm'))

    def _eval_(self, n, m, theta, phi, **kwargs):
        r"""
        TESTS::

            sage: # needs sage.symbolic
            sage: x, y = var('x y')
            sage: spherical_harmonic(1, 2, x, y)
            0
            sage: spherical_harmonic(1, -2, x, y)
            0
            sage: spherical_harmonic(1/2, 2, x, y)
            spherical_harmonic(1/2, 2, x, y)
            sage: spherical_harmonic(3, 2, x, y)
            1/8*sqrt(30)*sqrt(7)*cos(x)*e^(2*I*y)*sin(x)^2/sqrt(pi)
            sage: spherical_harmonic(3, 2, 1, 2)
            1/8*sqrt(30)*sqrt(7)*cos(1)*e^(4*I)*sin(1)^2/sqrt(pi)
            sage: spherical_harmonic(3 + I, 2., 1, 2)
            -0.351154337307488 - 0.415562233975369*I

        Check that :issue:`20939` is fixed::

            sage: ex = spherical_harmonic(3, 2, 1, 2*pi/3)                              # needs sage.symbolic
            sage: QQbar(ex * sqrt(pi)/cos(1)/sin(1)^2).minpoly()                        # needs sage.rings.number_field sage.symbolic
            x^4 + 105/32*x^2 + 11025/1024

        Check whether Sage yields correct results compared to Maxima,
        up to the Condon-Shortley phase factor `(-1)^m`
        (see :issue:`25034` and :issue:`33117`)::

            sage: # needs sage.symbolic
            sage: spherical_harmonic(1, 1, pi/3, pi/6).n()  # abs tol 1e-14
            -0.259120612103502 - 0.149603355150537*I
            sage: maxima.spherical_harmonic(1, 1, pi/3, pi/6).n()  # abs tol 1e-14
            0.259120612103502 + 0.149603355150537*I
            sage: spherical_harmonic(1, -1, pi/3, pi/6).n()  # abs tol 1e-14
            0.259120612103502 - 0.149603355150537*I
            sage: maxima.spherical_harmonic(1, -1, pi/3, pi/6).n()  # abs tol 1e-14
            -0.259120612103502 + 0.149603355150537*I

        Check that :issue:`33501` is fixed::

            sage: spherical_harmonic(2, 1, x, y)                                        # needs sage.symbolic
            -1/4*sqrt(6)*sqrt(5)*cos(x)*e^(I*y)*sin(x)/sqrt(pi)
            sage: spherical_harmonic(5, -3, x, y)                                       # needs sage.symbolic
            -1/32*(9*sqrt(385)*sin(x)^4 - 8*sqrt(385)*sin(x)^2)*e^(-3*I*y)*sin(x)/sqrt(pi)
        """
        if n in ZZ and m in ZZ and n > -1:
            if abs(m) > n:
                return ZZ(0)
            if m == 0 and theta.is_zero():
                return sqrt((2*n+1)/4/pi)
            from sage.arith.misc import factorial
            from sage.functions.trig import cos
            from sage.functions.orthogonal_polys import gen_legendre_P
            res = (sqrt(factorial(n-m) * (2*n+1) / (4*pi * factorial(n+m)))
                   * gen_legendre_P(n, m, cos(theta))
                   * exp(I*m*phi)).simplify_trig()
            res = res.substitute({sqrt(sin(theta)**2): sin(theta)})
            return res

    def _evalf_(self, n, m, theta, phi, parent, **kwds):
        r"""
        TESTS::

            sage: spherical_harmonic(3 + I, 2, 1, 2).n(100)                             # needs sage.symbolic
            -0.35115433730748836508201061672 - 0.41556223397536866209990358597*I
            sage: spherical_harmonic(I, I, I, I).n()                                    # needs sage.symbolic
            7.66678546069894 - 0.265754432549751*I

        Consistency with ``_eval_``::

            sage: d = lambda a, b: abs(spherical_harmonic(a, b, 1., 2.)
            ....:                      - spherical_harmonic(a, b, 1, 2).n())
            sage: ab = [(0, 0), (1, -1), (1, 0), (1, 1), (3, 2), (3, 3)]
            sage: all(d(a, b) < 1e-14 for a, b in ab)                                   # needs sage.symbolic
            True
        """
        return _mpmath_utils_call(_mpmath_spherharm, n, m, theta, phi, parent=parent)

    def _derivative_(self, n, m, theta, phi, diff_param):
        r"""
        TESTS::

            sage: # needs sage.symbolic
            sage: n, m, theta, phi = var('n m theta phi')
            sage: Ynm = spherical_harmonic(n, m, theta, phi)
            sage: DY_theta = Ynm.diff(theta); DY_theta
            m*cot(theta)*spherical_harmonic(n, m, theta, phi)
             + sqrt(-(m + n + 1)*(m - n))*e^(-I*phi)*spherical_harmonic(n, m + 1, theta, phi)
            sage: Ynm.diff(phi)
            I*m*spherical_harmonic(n, m, theta, phi)

        Check that :issue:`33117` is fixed::

            sage: # needs sage.symbolic
            sage: DY_theta.subs({n: 1, m: 0})
            -1/2*sqrt(3)*sin(theta)/sqrt(pi)
            sage: Ynm.subs({n: 1, m: 0}).diff(theta)
            -1/2*sqrt(3)*sin(theta)/sqrt(pi)
            sage: bool(DY_theta.subs({n: 1, m: 0}) == Ynm.subs({n: 1, m: 0}).diff(theta))
            True
            sage: bool(DY_theta.subs({n: 1, m: 1}) == Ynm.subs({n: 1, m: 1}).diff(theta))
            True
            sage: bool(DY_theta.subs({n: 1, m: -1}) == Ynm.subs({n: 1, m: -1}).diff(theta))
            True
        """
        if diff_param == 2:
            return (m * cot(theta) * spherical_harmonic(n, m, theta, phi) +
                    sqrt((n - m) * (n + m + 1)) * exp(-I * phi) *
                    spherical_harmonic(n, m + 1, theta, phi))
        if diff_param == 3:
            return I * m * spherical_harmonic(n, m, theta, phi)

        raise ValueError('only derivative with respect to theta or phi'
                         ' supported')

    def _latex_(self):
        r"""
        TESTS::

            sage: latex(spherical_harmonic)
            Y_n^m
        """
        return r"Y_n^m"

    def _print_latex_(self, n, m, theta, phi):
        r"""
        TESTS::

            sage: y = var('y')                                                          # needs sage.symbolic
            sage: latex(spherical_harmonic(3, 2, x, y, hold=True))                      # needs sage.symbolic
            Y_{3}^{2}\left(x, y\right)
        """
        return r"Y_{{{}}}^{{{}}}\left({}, {}\right)".format(
            latex(n), latex(m), latex(theta), latex(phi))


spherical_harmonic = SphericalHarmonic()


# elliptic functions and integrals

def elliptic_j(z, prec=53):
    r"""
    Return the elliptic modular `j`-function evaluated at `z`.

    INPUT:

    - ``z`` -- complex; a complex number with positive imaginary part

    - ``prec`` -- (default: 53) precision in bits for the complex field

    OUTPUT: (complex) the value of `j(z)`

    ALGORITHM:

    Calls the ``pari`` function ``ellj()``.

    AUTHOR:

    John Cremona

    EXAMPLES::

        sage: elliptic_j(CC(i))                                                         # needs sage.rings.real_mpfr
        1728.00000000000
        sage: elliptic_j(sqrt(-2.0))                                                    # needs sage.rings.complex_double
        8000.00000000000
        sage: z = ComplexField(100)(1, sqrt(11))/2                                      # needs sage.rings.real_mpfr sage.symbolic
        sage: elliptic_j(z)                                                             # needs sage.rings.real_mpfr sage.symbolic
        -32768.000...
        sage: elliptic_j(z).real().round()                                              # needs sage.rings.real_mpfr sage.symbolic
        -32768

    ::

        sage: tau = (1 + sqrt(-163))/2                                                  # needs sage.symbolic
        sage: (-elliptic_j(tau.n(100)).real().round())^(1/3)                            # needs sage.symbolic
        640320

    This example shows the need for higher precision than the default one of
    the `ComplexField`, see :issue:`28355`::

        sage: # needs sage.symbolic
        sage: -elliptic_j(tau)  # rel tol 1e-2
        2.62537412640767e17 - 732.558854258998*I
        sage: -elliptic_j(tau, 75)  # rel tol 1e-2
        2.625374126407680000000e17 - 0.0001309913593909879441262*I
        sage: -elliptic_j(tau, 100)  # rel tol 1e-2
        2.6253741264076799999999999999e17 - 1.3012822400356887122945119790e-12*I
        sage: (-elliptic_j(tau, 100).real().round())^(1/3)
        640320
    """
    CC = z.parent()
    if not isinstance(CC, sage.rings.abc.ComplexField):
        from sage.rings.complex_mpfr import ComplexField
        CC = ComplexField(prec)
        try:
            z = CC(z)
        except ValueError:
            raise ValueError("elliptic_j only defined for complex arguments.")
    from sage.libs.pari import pari
    return CC(pari(z).ellj())


# elliptic integrals

class EllipticE(BuiltinFunction):
    r"""
    Return the incomplete elliptic integral of the
    second kind:

    .. MATH::

        E(\varphi\,|\,m)=\int_0^\varphi \sqrt{1 - m\sin(x)^2}\, dx.

    EXAMPLES::

        sage: z = var("z")                                                              # needs sage.symbolic
        sage: elliptic_e(z, 1)                                                          # needs sage.symbolic
        elliptic_e(z, 1)
        sage: elliptic_e(z, 1).simplify()       # not tested                            # needs sage.symbolic
        2*round(z/pi) - sin(pi*round(z/pi) - z)
        sage: elliptic_e(z, 0)                                                          # needs sage.symbolic
        z
        sage: elliptic_e(0.5, 0.1)  # abs tol 2e-15                                     # needs mpmath
        0.498011394498832
        sage: elliptic_e(1/2, 1/10).n(200)                                              # needs sage.symbolic
        0.4980113944988315331154610406...

    .. SEEALSO::

        - Taking `\varphi = \pi/2` gives
          :func:`elliptic_ec()<sage.functions.special.EllipticEC>`.

        - Taking `\varphi = \operatorname{arc\,sin}(\operatorname{sn}(u,m))`
          gives :func:`elliptic_eu()<sage.functions.special.EllipticEU>`.

    REFERENCES:

    - :wikipedia:`Elliptic_integral#Incomplete_elliptic_integral_of_the_second_kind`

    - :wikipedia:`Jacobi_elliptic_functions`
    """
    def __init__(self):
        r"""
        TESTS::

            sage: loads(dumps(elliptic_e))
            elliptic_e
            sage: elliptic_e(x, x)._sympy_()                                            # needs sympy sage.symbolic
            elliptic_e(x, x)

        Check that :issue:`34085` is fixed::

            sage: _ = var("x y")                                                        # needs sage.symbolic
            sage: fricas(elliptic_e(x, y))                                          # optional - fricas, needs sage.symbolic
            ellipticE(sin(x),y)

        However, the conversion is only correct in the interval
        `[-\pi/2, \pi/2]`::

            sage: fricas(elliptic_e(x, y)).D(x).sage()/elliptic_e(x, y).diff(x)     # optional - fricas, needs sage.symbolic
            cos(x)/sqrt(-sin(x)^2 + 1)

        Numerically::

            sage: f = lambda x, y: elliptic_e(arcsin(x), y).subs(x=x, y=y)
            sage: g = lambda x, y: fricas.ellipticE(x, y).sage()
            sage: d = lambda x, y: f(x, y) - g(x, y)
            sage: [d(N(-pi/2 + x), y)                           # abs tol 1e-8      # optional - fricas, needs sage.symbolic
            ....:  for x in range(1, 3) for y in range(-2, 2)]
            [0.000000000000000,
             0.000000000000000,
             0.000000000000000,
             0.000000000000000,
             5.55111512312578e-17,
             0.000000000000000,
             0.000000000000000,
             0.000000000000000]
        """
        BuiltinFunction.__init__(self, 'elliptic_e', nargs=2,
                                 # Maple conversion left out since it uses
                                 # k instead of m as the second argument
                                 conversions=dict(mathematica='EllipticE',
                                                  maxima='elliptic_e',
                                                  sympy='elliptic_e',
                                                  fricas='((x,y)+->ellipticE(sin(x), y))'))

    def _eval_(self, z, m):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: z = var("z")
            sage: elliptic_e(0, x)
            0
            sage: elliptic_e(pi/2, x)
            elliptic_ec(x)
            sage: elliptic_e(z, 0)
            z
            sage: elliptic_e(z, 1)
            elliptic_e(z, 1)

        Here arccoth doesn't have 1 in its domain, so we just hold the expression::

            sage: elliptic_e(arccoth(1), x^2*e)                                         # needs sage.symbolic
            elliptic_e(+Infinity, x^2*e)
        """
        if z == 0:
            return Integer(0)
        elif z == pi / 2:
            return elliptic_ec(m)
        elif m == 0:
            return z

    def _evalf_(self, z, m, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: elliptic_e(0.5, 0.1)                                                  # needs mpmath
            0.498011394498832
            sage: elliptic_e(1/2, 1/10).n(200)                                          # needs sage.symbolic
            0.4980113944988315331154610406...
            sage: elliptic_e(I, I).n()                                                  # needs sage.symbolic
            -0.189847437084712 + 1.03209769372160*I

        TESTS:

        This gave an error in Maxima (:issue:`15046`)::

            sage: elliptic_e(2.5, 2.5)                                                  # needs mpmath
            0.535647771608740 + 1.63996015168665*I
        """
        R = parent or s_parent(z)
        return _mpmath_utils_call(_mpmath_ellipe, z, m, parent=R)

    def _derivative_(self, z, m, diff_param):
        """
        EXAMPLES::

            sage: x, z = var('x,z')                                                     # needs sage.symbolic
            sage: elliptic_e(z, x).diff(z, 1)                                           # needs sage.symbolic
            sqrt(-x*sin(z)^2 + 1)
            sage: elliptic_e(z, x).diff(x, 1)                                           # needs sage.symbolic
            1/2*(elliptic_e(z, x) - elliptic_f(z, x))/x
        """
        if diff_param == 0:
            return sqrt(Integer(1) - m * sin(z) ** Integer(2))
        elif diff_param == 1:
            return (elliptic_e(z, m) - elliptic_f(z, m)) / (Integer(2) * m)

    def _print_latex_(self, z, m):
        r"""
        EXAMPLES::

            sage: latex(elliptic_e(pi, x))                                              # needs sage.symbolic
            E(\pi\,|\,x)
        """
        return r"E(%s\,|\,%s)" % (latex(z), latex(m))


elliptic_e = EllipticE()


class EllipticEC(BuiltinFunction):
    r"""
    Return the complete elliptic integral of the second kind:

    .. MATH::

        E(m)=\int_0^{\pi/2} \sqrt{1 - m\sin(x)^2}\, dx.

    EXAMPLES::

        sage: elliptic_ec(0.1)                                                          # needs mpmath
        1.53075763689776
        sage: elliptic_ec(x).diff()                                                     # needs sage.symbolic
        1/2*(elliptic_ec(x) - elliptic_kc(x))/x

    .. SEEALSO::

        - :func:`elliptic_e()<sage.functions.special.EllipticE>`.

    REFERENCES:

    - :wikipedia:`Elliptic_integral#Complete_elliptic_integral_of_the_second_kind`
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: loads(dumps(elliptic_ec))
            elliptic_ec
            sage: elliptic_ec(x)._sympy_()                                              # needs sage.symbolic
            elliptic_e(x)

        TESTS::

            sage: fricas(elliptic_ec(x))                                        # optional - fricas, needs sage.symbolic
            ellipticE(x)

            sage: elliptic_ec(0.5)  # abs tol 1e-8                                      # needs sage.symbolic
            1.35064388104768
            sage: fricas.ellipticE(0.5).sage()  # abs tol 1e-8                  # optional - fricas, needs sage.symbolic
            1.3506438810476755025201749
        """
        BuiltinFunction.__init__(self, 'elliptic_ec', nargs=1, latex_name='E',
                                 conversions=dict(mathematica='EllipticE',
                                                  maxima='elliptic_ec',
                                                  sympy='elliptic_e',
                                                  fricas='ellipticE'))

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: elliptic_ec(0)                                                        # needs sage.symbolic
            1/2*pi
            sage: elliptic_ec(1)                                                        # needs sage.symbolic
            1
            sage: elliptic_ec(x)                                                        # needs sage.symbolic
            elliptic_ec(x)
        """
        if x == 0:
            return pi / Integer(2)
        elif x == 1:
            return Integer(1)

    def _evalf_(self, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: elliptic_ec(sqrt(2)/2).n()                                            # needs sage.symbolic
            1.23742252487318
            sage: elliptic_ec(sqrt(2)/2).n(200)                                         # needs sage.symbolic
            1.237422524873181672854746084083...
            sage: elliptic_ec(I).n()                                                    # needs sage.symbolic
            1.63241178144043 - 0.369219492375499*I
        """
        R = parent or s_parent(x)
        return _mpmath_utils_call(_mpmath_ellipe, x, parent=R)

    def _derivative_(self, x, diff_param):
        """
        EXAMPLES::

            sage: elliptic_ec(x).diff()                                                 # needs sage.symbolic
            1/2*(elliptic_ec(x) - elliptic_kc(x))/x
        """
        return (elliptic_ec(x) - elliptic_kc(x)) / (Integer(2) * x)


elliptic_ec = EllipticEC()


class EllipticEU(BuiltinFunction):
    r"""
    Return Jacobi's form of the incomplete elliptic integral of the second kind:

    .. MATH::

        E(u,m)=
        \int_0^u \mathrm{dn}(x,m)^2\, dx = \int_0^\tau
        \frac{\sqrt{1-m x^2}}{\sqrt{1-x^2}}\, dx.

    where `\tau = \mathrm{sn}(u, m)`.

    Also, ``elliptic_eu(u, m) = elliptic_e(asin(sn(u,m)),m)``.

    EXAMPLES::

        sage: elliptic_eu(0.5, 0.1)                                                     # needs mpmath
        0.496054551286597

    .. SEEALSO::

        - :func:`elliptic_e()<sage.functions.special.EllipticE>`.

    REFERENCES:

    - :wikipedia:`Elliptic_integral#Incomplete_elliptic_integral_of_the_second_kind`

    - :wikipedia:`Jacobi_elliptic_functions`
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(elliptic_eu))
            elliptic_eu
        """
        BuiltinFunction.__init__(self, 'elliptic_eu', nargs=2,
                                 conversions=dict(maxima='elliptic_eu'))

    def _eval_(self, u, m):
        """
        EXAMPLES::

            sage: elliptic_eu(1, 1)                                                     # needs sage.symbolic
            elliptic_eu(1, 1)
        """
        pass

    def _evalf_(self, u, m, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: elliptic_eu(1, 1).n()                                                 # needs sage.symbolic
            0.761594155955765
            sage: elliptic_eu(1, 1).n(200)                                              # needs sage.symbolic
            0.7615941559557648881194582...
        """
        R = parent or s_parent(u)
        return _mpmath_utils_call(elliptic_eu_f, u, m, parent=R)

    def _derivative_(self, u, m, diff_param):
        """
        EXAMPLES::

            sage: x, m = var('x,m')                                                     # needs sage.symbolic
            sage: elliptic_eu(x, m).diff(x)                                             # needs sage.symbolic
            sqrt(-m*jacobi_sn(x, m)^2 + 1)*jacobi_dn(x, m)
            sage: elliptic_eu(x, m).diff(m)                                             # needs sage.symbolic
            1/2*(elliptic_eu(x, m)
             - elliptic_f(jacobi_am(x, m), m))/m
             - 1/2*(m*jacobi_cn(x, m)*jacobi_sn(x, m)
             - (m - 1)*x
             - elliptic_eu(x, m)*jacobi_dn(x, m))*sqrt(-m*jacobi_sn(x, m)^2 + 1)/((m - 1)*m)
        """
        from sage.functions.jacobi import jacobi, jacobi_am
        if diff_param == 0:
            return (sqrt(-m * jacobi('sn', u, m) ** Integer(2) +
                         Integer(1)) * jacobi('dn', u, m))
        elif diff_param == 1:
            return (Integer(1) / Integer(2) *
                    (elliptic_eu(u, m) - elliptic_f(jacobi_am(u, m), m)) / m -
                    Integer(1) / Integer(2) * sqrt(-m * jacobi('sn', u, m) **
                    Integer(2) + Integer(1)) * (m * jacobi('sn', u, m) *
                    jacobi('cn', u, m) - (m - Integer(1)) * u -
                    elliptic_eu(u, m) * jacobi('dn', u, m)) /
                    ((m - Integer(1)) * m))

    def _print_latex_(self, u, m):
        """
        EXAMPLES::

            sage: latex(elliptic_eu(1, x))                                              # needs sage.symbolic
            E(1;x)
        """
        return r"E(%s;%s)" % (latex(u), latex(m))


def elliptic_eu_f(u, m):
    r"""
    Internal function for numeric evaluation of ``elliptic_eu``, defined as
    `E\left(\operatorname{am}(u, m)|m\right)`, where `E` is the incomplete
    elliptic integral of the second kind and `\operatorname{am}` is the Jacobi
    amplitude function.

    EXAMPLES::

        sage: from sage.functions.special import elliptic_eu_f
        sage: elliptic_eu_f(0.5, 0.1)                                                   # needs mpmath
        mpf('0.49605455128659691')
    """
    from mpmath import mp as ctx
    prec = ctx.prec
    try:
        u = ctx.convert(u)
        m = ctx.convert(m)
        ctx.prec += 10
        return ctx.ellipe(jacobi_am_f(u, m), m)
    finally:
        ctx.prec = prec


elliptic_eu = EllipticEU()


class EllipticF(BuiltinFunction):
    r"""
    Return the incomplete elliptic integral of the first kind.

    .. MATH::

        F(\varphi\,|\,m)=\int_0^\varphi \frac{dx}{\sqrt{1 - m\sin(x)^2}},

    Taking `\varphi = \pi/2` gives
    :func:`elliptic_kc()<sage.functions.special.EllipticKC>`.

    EXAMPLES::

        sage: z = var("z")                                                              # needs sage.symbolic
        sage: elliptic_f(z, 0)                                                          # needs sage.symbolic
        z
        sage: elliptic_f(z, 1).simplify()                                               # needs sage.symbolic
        log(tan(1/4*pi + 1/2*z))
        sage: elliptic_f(0.2, 0.1)                                                      # needs mpmath
        0.200132506747543

    .. SEEALSO::

        - :func:`elliptic_e()<sage.functions.special.EllipticE>`.

    REFERENCES:

    - :wikipedia:`Elliptic_integral#Incomplete_elliptic_integral_of_the_first_kind`
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(elliptic_f))
            elliptic_f
            sage: elliptic_f(x, 2)._sympy_()                                            # needs sympy sage.symbolic
            elliptic_f(x, 2)

        Check that :issue:`34186` is fixed::

            sage: _ = var("x y")                                                        # needs sage.symbolic
            sage: fricas(elliptic_f(x, y))                                          # optional - fricas, needs sage.symbolic
            ellipticF(sin(x),y)

        However, the conversion is only correct in the interval
        `[-\pi/2, \pi/2]`::

            sage: fricas(elliptic_f(x, y)).D(x).sage()/elliptic_f(x, y).diff(x)     # optional - fricas, needs sage.symbolic
            cos(x)/sqrt(-sin(x)^2 + 1)

        Numerically::

            sage: f = lambda x, y: elliptic_f(arcsin(x), y).subs(x=x, y=y)
            sage: g = lambda x, y: fricas.ellipticF(x, y).sage()
            sage: d = lambda x, y: f(x, y) - g(x, y)
            sage: [d(N(-pi/2 + x), y)                          # abs tol 1e-8       # optional - fricas, needs sage.symbolic
            ....:  for x in range(1, 3) for y in range(-2,2)]
            [0.000000000000000,
             0.000000000000000,
             0.000000000000000,
             0.000000000000000,
             5.55111512312578e-17,
             0.000000000000000,
             0.000000000000000,
             0.000000000000000]
        """
        BuiltinFunction.__init__(self, 'elliptic_f', nargs=2,
                                 conversions=dict(mathematica='EllipticF',
                                                  maxima='elliptic_f',
                                                  fricas='((x,y)+->ellipticF(sin(x), y))',
                                                  sympy='elliptic_f'))

    def _eval_(self, z, m):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: elliptic_f(x, 1)
            elliptic_f(x, 1)
            sage: elliptic_f(x, 0)
            x
            sage: elliptic_f(0, 1)
            0
            sage: elliptic_f(pi/2, x)
            elliptic_kc(x)
        """
        if m == 0:
            return z
        elif z == 0:
            return Integer(0)
        elif z == pi / 2:
            return elliptic_kc(m)

    def _evalf_(self, z, m, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: elliptic_f(1, 1).n()                                                  # needs sage.symbolic
            1.22619117088352
            sage: elliptic_f(1, 1).n(200)                                               # needs sage.symbolic
            1.22619117088351707081306096...
            sage: elliptic_f(I, I).n()                                                  # needs sage.symbolic
            0.149965060031782 + 0.925097284105771*I
        """
        R = parent or s_parent(z)
        return _mpmath_utils_call(_mpmath_ellipf, z, m, parent=R)

    def _derivative_(self, z, m, diff_param):
        """
        EXAMPLES::

            sage: x, m = var('x,m')                                                     # needs sage.symbolic
            sage: elliptic_f(x, m).diff(x)                                              # needs sage.symbolic
            1/sqrt(-m*sin(x)^2 + 1)
            sage: elliptic_f(x, m).diff(m)                                              # needs sage.symbolic
            -1/2*elliptic_f(x, m)/m
            + 1/4*sin(2*x)/(sqrt(-m*sin(x)^2 + 1)*(m - 1))
            - 1/2*elliptic_e(x, m)/((m - 1)*m)
        """
        if diff_param == 0:
            return Integer(1) / sqrt(Integer(1) - m * sin(z) ** Integer(2))
        elif diff_param == 1:
            return (elliptic_e(z, m) / (Integer(2) * (Integer(1) - m) * m) -
                    elliptic_f(z, m) / (Integer(2) * m) -
                    (sin(Integer(2) * z) /
                     (Integer(4) * (Integer(1) - m) *
                      sqrt(Integer(1) - m * sin(z) ** Integer(2)))))

    def _print_latex_(self, z, m):
        r"""
        EXAMPLES::

            sage: latex(elliptic_f(x, pi))                                              # needs sage.symbolic
            F(x\,|\,\pi)
        """
        return r"F(%s\,|\,%s)" % (latex(z), latex(m))


elliptic_f = EllipticF()


class EllipticKC(BuiltinFunction):
    r"""
    Return the complete elliptic integral of the first kind:

    .. MATH::

        K(m)=\int_0^{\pi/2} \frac{dx}{\sqrt{1 - m\sin(x)^2}}.

    EXAMPLES::

        sage: elliptic_kc(0.5)                                                          # needs mpmath
        1.85407467730137

    .. SEEALSO::

        - :func:`elliptic_f()<sage.functions.special.EllipticF>`.

        - :func:`elliptic_ec()<sage.functions.special.EllipticEC>`.

    REFERENCES:

    - :wikipedia:`Elliptic_integral#Complete_elliptic_integral_of_the_first_kind`

    - :wikipedia:`Elliptic_integral#Incomplete_elliptic_integral_of_the_first_kind`
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: loads(dumps(elliptic_kc))
            elliptic_kc
            sage: elliptic_kc(x)._sympy_()                                              # needs sage.symbolic
            elliptic_k(x)

        TESTS::

            sage: fricas(elliptic_kc(x))                                        # optional - fricas, needs sage.symbolic
            ellipticK(x)

            sage: elliptic_kc(0.3)  # abs tol 1e-8                                      # needs mpmath
            1.71388944817879
            sage: fricas.ellipticK(0.3).sage()  # abs tol 1e-3                  # optional - fricas, needs sage.symbolic
            1.7138894481787910555457043
        """
        BuiltinFunction.__init__(self, 'elliptic_kc', nargs=1, latex_name='K',
                                 conversions=dict(mathematica='EllipticK',
                                                  maxima='elliptic_kc',
                                                  sympy='elliptic_k',
                                                  fricas='ellipticK'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: elliptic_kc(0)                                                        # needs sage.symbolic
            1/2*pi
            sage: elliptic_kc(1/2)                                                      # needs sage.symbolic
            elliptic_kc(1/2)

        TESTS:

        Check if complex numbers in the arguments are converted to maxima
        correctly (see :issue:`7557`)::

            sage: t = jacobi_sn(1.2 + 2*I*elliptic_kc(1 - .5), .5)                      # needs sage.symbolic
            sage: maxima(t)  # abs tol 1e-13                                            # needs sage.symbolic
            0.88771548861928029 - 1.7301614091485560e-15*%i
            sage: t.n()  # abs tol 1e-13                                                # needs sage.symbolic
            0.887715488619280 - 1.73016140914856e-15*I
        """
        if z == 0:
            return pi / 2
        else:
            return None

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: elliptic_kc(1/2).n()                                                  # needs sage.symbolic
            1.85407467730137
            sage: elliptic_kc(1/2).n(200)                                               # needs sage.symbolic
            1.85407467730137191843385034...
            sage: elliptic_kc(I).n()                                                    # needs sage.symbolic
            1.42127228104504 + 0.295380284214777*I
        """
        R = parent or s_parent(z)
        return _mpmath_utils_call(_mpmath_ellipk, z, parent=R)

    def _derivative_(self, z, diff_param):
        """
        EXAMPLES::

            sage: elliptic_kc(x).diff(x)                                                # needs sage.symbolic
            -1/2*((x - 1)*elliptic_kc(x)
            + elliptic_ec(x))/((x - 1)*x)
        """
        return ((elliptic_ec(z) - (Integer(1) - z) * elliptic_kc(z)) /
                (Integer(2) * (Integer(1) - z) * z))


elliptic_kc = EllipticKC()


class EllipticPi(BuiltinFunction):
    r"""
    Return the incomplete elliptic integral of the third kind:

    .. MATH::

        \Pi(n, t, m) = \int_0^t \frac{dx}{(1 - n \sin(x)^2)\sqrt{1 - m \sin(x)^2}}.

    INPUT:

    - ``n`` -- a real number, called the "characteristic"

    - ``t`` -- a real number, called the "amplitude"

    - ``m`` -- a real number, called the "parameter"

    EXAMPLES::

        sage: N(elliptic_pi(1, pi/4, 1))                                                # needs sage.symbolic
        1.14779357469632

    Compare the value computed by Maxima to the definition as a definite integral
    (using GSL)::

        sage: elliptic_pi(0.1, 0.2, 0.3)                                                # needs mpmath
        0.200665068220979
        sage: numerical_integral(1/(1-0.1*sin(x)^2)/sqrt(1-0.3*sin(x)^2), 0.0, 0.2)     # needs sage.symbolic
        (0.2006650682209791, 2.227829789769088e-15)

    REFERENCES:

    - :wikipedia:`Elliptic_integral#Incomplete_elliptic_integral_of_the_third_kind`
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: loads(dumps(elliptic_pi))
            elliptic_pi
            sage: elliptic_pi(x, pi/4, 1)._sympy_()                                     # needs sympy sage.symbolic
            elliptic_pi(x, pi/4, 1)
        """
        BuiltinFunction.__init__(self, 'elliptic_pi', nargs=3,
                                 conversions=dict(mathematica='EllipticPi',
                                                  maxima='EllipticPi',
                                                  # fricas='ellipticPi', doubt
                                                  sympy='elliptic_pi'))

    def _eval_(self, n, z, m):
        """
        EXAMPLES::

            sage: elliptic_pi(x, x, pi)                                                 # needs sympy sage.symbolic
            elliptic_pi(x, x, pi)
            sage: elliptic_pi(0, x, pi)                                                 # needs sympy sage.symbolic
            elliptic_f(x, pi)
        """
        if n == 0:
            return elliptic_f(z, m)

    def _evalf_(self, n, z, m, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: elliptic_pi(pi, 1/2, 1).n()
            0.795062820631931
            sage: elliptic_pi(pi, 1/2, 1).n(200)
            0.79506282063193125292514098445...
            sage: elliptic_pi(pi, 1, 1).n()
            0.0991592574231369 - 1.30004368185937*I
            sage: elliptic_pi(pi, I, I).n()
            0.0542471560940594 + 0.552096453413081*I
        """
        R = parent or s_parent(z)
        return _mpmath_utils_call(_mpmath_ellippi, n, z, m, parent=R)

    def _derivative_(self, n, z, m, diff_param):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: n, z, m = var('n,z,m')
            sage: elliptic_pi(n, z, m).diff(n)
            1/4*(sqrt(-m*sin(z)^2 + 1)*n*sin(2*z)/(n*sin(z)^2 - 1)
            + 2*(m - n)*elliptic_f(z, m)/n
            + 2*(n^2 - m)*elliptic_pi(n, z, m)/n
            + 2*elliptic_e(z, m))/((m - n)*(n - 1))
            sage: elliptic_pi(n, z, m).diff(z)
            -1/(sqrt(-m*sin(z)^2 + 1)*(n*sin(z)^2 - 1))
            sage: elliptic_pi(n, z, m).diff(m)
            1/4*(m*sin(2*z)/(sqrt(-m*sin(z)^2 + 1)*(m - 1))
            - 2*elliptic_e(z, m)/(m - 1)
            - 2*elliptic_pi(n, z, m))/(m - n)
        """
        if diff_param == 0:
            return ((Integer(1) / (Integer(2) * (m - n) * (n - Integer(1)))) *
                    (elliptic_e(z, m) + ((m - n) / n) * elliptic_f(z, m) +
                    ((n ** Integer(2) - m) / n) * elliptic_pi(n, z, m) -
                    (n * sqrt(Integer(1) - m * sin(z) ** Integer(2)) *
                     sin(Integer(2) * z)) /
                    (Integer(2) * (Integer(1) - n * sin(z) ** Integer(2)))))
        elif diff_param == 1:
            return (Integer(1) /
                    (sqrt(Integer(1) - m * sin(z) ** Integer(Integer(2))) *
                     (Integer(1) - n * sin(z) ** Integer(2))))
        elif diff_param == 2:
            return ((Integer(1) / (Integer(2) * (n - m))) *
                    (elliptic_e(z, m) / (m - Integer(1)) +
                     elliptic_pi(n, z, m) - (m * sin(Integer(2) * z)) /
                     (Integer(2) * (m - Integer(1)) *
                     sqrt(Integer(1) - m * sin(z) ** Integer(2)))))

    def _print_latex_(self, n, z, m):
        r"""
        EXAMPLES::

            sage: latex(elliptic_pi(x, pi, 0))                                          # needs sage.symbolic
            \Pi(x,\pi,0)
        """
        return r"\Pi(%s,%s,%s)" % (latex(n), latex(z), latex(m))


elliptic_pi = EllipticPi()
