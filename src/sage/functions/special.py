r"""
Special Functions

AUTHORS:

- David Joyner (2006-13-06): initial version

- David Joyner (2006-30-10): bug fixes to pari wrappers of Bessel
  functions, hypergeometric_U

- William Stein (2008-02): Impose some sanity checks.

- David Joyner (2008-04-23): addition of elliptic integrals

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


-  Spherical harmonics: Laplace's equation in spherical coordinates
   is:

   .. math::

       {\frac{1}{r^2}}{\frac{\partial}{\partial r}}   \left(r^2 {\frac{\partial f}{\partial r}}\right) +   {\frac{1}{r^2}\sin\theta}{\frac{\partial}{\partial \theta}}   \left(\sin\theta {\frac{\partial f}{\partial \theta}}\right) +   {\frac{1}{r^2\sin^2\theta}}{\frac{\partial^2 f}{\partial \varphi^2}} = 0.


   Note that the spherical coordinates `\theta` and
   `\varphi` are defined here as follows: `\theta` is
   the colatitude or polar angle, ranging from
   `0\leq\theta\leq\pi` and `\varphi` the azimuth or
   longitude, ranging from `0\leq\varphi<2\pi`.

   The general solution which remains finite towards infinity is a
   linear combination of functions of the form

   .. math::

         r^{-1-\ell} \cos (m \varphi) P_\ell^m (\cos{\theta} )


   and

   .. math::

         r^{-1-\ell} \sin (m \varphi) P_\ell^m (\cos{\theta} )


   where `P_\ell^m` are the associated Legendre polynomials,
   and with integer parameters `\ell \ge 0` and `m`
   from `0` to `\ell`. Put in another way, the
   solutions with integer parameters `\ell \ge 0` and
   `- \ell\leq m\leq \ell`, can be written as linear
   combinations of:

   .. math::

         U_{\ell,m}(r,\theta , \varphi ) = r^{-1-\ell} Y_\ell^m( \theta , \varphi )


   where the functions `Y` are the spherical harmonic
   functions with parameters `\ell`, `m`, which can be
   written as:

   .. math::

         Y_\ell^m( \theta , \varphi )     = \sqrt{{\frac{(2\ell+1)}{4\pi}}{\frac{(\ell-m)!}{(\ell+m)!}}}       \cdot e^{i m \varphi } \cdot P_\ell^m ( \cos{\theta} ) .



   The spherical harmonics obey the normalisation condition


   .. math::

     \int_{\theta=0}^\pi\int_{\varphi=0}^{2\pi} Y_\ell^mY_{\ell'}^{m'*}\,d\Omega =\delta_{\ell\ell'}\delta_{mm'}\quad\quad d\Omega =\sin\theta\,d\varphi\,d\theta .



-  When solving for separable solutions of Laplace's equation in
   spherical coordinates, the radial equation has the form:

   .. math::

         x^2 \frac{d^2 y}{dx^2} + 2x \frac{dy}{dx} + [x^2 - n(n+1)]y = 0.


   The spherical Bessel functions `j_n` and `y_n`,
   are two linearly independent solutions to this equation. They are
   related to the ordinary Bessel functions `J_n` and
   `Y_n` by:

   .. math::

         j_n(x) = \sqrt{\frac{\pi}{2x}} J_{n+1/2}(x),



   .. math::

         y_n(x) = \sqrt{\frac{\pi}{2x}} Y_{n+1/2}(x)     = (-1)^{n+1} \sqrt{\frac{\pi}{2x}} J_{-n-1/2}(x).



-  For `x>0`, the confluent hypergeometric function
   `y = U(a,b,x)` is defined to be the solution to Kummer's
   differential equation


   .. math::

     xy'' + (b-x)y' - ay = 0,

   which satisfies `U(a,b,x) \sim x^{-a}`, as
   `x\rightarrow \infty`. (There is a linearly independent
   solution, called Kummer's function `M(a,b,x)`, which is not
   implemented.)

   -  The incomplete elliptic integrals (of the first kind, etc.) are:

      .. math::

         \begin{array}{c} \displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2}}\, dx,\\ \displaystyle\int_0^\phi \sqrt{1 - m\sin(x)^2}\, dx,\\ \displaystyle\int_0^\phi \frac{\sqrt{1-mt^2}}{\sqrt(1 - t^2)}\, dx,\\ \displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2\sqrt{1 - n\sin(x)^2}}}\, dx, \end{array}

      and the complete ones are obtained by taking `\phi =\pi/2`.


REFERENCES:

- Abramowitz and Stegun: Handbook of Mathematical Functions,
  http://www.math.sfu.ca/~cbm/aands/

- http://en.wikipedia.org/wiki/Spherical_harmonics

- http://en.wikipedia.org/wiki/Helmholtz_equation

- Online Encyclopedia of Special Function
  http://algo.inria.fr/esf/index.html

TODO: Resolve weird bug in commented out code in hypergeometric_U
below.

AUTHORS:

- David Joyner and William Stein

Added 16-02-2008 (wdj): optional calls to scipy and replace all
'#random' by '...' (both at the request of William Stein)

.. warning::

   SciPy's versions are poorly documented and seem less
   accurate than the Maxima and PARI versions; typically they are limited
   by hardware floats precision.
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.real_mpfr import RealField
from sage.rings.complex_field import ComplexField
from sage.misc.latex import latex
from sage.rings.all import ZZ, RR, RDF, CDF
from sage.structure.parent import Parent
from sage.functions.other import log_gamma
from sage.symbolic.function import BuiltinFunction
from sage.calculus.calculus import maxima
from sage.libs.mpmath import utils as mpmath_utils
from sage.functions.all import sqrt, cot, exp
from sage.symbolic.all import I

_done = False
def _init():
    """
    Internal function which checks if Maxima has loaded the
    "orthopoly" package.  All functions using this in this
    file should call this function first.

    TEST:

    The global starts ``False``::

        sage: sage.functions.special._done
        False

    Then after using one of the MaximaFunctions, it changes::

        sage: from sage.functions.special import elliptic_ec
        sage: elliptic_ec(0.1)
        1.53075763689776

        sage: sage.functions.special._done
        True
    """
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True

def meval(x):
    """
    Return ``x`` evaluated in Maxima, then returned to Sage.

    This is used to evaluate several of these special functions.

    TEST::

        sage: from sage.functions.special import spherical_bessel_J
        sage: spherical_bessel_J(2.,3.)      # rel tol 1e-10
        0.2986374970757335
    """
    return maxima(x).sage()


class MaximaFunction(BuiltinFunction):
    """
    EXAMPLES::

        sage: from sage.functions.special import MaximaFunction
        sage: f = MaximaFunction("jacobi_sn")
        sage: f(1,1)
        tanh(1)
        sage: f(1/2,1/2).n()
        0.470750473655657
    """
    def __init__(self, name, nargs=2, conversions={}):
        """
        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f(1,1)
            tanh(1)
            sage: f(1/2,1/2).n()
            0.470750473655657
        """
        c = dict(maxima=name)
        c.update(conversions)
        BuiltinFunction.__init__(self, name=name, nargs=nargs,
                                   conversions=c)

    def _maxima_init_evaled_(self, *args):
        """
        Returns a string which represents this function evaluated at
        *args* in Maxima.

        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f._maxima_init_evaled_(1/2, 1/2)
            'jacobi_sn(1/2, 1/2)'

        TESTS:

        Check if complex numbers in the arguments are converted to maxima
        correctly (see :trac:`7557`)::

            sage: t = f(1.2+2*I*elliptic_kc(1-.5),.5)
            sage: maxima(t) # abs tol 1e-13
            0.88771548861928029 - 1.7301614091485560e-15*%i
            sage: t.n() # abs tol 1e-13
            0.887715488619280 - 1.73016140914856e-15*I
        """
        args_maxima = []
        for a in args:
            if isinstance(a, str):
                args_maxima.append(a)
            elif hasattr(a, '_maxima_init_'):
                args_maxima.append(a._maxima_init_())
            else:
                args_maxima.append(str(a))
        return "%s(%s)"%(self.name(), ', '.join(args_maxima))

    def _evalf_(self, *args, **kwds):
        """
        Returns a numerical approximation of this function using
        Maxima.  Currently, this is limited to 53 bits of precision.

        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f(1/2, 1/2)
            jacobi_sn(1/2, 1/2)
            sage: f(1/2, 1/2).n()
            0.470750473655657
            sage: f(1/2, 1/2).n(20)
            0.47075
            sage: f(1, I).n()
            0.848379519751901 - 0.0742924572771414*I

        TESTS::

            sage: f(1/2, 1/2).n(150)
            Traceback (most recent call last):
            ...
            NotImplementedError: Maxima function jacobi_sn not implemented for Real Field with 150 bits of precision
            sage: f._evalf_(1/2, 1/2, parent=int)
            Traceback (most recent call last):
            ...
            NotImplementedError: Maxima function jacobi_sn not implemented for <type 'int'>
            sage: f._evalf_(1/2, 1/2, parent=complex)
            (0.4707504736556572+0j)
            sage: f._evalf_(1/2, 1/2, parent=RDF)
            0.4707504736556572
            sage: f._evalf_(1, I, parent=CDF)  # abs tol 1e-16
            0.8483795707591759 - 0.07429247342160791*I
            sage: f._evalf_(1, I, parent=RR)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='0.848379570759176-0.0742924734216079*I') to real number.
        """
        parent = kwds['parent']
        # The result from maxima is a machine double, which corresponds
        # to RDF (or CDF). Therefore, before converting, we check that
        # we can actually coerce RDF into our parent.
        if parent is not float and parent is not complex:
            if not isinstance(parent, Parent) or not parent.has_coerce_map_from(RDF):
                raise NotImplementedError("Maxima function %s not implemented for %r"%(self.name(), parent))
        _init()
        return parent(maxima("%s, numer"%self._maxima_init_evaled_(*args)))

    def _eval_(self, *args):
        """
        Try to evaluate this function at ``*args``, return ``None`` if
        Maxima did not compute a numerical evaluation.

        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f(1,1)
            tanh(1)

            sage: f._eval_(1,1)
            tanh(1)

        Here arccoth doesn't have 1 in its domain, so we just hold the expression:

            sage: elliptic_e(arccoth(1), x^2*e)
            elliptic_e(arccoth(1), x^2*e)

        Since Maxima works only with double precision, numerical
        results are in ``RDF``, no matter what the input precision is::

            sage: R = RealField(300)
            sage: r = elliptic_eu(R(1/2), R(1/8)); r
            0.4950737320232015
            sage: parent(r)
            Real Double Field
        """
        _init()
        try:
            s = maxima(self._maxima_init_evaled_(*args))
        except TypeError:
            return None

        if self.name() in s.__repr__():  # Avoid infinite recursion
            return None
        else:
            return s.sage()

from sage.misc.cachefunc import cached_function

@cached_function
def maxima_function(name):
    """
    Returns a function which is evaluated both symbolically and
    numerically via Maxima.  In particular, it returns an instance
    of :class:`MaximaFunction`.

    .. note::

       This function is cached so that duplicate copies of the same
       function are not created.

    EXAMPLES::

        sage: spherical_hankel2(2,i)
        -e
    """
    # The superclass of MaximaFunction, BuiltinFunction, assumes that there
    # will be only one symbolic function with the same name and class.
    # We create a new class for each Maxima function wrapped.
    class NewMaximaFunction(MaximaFunction):
        def __init__(self):
            """
            Constructs an object that wraps a Maxima function.

            TESTS::

                sage: spherical_hankel2(2,x)
                (-I*x^2 - 3*x + 3*I)*e^(-I*x)/x^3
            """
            MaximaFunction.__init__(self, name)

    return NewMaximaFunction()


def hypergeometric_U(alpha,beta,x,algorithm="pari",prec=53):
    r"""
    Default is a wrap of PARI's hyperu(alpha,beta,x) function.
    Optionally, algorithm = "scipy" can be used.

    The confluent hypergeometric function `y = U(a,b,x)` is
    defined to be the solution to Kummer's differential equation

    .. math::

             xy'' + (b-x)y' - ay = 0.

    This satisfies `U(a,b,x) \sim x^{-a}`, as
    `x\rightarrow \infty`, and is sometimes denoted
    ``x^{-a}2_F_0(a,1+a-b,-1/x)``. This is not the same as Kummer's
    `M`-hypergeometric function, denoted sometimes as
    ``_1F_1(alpha,beta,x)``, though it satisfies the same DE that
    `U` does.

    .. warning::

       In the literature, both are called "Kummer confluent
       hypergeometric" functions.

    EXAMPLES::

        sage: hypergeometric_U(1,1,1,"scipy")
        0.596347362323...
        sage: hypergeometric_U(1,1,1)
        0.59634736232319...
        sage: hypergeometric_U(1,1,1,"pari",70)
        0.59634736232319407434...
    """
    if algorithm == "scipy":
        if prec != 53:
            raise ValueError("for the scipy algorithm the precision must be 53")
        import scipy.special
        return RDF(scipy.special.hyperu(float(alpha), float(beta), float(x)))
    elif algorithm == 'pari':
        from sage.libs.pari.all import pari
        R = RealField(prec)
        return R(pari(R(alpha)).hyperu(R(beta), R(x), precision=prec))
    else:
        raise ValueError("unknown algorithm '%s'" % algorithm)

def spherical_bessel_J(n, var, algorithm="maxima"):
    r"""
    Returns the spherical Bessel function of the first kind for
    integers n >= 1.

    Reference: AS 10.1.8 page 437 and AS 10.1.15 page 439.

    EXAMPLES::

        sage: spherical_bessel_J(2,x)
        ((3/x^2 - 1)*sin(x) - 3*cos(x)/x)/x
        sage: spherical_bessel_J(1, 5.2, algorithm='scipy')
        -0.12277149950007...
        sage: spherical_bessel_J(1, 3, algorithm='scipy')
        0.345677499762355...
    """
    if algorithm == "scipy":
        from scipy.special.specfun import sphj
        return CDF(sphj(int(n), float(var))[1][-1])
    elif algorithm == 'maxima':
        _init()
        return meval("spherical_bessel_j(%s,%s)"%(ZZ(n),var))
    else:
        raise ValueError("unknown algorithm '%s'"%algorithm)

def spherical_bessel_Y(n,var, algorithm="maxima"):
    r"""
    Returns the spherical Bessel function of the second kind for
    integers n -1.

    Reference: AS 10.1.9 page 437 and AS 10.1.15 page 439.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: spherical_bessel_Y(2,x)
        -((3/x^2 - 1)*cos(x) + 3*sin(x)/x)/x
    """
    if algorithm == "scipy":
        import scipy.special
        return CDF(scipy.special.sph_yn(int(n),float(var)))
    elif algorithm == 'maxima':
        _init()
        return meval("spherical_bessel_y(%s,%s)"%(ZZ(n),var))
    else:
        raise ValueError("unknown algorithm '%s'"%algorithm)

def spherical_hankel1(n, var):
    r"""
    Returns the spherical Hankel function of the first kind for
    integers `n > -1`, written as a string. Reference: AS
    10.1.36 page 439.

    EXAMPLES::

        sage: spherical_hankel1(2, x)
        (I*x^2 - 3*x - 3*I)*e^(I*x)/x^3
    """
    return maxima_function("spherical_hankel1")(ZZ(n), var)

def spherical_hankel2(n,x):
    r"""
    Returns the spherical Hankel function of the second kind for
    integers `n > -1`, written as a string. Reference: AS 10.1.17 page
    439.

    EXAMPLES::

        sage: spherical_hankel2(2, x)
        (-I*x^2 - 3*x + 3*I)*e^(-I*x)/x^3

    Here I = sqrt(-1).
    """
    return maxima_function("spherical_hankel2")(ZZ(n), x)


class SphericalHarmonic(BuiltinFunction):
    r"""
    Returns the spherical harmonic function `Y_n^m(\theta, \varphi)`.

    For integers `n > -1`, `|m| \leq n`, simplification is done automatically.
    Numeric evaluation is supported for complex `n` and `m`.

    EXAMPLES::

        sage: x, y = var('x, y')
        sage: spherical_harmonic(3, 2, x, y)
        15/4*sqrt(7/30)*cos(x)*e^(2*I*y)*sin(x)^2/sqrt(pi)
        sage: spherical_harmonic(3, 2, 1, 2)
        15/4*sqrt(7/30)*cos(1)*e^(4*I)*sin(1)^2/sqrt(pi)
        sage: spherical_harmonic(3 + I, 2., 1, 2)
        -0.351154337307488 - 0.415562233975369*I
        sage: latex(spherical_harmonic(3, 2, x, y, hold=True))
        Y_{3}^{2}\left(x, y\right)
        sage: spherical_harmonic(1, 2, x, y)
        0
    """
    def __init__(self):
        r"""
        TESTS::

            sage: n, m, theta, phi = var('n m theta phi')
            sage: spherical_harmonic(n, m, theta, phi)._sympy_()
            Ynm(n, m, theta, phi)
        """
        BuiltinFunction.__init__(self, 'spherical_harmonic', nargs=4,
                                 conversions=dict(
                                    maple='SphericalY',
                                    mathematica= 'SphericalHarmonicY',
                                    maxima='spherical_harmonic',
                                    sympy='Ynm'))

    def _eval_(self, n, m, theta, phi, **kwargs):
        r"""
        TESTS::

            sage: x, y = var('x y')
            sage: spherical_harmonic(1, 2, x, y)
            0
            sage: spherical_harmonic(1, -2, x, y)
            0
            sage: spherical_harmonic(1/2, 2, x, y)
            spherical_harmonic(1/2, 2, x, y)
            sage: spherical_harmonic(3, 2, x, y)
            15/4*sqrt(7/30)*cos(x)*e^(2*I*y)*sin(x)^2/sqrt(pi)
            sage: spherical_harmonic(3, 2, 1, 2)
            15/4*sqrt(7/30)*cos(1)*e^(4*I)*sin(1)^2/sqrt(pi)
            sage: spherical_harmonic(3 + I, 2., 1, 2)
            -0.351154337307488 - 0.415562233975369*I
        """
        if n in ZZ and m in ZZ and n > -1:
            if abs(m) > n:
                return ZZ(0)
            return meval("spherical_harmonic({},{},{},{})".format(
                ZZ(n), ZZ(m), maxima(theta), maxima(phi)))

    def _evalf_(self, n, m, theta, phi, parent, **kwds):
        r"""
        TESTS::

            sage: spherical_harmonic(3 + I, 2, 1, 2).n(100)
            -0.35115433730748836508201061672 - 0.41556223397536866209990358597*I
            sage: spherical_harmonic(I, I, I, I).n()
            7.66678546069894 - 0.265754432549751*I
        """
        from mpmath import spherharm
        return mpmath_utils.call(spherharm, n, m, theta, phi, parent=parent)

    def _derivative_(self, n, m, theta, phi, diff_param):
        r"""
        TESTS::

            sage: n, m, theta, phi = var('n m theta phi')
            sage: spherical_harmonic(n, m, theta, phi).diff(theta)
            m*cot(theta)*spherical_harmonic(n, m, theta, phi)
             + sqrt(-(m + n + 1)*(m - n))*e^(-I*phi)*spherical_harmonic(n, m + 1, theta, phi)
            sage: spherical_harmonic(n, m, theta, phi).diff(phi)
            I*m*spherical_harmonic(n, m, theta, phi)
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

            sage: y = var('y')
            sage: latex(spherical_harmonic(3, 2, x, y, hold=True))
            Y_{3}^{2}\left(x, y\right)
        """
        return r"Y_{{{}}}^{{{}}}\left({}, {}\right)".format(
                 latex(n), latex(m), latex(theta), latex(phi))

spherical_harmonic = SphericalHarmonic()

####### elliptic functions and integrals

def elliptic_j(z):
   r"""
   Returns the elliptic modular `j`-function evaluated at `z`.

   INPUT:

   - ``z`` (complex) -- a complex number with positive imaginary part.

   OUTPUT:

   (complex) The value of `j(z)`.

   ALGORITHM:

   Calls the ``pari`` function ``ellj()``.

   AUTHOR:

   John Cremona

   EXAMPLES::

       sage: elliptic_j(CC(i))
       1728.00000000000
       sage: elliptic_j(sqrt(-2.0))
       8000.00000000000
       sage: z = ComplexField(100)(1,sqrt(11))/2
       sage: elliptic_j(z)
       -32768.000...
       sage: elliptic_j(z).real().round()
       -32768

    ::

       sage: tau = (1 + sqrt(-163))/2
       sage: (-elliptic_j(tau.n(100)).real().round())^(1/3)
       640320

   """
   CC = z.parent()
   from sage.rings.complex_field import is_ComplexField
   if not is_ComplexField(CC):
      CC = ComplexField()
      try:
         z = CC(z)
      except ValueError:
         raise ValueError("elliptic_j only defined for complex arguments.")
   from sage.libs.all import pari
   return CC(pari(z).ellj())

#### elliptic integrals

class EllipticE(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    second kind,"

    .. math::

        \int_0^\phi \sqrt{1 - m\sin(x)^2}\, dx,

    i.e., ``integrate(sqrt(1 - m*sin(x)^2), x, 0, phi)``.  Taking `\phi
    = \pi/2` gives ``elliptic_ec``.

    EXAMPLES::

        sage: z = var("z")
        sage: # this is still wrong: must be abs(sin(z)) + 2*round(z/pi)
        sage: elliptic_e(z, 1)
        2*round(z/pi) + sin(z)
        sage: elliptic_e(z, 0)
        z
        sage: elliptic_e(0.5, 0.1)  # abs tol 2e-15
        0.498011394498832
    """
    def __init__(self):
        """
        TESTS::

            sage: loads(dumps(elliptic_e))
            elliptic_e
        """
        MaximaFunction.__init__(self, "elliptic_e")

elliptic_e = EllipticE()

class EllipticEC(MaximaFunction):
    """
    This returns the value of the "complete elliptic integral of the
    second kind,"

    .. math::

        \int_0^{\pi/2} \sqrt{1 - m\sin(x)^2}\, dx.

    EXAMPLES::

        sage: elliptic_ec(0.1)
        1.53075763689776
        sage: elliptic_ec(x).diff()
        1/2*(elliptic_ec(x) - elliptic_kc(x))/x

        sage: loads(dumps(elliptic_ec))
        elliptic_ec
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: elliptic_ec(0.1)
            1.53075763689776

        TESTS::

            sage: elliptic_ec(x)._sympy_()
            elliptic_e(x)
        """
        MaximaFunction.__init__(self, "elliptic_ec", nargs=1,
                                conversions=dict(sympy='elliptic_e'))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: elliptic_ec(x).diff()
            1/2*(elliptic_ec(x) - elliptic_kc(x))/x
        """
        diff_param = kwds['diff_param']
        assert diff_param == 0
        x = args[diff_param]
        return (elliptic_ec(x) - elliptic_kc(x))/(2*x)

elliptic_ec = EllipticEC()


class EllipticEU(MaximaFunction):
    r"""
    Return the value of the "incomplete elliptic integral of the
    second kind,"

    .. math::

        \int_0^u \mathrm{dn}(x,m)^2\, dx = \int_0^\tau
        {\sqrt{1-m x^2}\over\sqrt{1-x^2}}\, dx.

    where `\tau = \mathrm{sn}(u, m)`.

    EXAMPLES::

        sage: elliptic_eu (0.5, 0.1)
        0.496054551286597
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_eu (0.5, 0.1)
            0.496054551286597
        """
        MaximaFunction.__init__(self, "elliptic_eu")

elliptic_eu = EllipticEU()

class EllipticF(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    first kind,"

    .. math::

        \int_0^\phi \frac{dx}{\sqrt{1 - m\sin(x)^2}},

    i.e., ``integrate(1/sqrt(1 - m*sin(x)^2), x, 0, phi)``.  Taking
    `\phi = \pi/2` gives ``elliptic_kc``.

    EXAMPLES::

        sage: z = var("z")
        sage: elliptic_f (z, 0)
        z
        sage: elliptic_f (z, 1)
        log(tan(1/4*pi + 1/2*z))
        sage: elliptic_f (0.2, 0.1)
        0.200132506747543
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_f (0.2, 0.1)
            0.200132506747543

        TESTS::

            sage: elliptic_f(x, 2)._sympy_()
            elliptic_f(x, 2)
        """
        MaximaFunction.__init__(self, "elliptic_f",
                                conversions=dict(sympy='elliptic_f'))

elliptic_f = EllipticF()

class EllipticKC(MaximaFunction):
    r"""
    This returns the value of the "complete elliptic integral of the
    first kind,"

    .. math::

        \int_0^{\pi/2} \frac{dx}{\sqrt{1 - m\sin(x)^2}}.

    EXAMPLES::

        sage: elliptic_kc(0.5)
        1.85407467730137
        sage: elliptic_f(RR(pi/2), 0.5)
        1.85407467730137
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_kc(0.5)
            1.85407467730137
            sage: elliptic_f(RR(pi/2), 0.5)
            1.85407467730137

        TESTS::

            sage: elliptic_kc(x)._sympy_()
            elliptic_k(x)
        """
        MaximaFunction.__init__(self, "elliptic_kc", nargs=1,
                                conversions=dict(sympy='elliptic_k'))

elliptic_kc = EllipticKC()

class EllipticPi(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    third kind,"

    .. math::

        \text{elliptic\_pi}(n, t, m) = \int_0^t \frac{dx}{(1 - n \sin(x)^2)
        \sqrt{1 - m \sin(x)^2}}.

    INPUT:

    - ``n`` -- a real number, called the "characteristic"

    - ``t`` -- a real number, called the "amplitude"

    - ``m`` -- a real number, called the "parameter"

    EXAMPLES::

        sage: N(elliptic_pi(1, pi/4, 1))
        1.14779357469632

    Compare the value computed by Maxima to the definition as a definite integral
    (using GSL)::

        sage: elliptic_pi(0.1, 0.2, 0.3)
        0.200665068220979
        sage: numerical_integral(1/(1-0.1*sin(x)^2)/sqrt(1-0.3*sin(x)^2), 0.0, 0.2)
        (0.2006650682209791, 2.227829789769088e-15)

    ALGORITHM:

    Numerical evaluation and symbolic manipulation are provided by `Maxima`_.

    REFERENCES:

    - Abramowitz and Stegun: Handbook of Mathematical Functions, section 17.7
      http://www.math.sfu.ca/~cbm/aands/
    - Elliptic Functions in `Maxima`_

    .. _`Maxima`: http://maxima.sourceforge.net/docs/manual/en/maxima_16.html#SEC91
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_pi(0.1, 0.2, 0.3)
            0.200665068220979

        TESTS::

            sage: elliptic_pi(x, pi/4, 1)._sympy_()
            elliptic_pi(x, pi/4, 1)
        """
        MaximaFunction.__init__(self, "elliptic_pi", nargs=3,
                                conversions=dict(sympy='elliptic_pi'))

elliptic_pi = EllipticPi()

def error_fcn(t):
    r"""
    The complementary error function
    `\frac{2}{\sqrt{\pi}}\int_t^\infty e^{-x^2} dx` (t belongs
    to RR).  This function is currently always
    evaluated immediately.

    EXAMPLES::

        sage: error_fcn(6)
        2.15197367124989e-17
        sage: error_fcn(RealField(100)(1/2))
        0.47950012218695346231725334611

    Note this is literally equal to `1 - erf(t)`::

        sage: 1 - error_fcn(0.5)
        0.520499877813047
        sage: erf(0.5)
        0.520499877813047
    """
    try:
        return t.erfc()
    except AttributeError:
        try:
            return RR(t).erfc()
        except Exception:
            raise NotImplementedError



