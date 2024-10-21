# sage_setup: distribution = sagemath-repl
# sage.doctest: needs sage.symbolic
"""
Benchmarks

COMMENTS:

Taken as a whole these benchmarks suggest that by far the fastest math
software is MAGMA, Mathematica, and Sage (with appropriate tuning and
choices -- it can also be very slow !).  Maxima is very slow, at least
in the form that comes with Sage (perhaps this is because of using
clisp, at least partly).  GP is slow at some thing and very fast at
others.

TESTS::

    sage: import sage.tests.benchmark
"""
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
from sage.combinat.combinat import fibonacci
from sage.functions.other import factorial
from sage.interfaces.gp import gp
from sage.interfaces.macaulay2 import macaulay2
from sage.interfaces.magma import magma, Magma
from sage.interfaces.maple import maple
from sage.interfaces.mathematica import mathematica
from sage.interfaces.maxima import maxima
from sage.interfaces.singular import singular
from sage.libs.pari import pari
from sage.libs.gap.libgap import libgap
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.functional import log
from sage.misc.timing import cputime, walltime
from sage.modular.modsym.modsym import ModularSymbols
from sage.rings.complex_mpfr import ComplexField
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.constructor import EllipticCurve


def avg(X):
    """
    Return the average of the list X.

    EXAMPLES::

        sage: from sage.tests.benchmark import avg
        sage: avg([1,2,3])
        2.0
    """
    s = sum(X, 0)
    return s / float(len(X))


STD_SYSTEMS = ['sage', 'maxima', 'libgap', 'gp', 'pari', 'python']
OPT_SYSTEMS = ['magma', 'macaulay2', 'maple', 'mathematica']


class Benchmark:
    """
    A class for running specific benchmarks against different systems.

    In order to implement an extension of this class, one must write
    functions named after the different systems one wants to test. These
    functions must perform the same task for each function. Calling
    the run command with a list of systems will then show the timings.

    EXAMPLES::

        sage: from sage.tests.benchmark import Benchmark
        sage: B = Benchmark()
        sage: def python():
        ....:     t = cputime()
        ....:     2+2
        ....:     return cputime(t)
        sage: B.python = python
        sage: B.run(systems=['python'])
        sage.tests.benchmark.Benchmark instance
          System      min         avg         max         trials          cpu or wall
        * python ...
    """
    def run(self, systems=None, timeout=60, trials=1, sort=False, optional=False):
        """
        Run the benchmarking functions for the current benchmark on the systems
        given.

        INPUT:

        - ``systems`` -- (optional) list of strings of which systems to run tests on;
          if ``None``, runs the standard systems
        - ``timeout`` -- integer (default: 60); how long (in seconds)
          to run each test for
        - ``trials`` -- integer (default: 1); number of trials
        - ``sort`` -- boolean (default: ``False``); whether to sort system names
        - ``optional`` -- boolean (default: ``False``); if systems is ``None``,
          whether to test optional systems

        EXAMPLES::

            sage: from sage.tests.benchmark import PolyFactor
            sage: PolyFactor(10, QQ).run()
            Factor a product of 2 polynomials of degree 10 over Rational Field.
              System      min         avg         max         trials          cpu or wall
            * sage        ...
            * gp          ...
        """
        if sort:
            systems.sort()
        print('\n\n\n' + str(self))
        print('  %-12s%-12s%-12s%-12s%-12s%15s' % ('System', 'min',
                                                   'avg', 'max',
                                                   'trials', 'cpu or wall'))
        if systems is None:
            systems = STD_SYSTEMS
            if optional:
                systems += OPT_SYSTEMS
        for S in systems:
            try:
                X = []
                wall = False
                for i in range(trials):
                    alarm(timeout)
                    t = getattr(self, S)()
                    cancel_alarm()
                    if isinstance(t, tuple):
                        wall = True
                        t = t[1]
                    X.append(t)
                mn = min(X)
                mx = max(X)
                av = avg(X)
                s = '* %-12s%-12f%-12f%-12f%-12s' % (S, mn, av,
                                                     mx, trials)
                if wall:
                    s += '%15fw' % t
                else:
                    s += '%15fc' % t
                print(s)
            except AlarmInterrupt:
                print('%-12sinterrupted (timeout: %s seconds wall time)' %
                      (S, timeout))
            except AttributeError:
                pass
            except Exception as msg:
                print(msg)

    bench = run

    def __repr__(self):
        """
        Print representation of self, simply coming from self.repr_str.

        EXAMPLES::

            sage: from sage.tests.benchmark import Benchmark
            sage: B = Benchmark()
            sage: B.repr_str = 'spam'
            sage: B
            spam
        """
        try:
            return self.repr_str
        except AttributeError:
            return 'sage.tests.benchmark.Benchmark instance'


class Divpoly(Benchmark):
    def __init__(self, n):
        """
        Class for benchmarking computation of the division polynomial
        of the following elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.division_polynomial(3)
            3*x^4 + 9*x^3 + 33*x^2 + 87*x + 35
            sage: E.division_polynomial(5)
            5*x^12 + 45*x^11 + 422*x^10 + ... - 426371

            sage: from sage.tests.benchmark import Divpoly
            sage: B = Divpoly(99)
            sage: B
            99-Division polynomial
        """
        self.__n = n
        self.repr_str = "%s-Division polynomial" % self.__n

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import Divpoly
            sage: B = Divpoly(3)
            sage: isinstance(B.sage(), float)
            True
        """
        n = self.__n
        t = cputime()
        E = EllipticCurve([1, 2, 3, 4, 5])
        E.division_polynomial(n)
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import Divpoly
            sage: B = Divpoly(3)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        n = self.__n
        t = magma.cputime()
        magma('DivisionPolynomial(EllipticCurve([1,2,3,4,5]), %s)' % n)
        return magma.cputime(t)


class PolySquare(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R
        self.repr_str = 'Square a polynomial of degree %s over %s' % (self.__n, self.__R)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import PolySquare
            sage: B = PolySquare(3, QQ)
            sage: isinstance(B.sage(), float)
            True
        """
        R = self.__R
        n = self.__n
        f = R['x'](range(1, n + 1))
        t = cputime()
        f**2
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import PolySquare
            sage: B = PolySquare(3, QQ)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma(self.__R)
        f = magma('PolynomialRing(%s)![1..%s]' % (R.name(), self.__n))
        t = magma.cputime()
        f * f
        return magma.cputime(t)

    def maple(self):
        """
        Time the computation in Maple.

        EXAMPLES::

            sage: from sage.tests.benchmark import PolySquare
            sage: B = PolySquare(3, QQ)
            sage: isinstance(B.maple()[1], float) # optional - maple
            True
        """
        R = self.__R
        if not (R == ZZ or R == QQ):
            raise NotImplementedError
        n = self.__n
        f = maple(str(R['x'](range(1, n + 1))))
        t = walltime()
        f * f
        return False, walltime(t)


class MPolynomialPower(Benchmark):
    def __init__(self, nvars=2, exp=10, base=QQ, allow_singular=True):
        self.nvars = nvars
        self.exp = exp
        self.base = base
        self.allow_singular = allow_singular
        s = 'Compute (x_0 + ... + x_%s)^%s over %s' % (
            self.nvars - 1, self.exp, self.base)
        if self.allow_singular:
            s += ' (use singular for Sage mult.)'
        self.repr_str = s

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialPower
            sage: B = MPolynomialPower()
            sage: isinstance(B.sage()[1], float)
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        z = sum(R.gens())
        if self.allow_singular:
            z = singular(sum(R.gens()))
            t = walltime()
            z**self.exp
            return False, walltime(t)
        t = cputime()
        z**self.exp
        return cputime(t)

    def macaulay2(self):
        """
        Time the computation in Macaulay2.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialPower
            sage: B = MPolynomialPower()
            sage: isinstance(B.macaulay2()[1], float) # optional - macaulay2
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        z = macaulay2(sum(R.gens()))
        t = walltime()
        z**self.exp
        return False, walltime(t)

    def maxima(self):
        """
        Time the computation in Maxima.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialPower
            sage: B = MPolynomialPower()
            sage: isinstance(B.maxima()[1], float)
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        z = maxima(str(sum(R.gens())))
        w = walltime()
        (z**self.exp).expand()
        return False, walltime(w)

    def maple(self):
        """
        Time the computation in Maple.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialPower
            sage: B = MPolynomialPower()
            sage: isinstance(B.maple()[1], float)  # optional - maple
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        z = maple(str(sum(R.gens())))
        w = walltime()
        (z**self.exp).expand()
        return False, walltime(w)

    def mathematica(self):
        """
        Time the computation in Mathematica.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialPower
            sage: B = MPolynomialPower()
            sage: isinstance(B.mathematica()[1], float) # optional - mathematica
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        z = mathematica(str(sum(R.gens())))
        w = walltime()
        (z**self.exp).Expand()
        return False, walltime(w)

# this doesn't really expand out -- pari has no function to do so,
# as far as I know.
#     def gp(self):
#         R = PolynomialRing(self.base, self.nvars)
#         z = gp(str(sum(R.gens())))
#         gp.eval('gettime')
#         z**self.exp
#         return float(gp.eval('gettime/1000.0'))

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialPower
            sage: B = MPolynomialPower()
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma.PolynomialRing(self.base, self.nvars)
        z = R.gen(1)
        for i in range(2, self.nvars + 1):
            z += R.gen(i)
        t = magma.cputime()
        z**magma(self.exp)
        return magma.cputime(t)


class MPolynomialMult(Benchmark):
    def __init__(self, nvars=2, base=QQ, allow_singular=True):
        if nvars % 2:
            nvars += 1
        self.nvars = nvars
        self.base = base
        self.allow_singular = allow_singular
        s = 'Compute (x_0 + ... + x_%s) * (x_%s + ... + x_%s) over %s' % (
            self.nvars // 2 - 1, self.nvars // 2, self.nvars, self.base)
        if self.allow_singular:
            s += ' (use singular for Sage mult.)'
        self.repr_str = s

    def maxima(self):
        """
        Time the computation in Maxima.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult
            sage: B = MPolynomialMult()
            sage: isinstance(B.maxima()[1], float)
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = maxima(str(sum(R.gens()[:k])))
        z1 = maxima(str(sum(R.gens()[k:])))
        w = walltime()
        (z0 * z1).expand()
        return False, walltime(w)

    def maple(self):
        """
        Time the computation in Maple.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult
            sage: B = MPolynomialMult()
            sage: isinstance(B.maple()[1], float)  # optional - maple
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = maple(str(sum(R.gens()[:k])))
        z1 = maple(str(sum(R.gens()[k:])))
        w = walltime()
        (z0 * z1).expand()
        return False, walltime(w)

    def mathematica(self):
        """
        Time the computation in Mathematica.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult
            sage: B = MPolynomialMult()
            sage: isinstance(B.mathematica()[1], float) # optional - mathematica
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = mathematica(str(sum(R.gens()[:k])))
        z1 = mathematica(str(sum(R.gens()[k:])))
        w = walltime()
        (z0 * z1).Expand()
        return False, walltime(w)

#     def gp(self):
#         R = PolynomialRing(self.base, self.nvars)
#         k = self.nvars // 2
#         z0 = gp(str(sum(R.gens()[:k])))
#         z1 = gp(str(sum(R.gens()[k:])))
#         gp.eval('gettime')
#         z0*z1
#         return float(gp.eval('gettime/1000.0'))

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult
            sage: B = MPolynomialMult()
            sage: isinstance(B.sage()[1], float)
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = sum(R.gens()[:k])
        z1 = sum(R.gens()[k:])
        if self.allow_singular:
            z0 = singular(z0)
            z1 = singular(z1)
            t = walltime()
            z0 * z1
            return False, walltime(t)
        else:
            t = cputime()
            z0 * z1
            return cputime(t)

    def macaulay2(self):
        """
        Time the computation in Macaulay2.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult
            sage: B = MPolynomialMult()
            sage: isinstance(B.macaulay2()[1], float) # optional - macaulay2
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = macaulay2(sum(R.gens()[:k]))
        z1 = macaulay2(sum(R.gens()[k:]))
        t = walltime()
        z0 * z1
        return False, walltime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult
            sage: B = MPolynomialMult()
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma.PolynomialRing(self.base, self.nvars)
        z0 = R.gen(1)
        k = self.nvars // 2
        for i in range(2, k + 1):
            z0 += R.gen(i)
        z1 = R.gen(k + 1)
        for i in range(k + 1, self.nvars + 1):
            z1 += R.gen(i)
        t = magma.cputime()
        z0 * z1
        return magma.cputime(t)


class MPolynomialMult2(Benchmark):
    def __init__(self, nvars=2, base=QQ, allow_singular=True):
        if nvars % 2:
            nvars += 1
        self.nvars = nvars
        self.base = base
        self.allow_singular = allow_singular
        s = 'Compute (x_1 + 2*x_2 + 3*x_3 + ... + %s*x_%s) * (%s * x_%s + ... + %s*x_%s) over %s' % (
            self.nvars // 2, self.nvars // 2, self.nvars // 2 + 1, self.nvars // 2 + 1,
            self.nvars + 1, self.nvars + 1, self.base)
        if self.allow_singular:
            s += ' (use singular for Sage mult.)'
        self.repr_str = s

#     def gp(self):
#         R = PolynomialRing(self.base, self.nvars)
#         k = self.nvars // 2
#         z0 = R(0)
#         z1 = R(0)
#         for i in range(k):
#             z0 += (i+1)*R.gen(i)
#         for i in range(k,self.nvars):
#             z1 += (i+1)*R.gen(i)
#         z0 = gp(str(z0))
#         z1 = gp(str(z1))
#         gp.eval('gettime')
#         z0*z1
#         return float(gp.eval('gettime/1000.0'))

    def maxima(self):
        """
        Time the computation in Maxima.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult2
            sage: B = MPolynomialMult2()
            sage: isinstance(B.maxima()[1], float)
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i + 1) * R.gen(i)
        for i in range(k, self.nvars):
            z1 += (i + 1) * R.gen(i)
        z0 = maxima(str(z0))
        z1 = maxima(str(z1))
        w = walltime()
        (z0 * z1).expand()
        return False, walltime(w)

    def macaulay2(self):
        """
        Time the computation in Macaulay2.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult2
            sage: B = MPolynomialMult2()
            sage: isinstance(B.macaulay2()[1], float) # optional - macaulay2
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i + 1) * R.gen(i)
        for i in range(k, self.nvars):
            z1 += (i + 1) * R.gen(i)
        z0 = macaulay2(z0)
        z1 = macaulay2(z1)
        t = walltime()
        z0 * z1
        return False, walltime(t)

    def maple(self):
        """
        Time the computation in Maple.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult2
            sage: B = MPolynomialMult2()
            sage: isinstance(B.maple()[1], float) # optional - maple
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i + 1) * R.gen(i)
        for i in range(k, self.nvars):
            z1 += (i + 1) * R.gen(i)
        z0 = maple(str(z0))
        z1 = maple(str(z1))
        w = walltime()
        (z0 * z1).expand()
        return False, walltime(w)

    def mathematica(self):
        """
        Time the computation in Mathematica.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult2
            sage: B = MPolynomialMult2()
            sage: isinstance(B.mathematica()[1], float) # optional - mathematica
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i + 1) * R.gen(i)
        for i in range(k, self.nvars):
            z1 += (i + 1) * R.gen(i)
        z0 = mathematica(str(z0))
        z1 = mathematica(str(z1))
        w = walltime()
        (z0 * z1).Expand()
        return False, walltime(w)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult2
            sage: B = MPolynomialMult2()
            sage: isinstance(B.sage()[1], float)
            True
        """
        R = PolynomialRing(self.base, self.nvars, 'x')
        k = self.nvars // 2
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i + 1) * R.gen(i)
        for i in range(k, self.nvars):
            z1 += (i + 1) * R.gen(i)
        if self.allow_singular:
            z0 = singular(z0)
            z1 = singular(z1)
            t = walltime()
            z0 * z1
            return False, walltime(t)
        else:
            t = cputime()
            z0 * z1
            return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import MPolynomialMult2
            sage: B = MPolynomialMult2()
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma.PolynomialRing(self.base, self.nvars)
        z0 = R.gen(1)
        k = self.nvars // 2
        for i in range(2, k + 1):
            z0 += magma(i) * R.gen(i)
        z1 = R.gen(k + 1)
        for i in range(k + 1, self.nvars + 1):
            z1 += magma(i) * R.gen(i)
        t = magma.cputime()
        z0 * z1
        return magma.cputime(t)


class CharPolyTp(Benchmark):
    def __init__(self, N=37, k=2, p=2, sign=1):
        self.N = N
        self.k = k
        self.p = p
        self.sign = sign
        self.repr_str = "Compute the charpoly (given the matrix) of T_%s on S_%s(Gamma_0(%s)) with sign %s." % (self.p, self.k, self.N, self.sign)

    def matrix(self):
        try:
            return self._matrix
        except AttributeError:
            self._matrix = ModularSymbols(group=self.N, weight=self.k, sign=self.sign).T(self.p).matrix()
        return self._matrix

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import CharPolyTp
            sage: B = CharPolyTp()
            sage: isinstance(B.sage(), float)
            True
        """
        m = self.matrix()
        t = cputime()
        m.charpoly('x')
        return cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import CharPolyTp
            sage: B = CharPolyTp()
            sage: isinstance(B.gp(), float)
            True
        """
        m = gp(self.matrix())
        gp.eval('gettime')
        m.charpoly('x')
        return float(gp.eval('gettime/1000.0'))

    def pari(self):
        """
        Time the computation in Pari.

        EXAMPLES::

            sage: from sage.tests.benchmark import CharPolyTp
            sage: B = CharPolyTp()
            sage: isinstance(B.pari(), float)
            True
        """
        m = pari(self.matrix())
        t = cputime()
        m.charpoly('x')
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import CharPolyTp
            sage: B = CharPolyTp()
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        m = magma(self.matrix())
        t = magma.cputime()
        m.CharacteristicPolynomial()
        return magma.cputime(t)


class PolyFactor(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R
        self.repr_str = "Factor a product of 2 polynomials of degree %s over %s." % (self.__n, self.__R)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import PolyFactor
            sage: B = PolyFactor(3, QQ)
            sage: isinstance(B.sage(), float)
            True
        """
        R = PolynomialRing(self.__R, 'x')
        f = R(range(1, self.__n + 1))
        g = R(range(self.__n + 1, 2 * (self.__n + 1)))
        h = f * g
        t = cputime()
        h.factor()
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import PolyFactor
            sage: B = PolyFactor(3, QQ)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma(self.__R)
        f = magma('PolynomialRing(%s)![1..%s]' % (R.name(), self.__n))
        g = magma('PolynomialRing(%s)![%s+1..2*(%s+1)]' % (
            R.name(), self.__n, self.__n))
        h = f * g
        t = magma.cputime()
        h.Factorization()
        return magma.cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import PolyFactor
            sage: B = PolyFactor(3, QQ)
            sage: isinstance(B.gp(), float)
            True
        """
        R = PolynomialRing(self.__R, 'x')
        f = R(range(1, self.__n + 1))
        g = R(range(self.__n + 1, 2 * (self.__n + 1)))
        h = f * g
        f = gp(h)
        gp.eval('gettime')
        f.factor()
        return float(gp.eval('gettime/1000.0'))


class SquareInts(Benchmark):
    def __init__(self, base=10, ndigits=10**5):
        self.__ndigits = ndigits
        self.base = base
        self.repr_str = "Square the integer %s^%s" % (self.base, self.__ndigits)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.sage(), float)
            True
        """
        n = Integer(self.base)**self.__ndigits
        t = cputime()
        n**2
        return cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.gp(), float)
            True
        """
        n = gp('%s^%s' % (self.base, self.__ndigits))
        gp.eval('gettime')
        n**2
        return float(gp.eval('gettime/1000.0'))

    def maxima(self):
        """
        Time the computation in Maxima.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.maxima()[1], float)
            True
        """
        n = maxima('%s^%s' % (self.base, self.__ndigits))
        t = walltime()
        n**2
        return False, walltime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        n = magma('%s^%s' % (self.base, self.__ndigits))
        t = magma.cputime()
        n**2
        return magma.cputime(t)

    def python(self):
        """
        Time the computation in Python.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.python(), float)
            True
        """
        n = self.base**self.__ndigits
        t = cputime()
        n**2
        return cputime(t)

    def maple(self):
        """
        Time the computation in Maple.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.maple()[1], float) # optional - maple
            True
        """
        n = maple('%s^%s' % (self.base, self.__ndigits))
        t = walltime()
        n**2
        return False, walltime(t)

    def libgap(self):
        """
        Time the computation in libGAP.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.libgap()[1], float)
            True
        """
        n = libgap(self.base)**libgap(self.__ndigits)
        t = walltime()
        n**2
        return False, walltime(t)

    def mathematica(self):
        """
        Time the computation in Mathematica.

        EXAMPLES::

            sage: from sage.tests.benchmark import SquareInts
            sage: B = SquareInts()
            sage: isinstance(B.mathematica()[1], float) # optional - mathematica
            True
        """
        n = mathematica('%s^%s' % (self.base, self.__ndigits))
        t = walltime()
        n**2
        return False, walltime(t)


class MatrixSquare(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R
        self.repr_str = 'Square a matrix of degree %s over %s' % (self.__n, self.__R)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixSquare
            sage: B = MatrixSquare(3, QQ)
            sage: isinstance(B.sage(), float)
            True
        """
        R = self.__R
        n = self.__n
        f = MatrixSpace(R, n)(list(range(n * n)))
        t = cputime()
        f**2
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixSquare
            sage: B = MatrixSquare(3, QQ)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma(self.__R)
        f = magma('MatrixAlgebra(%s, %s)![0..%s^2-1]' % (
            R.name(), self.__n, self.__n))
        t = magma.cputime()
        f * f
        return magma.cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixSquare
            sage: B = MatrixSquare(3, QQ)
            sage: isinstance(B.gp(), float)
            True
        """
        n = self.__n
        m = gp('matrix(%s,%s,m,n,%s*(m-1)+(n-1))' % (n, n, n))
        gp('gettime')
        m * m
        return float(gp.eval('gettime/1000.0'))

    def libgap(self):
        """
        Time the computation in libGAP.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixSquare
            sage: B = MatrixSquare(3, QQ)
            sage: isinstance(B.libgap()[1], float)
            True
        """
        R = self.__R
        n = self.__n
        f = MatrixSpace(R, n)(list(range(n * n)))
        m = libgap(f)
        t = walltime()
        m * m
        return False, walltime(t)


class Factorial(Benchmark):
    def __init__(self, n):
        self.__n = n
        self.repr_str = "Compute the factorial of %s" % self.__n

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import Factorial
            sage: B = Factorial(10)
            sage: isinstance(B.sage(), float)
            True
        """
        t = cputime()
        factorial(self.__n)
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import Factorial
            sage: B = Factorial(10)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        t = magma.cputime()
        magma('&*[1..%s]' % self.__n)  # &* is way better than Factorial!!
        return magma.cputime(t)

    def maple(self):
        """
        Time the computation in Maple.

        EXAMPLES::

            sage: from sage.tests.benchmark import Factorial
            sage: B = Factorial(10)
            sage: isinstance(B.maple()[1], float) # optional - maple
            True
        """
        n = maple(self.__n)
        t = walltime()
        n.factorial()
        return False, walltime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import Factorial
            sage: B = Factorial(10)
            sage: isinstance(B.gp(), float)
            True
        """
        gp.eval('gettime')
        gp('%s!' % self.__n)
        return float(gp.eval('gettime/1000.0'))


class Fibonacci(Benchmark):
    def __init__(self, n):
        self.__n = n
        self.repr_str = "Compute the %s-th Fibonacci number" % self.__n

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import Fibonacci
            sage: B = Fibonacci(10)
            sage: isinstance(B.sage(), float)
            True
        """
        t = cputime()
        fibonacci(self.__n)
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import Fibonacci
            sage: B = Fibonacci(10)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        t = magma.cputime()
        magma('Fibonacci(%s)' % self.__n)
        return magma.cputime(t)

    def libgap(self):
        """
        Time the computation in libGAP.

        EXAMPLES::

            sage: from sage.tests.benchmark import Fibonacci
            sage: B = Fibonacci(10)
            sage: isinstance(B.libgap()[1], float)
            True
        """
        n = libgap(self.__n)
        t = walltime()
        n.Fibonacci()
        return False, walltime(t)

    def mathematica(self):
        """
        Time the computation in Mathematica.

        EXAMPLES::

            sage: from sage.tests.benchmark import Fibonacci
            sage: B = Fibonacci(10)
            sage: isinstance(B.mathematica()[1], float) # optional - mathematica
            True
        """
        n = mathematica(self.__n)
        t = walltime()
        n.Fibonacci()
        return False, walltime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import Fibonacci
            sage: B = Fibonacci(10)
            sage: isinstance(B.gp(), float)
            True
        """
        gp.eval('gettime')
        gp('fibonacci(%s)' % self.__n)
        return float(gp.eval('gettime/1000.0'))


class SEA(Benchmark):
    def __init__(self, p):
        self.__p = p
        self.repr_str = "Do SEA on an elliptic curve over GF(%s)" % self.__p

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import SEA
            sage: B = SEA(5)
            sage: isinstance(B.sage()[1], float)
            True
        """
        E = EllipticCurve([1, 2, 3, 4, 5])
        t = walltime()
        # Note that from pari 2.4.3, the SEA algorithm is used by the
        # pari library, but only for large primes, so for a better
        # test a prime > 2^30 should be used and not 5.  In fact
        # next_prime(2^100) works fine (<<1s).
        E.change_ring(GF(self.__p)).cardinality_pari()
        return False, walltime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import SEA
            sage: B = SEA(5)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        magma(0)
        t = magma.cputime()
        magma('#EllipticCurve([GF(%s)|1,2,3,4,5])' % (self.__p))
        return magma.cputime(t)


class MatrixKernel(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R
        self.repr_str = 'Kernel of a matrix of degree %s over %s' % (self.__n, self.__R)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixKernel
            sage: B = MatrixKernel(3, QQ)
            sage: isinstance(B.sage(), float)
            True
        """
        R = self.__R
        n = self.__n
        f = MatrixSpace(R, n, 2 * n)(list(range(n * (2 * n))))
        t = cputime()
        f.kernel()
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixKernel
            sage: B = MatrixKernel(3, QQ)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        R = magma(self.__R)
        f = magma('RMatrixSpace(%s, %s, %s)![0..(%s*2*%s)-1]' % (
            R.name(), self.__n, 2 * self.__n, self.__n, self.__n))
        t = magma.cputime()
        f.Kernel()
        return magma.cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import MatrixKernel
            sage: B = MatrixKernel(3, QQ)
            sage: isinstance(B.gp(), float)
            True
        """
        n = self.__n
        m = gp('matrix(%s,%s,m,n,%s*(m-1)+(n-1))' % (n, 2 * n, n))
        gp('gettime')
        m.matker()
        return float(gp.eval('gettime/1000.0'))


class ComplexMultiply(Benchmark):
    def __init__(self, bits_prec, times):
        self.__bits_prec = bits_prec
        self.__times = times
        self.repr_str = "List of multiplies of two complex numbers with %s bits of precision %s times" % (self.__bits_prec, self.__times)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import ComplexMultiply
            sage: B = ComplexMultiply(28, 2)
            sage: isinstance(B.sage(), float)
            True
        """
        CC = ComplexField(self.__bits_prec)
        s = CC(2).sqrt() + (CC.gen() * 2).sqrt()
        t = cputime()
        [s * s for _ in range(self.__times)]
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import ComplexMultiply
            sage: B = ComplexMultiply(28, 2)
            sage: isinstance(B.magma(), float) # optional - magma
            True

        .. NOTE::

            decimal digits (despite magma docs that say bits!!)
        """
        n = int(self.__bits_prec / log(10, 2)) + 1
        CC = magma.ComplexField(n)
        s = CC(2).Sqrt() + CC.gen(1).Sqrt()
        t = magma.cputime()
        magma.eval('s := %s;' % s.name())
        magma('[s*s : i in [1..%s]]' % self.__times)
        return magma.cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import ComplexMultiply
            sage: B = ComplexMultiply(28, 2)
            sage: isinstance(B.gp(), float)
            True
        """
        n = int(self.__bits_prec / log(10, 2)) + 1
        gp.set_real_precision(n)
        gp.eval('s = sqrt(2) + sqrt(2*I);')
        gp.eval('gettime;')
        gp('vector(%s,i,s*s)' % self.__times)
        return float(gp.eval('gettime/1000.0'))


class ModularSymbols1(Benchmark):
    def __init__(self, N, k=2):
        self.__N = N
        self.__k = k
        self.repr_str = 'Presentation for modular symbols on Gamma_0(%s) of weight %s' % (self.__N, self.__k)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import ModularSymbols1
            sage: B = ModularSymbols1(11)
            sage: isinstance(B.sage(), float)
            True
        """
        t = cputime()
        ModularSymbols(self.__N, self.__k)
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import ModularSymbols1
            sage: B = ModularSymbols1(11)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        magma = Magma()  # new instance since otherwise modsyms are cached, and cache can't be cleared
        t = magma.cputime()
        magma('ModularSymbols(%s, %s)' % (self.__N, self.__k))
        return magma.cputime(t)


class ModularSymbolsDecomp1(Benchmark):
    def __init__(self, N, k=2, sign=1, bnd=10):
        self.N = N
        self.k = k
        self.sign = sign
        self.bnd = bnd
        self.repr_str = 'Decomposition of modular symbols on Gamma_0(%s) of weight %s and sign %s' % (self.N, self.k, self.sign)

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import ModularSymbolsDecomp1
            sage: B = ModularSymbolsDecomp1(11)
            sage: isinstance(B.sage(), float)
            True
        """
        t = cputime()
        M = ModularSymbols(self.N, self.k, sign=self.sign, use_cache=False)
        M.decomposition(self.bnd)
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import ModularSymbolsDecomp1
            sage: B = ModularSymbolsDecomp1(11)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        m = Magma()  # new instance since otherwise modsyms are cached, and cache can't be cleared
        t = m.cputime()
        m.eval('Decomposition(ModularSymbols(%s, %s, %s),%s);' % (
            self.N, self.k, self.sign, self.bnd))
        return m.cputime(t)


class EllipticCurveTraces(Benchmark):
    def __init__(self, B):
        self.B = B
        self.repr_str = "Compute all a_p for the elliptic curve [1,2,3,4,5], for p < %s" % self.B

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurveTraces
            sage: B = EllipticCurveTraces(11)
            sage: isinstance(B.sage(), float)
            Traceback (most recent call last):
            ...
            TypeError: ...anlist() got an unexpected keyword argument 'pari_ints'
        """
        E = EllipticCurve([1, 2, 3, 4, 5])
        t = cputime()
        E.anlist(self.B, pari_ints=True)
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurveTraces
            sage: B = EllipticCurveTraces(11)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        E = magma.EllipticCurve([1, 2, 3, 4, 5])
        t = magma.cputime()
        E.TracesOfFrobenius(self.B)
        return magma.cputime(t)


class EllipticCurvePointMul(Benchmark):
    def __init__(self, n):
        self.n = n
        self.repr_str = "Compute %s*(0,0) on the elliptic curve [0, 0, 1, -1, 0] over QQ" % self.n

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurvePointMul
            sage: B = EllipticCurvePointMul(11)
            sage: isinstance(B.sage(), float)
            True
        """
        E = EllipticCurve([0, 0, 1, -1, 0])
        P = E([0, 0])
        t = cputime()
        self.n * P
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurvePointMul
            sage: B = EllipticCurvePointMul(11)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        E = magma.EllipticCurve('[0, 0, 1, -1, 0]')
        P = E('[0,0]')
        t = magma.cputime()
        magma(self.n) * P
        return magma.cputime(t)

    def gp(self):
        """
        Time the computation in GP.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurvePointMul
            sage: B = EllipticCurvePointMul(11)
            sage: isinstance(B.gp(), float)
            True
        """
        E = gp.ellinit('[0, 0, 1, -1, 0]')
        gp.eval('gettime')
        P = gp([0, 0])
        E.ellmul(P, self.n)
        return float(gp.eval('gettime/1000.0'))

    def pari(self):
        """
        Time the computation in Pari.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurvePointMul
            sage: B = EllipticCurvePointMul(11)
            sage: isinstance(B.pari(), float)
            True
        """
        E = pari('ellinit([0, 0, 1, -1, 0])')
        pari('gettime')
        P = pari([0, 0])
        E.ellmul(P, self.n)
        return float(pari('gettime/1000.0'))


class EllipticCurveMW(Benchmark):
    def __init__(self, ainvs):
        self.ainvs = ainvs
        self.repr_str = "Compute generators for the Mordell-Weil group of the elliptic curve %s over QQ" % self.ainvs

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurveMW
            sage: B = EllipticCurveMW([1,2,3,4,5])
            sage: isinstance(B.sage()[1], float)
            True
        """
        E = EllipticCurve(self.ainvs)
        t = walltime()
        E.gens()
        return False, walltime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import EllipticCurveMW
            sage: B = EllipticCurveMW([1,2,3,4,5])
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        E = magma.EllipticCurve(str(self.ainvs))
        t = magma.cputime()
        E.Generators()
        return magma.cputime(t)


class FiniteExtFieldMult(Benchmark):
    def __init__(self, field, times):
        self.__times = times
        self.field = field
        self.e = field.gen()**(field.cardinality() / 3)
        self.f = field.gen()**(2 * field.cardinality() / 3)
        self.repr_str = "Multiply a^(#K/3) with a^(2*#K/3) where a == K.gen()"

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import FiniteExtFieldMult
            sage: B = FiniteExtFieldMult(GF(9, 'x'), 2)
            sage: isinstance(B.sage(), float)
            True
        """
        e = self.e
        f = self.f
        t = cputime()
        [e * f for _ in range(self.__times)]
        return cputime(t)

    def pari(self):
        """
        Time the computation in Pari.

        EXAMPLES::

            sage: from sage.tests.benchmark import FiniteExtFieldMult
            sage: B = FiniteExtFieldMult(GF(9, 'x'), 2)
            sage: isinstance(B.pari(), float)
            True
        """
        e = self.e.__pari__()
        f = self.f.__pari__()
        t = cputime()
        [e * f for _ in range(self.__times)]
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import FiniteExtFieldMult
            sage: B = FiniteExtFieldMult(GF(9, 'x'), 2)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        magma.eval('F<a> := GF(%s)' % (self.field.cardinality()))
        magma.eval('e := a^Floor(%s/3);' % (self.field.cardinality()))
        magma.eval('f := a^Floor(2*%s/3);' % (self.field.cardinality()))
        t = magma.cputime()
        magma('[e*f : i in [1..%s]]' % self.__times)
        return magma.cputime(t)


class FiniteExtFieldAdd(Benchmark):
    def __init__(self, field, times):
        self.__times = times
        self.field = field
        self.e = field.gen()**(field.cardinality() / 3)
        self.f = field.gen()**(2 * field.cardinality() / 3)
        self.repr_str = "Add a^(#K/3) to a^(2*#K/3) where a == K.gen()"

    def sage(self):
        """
        Time the computation in Sage.

        EXAMPLES::

            sage: from sage.tests.benchmark import FiniteExtFieldAdd
            sage: B = FiniteExtFieldAdd(GF(9,'x'), 2)
            sage: isinstance(B.sage(), float)
            True
        """
        e = self.e
        f = self.f
        t = cputime()
        [e + f for _ in range(self.__times)]
        return cputime(t)

    def pari(self):
        """
        Time the computation in Pari.

        EXAMPLES::

            sage: from sage.tests.benchmark import FiniteExtFieldAdd
            sage: B = FiniteExtFieldAdd(GF(9,'x'), 2)
            sage: isinstance(B.pari(), float)
            True
        """
        e = self.e.__pari__()
        f = self.f.__pari__()
        t = cputime()
        [e + f for _ in range(self.__times)]
        return cputime(t)

    def magma(self):
        """
        Time the computation in Magma.

        EXAMPLES::

            sage: from sage.tests.benchmark import FiniteExtFieldAdd
            sage: B = FiniteExtFieldAdd(GF(9,'x'), 2)
            sage: isinstance(B.magma(), float) # optional - magma
            True
        """
        magma.eval('F<a> := GF(%s)' % (self.field.cardinality()))
        magma.eval('e := a^Floor(%s/3);' % (self.field.cardinality()))
        magma.eval('f := a^Floor(2*%s/3);' % (self.field.cardinality()))
        t = magma.cputime()
        magma('[e+f : i in [1..%s]]' % self.__times)
        return magma.cputime(t)


"""
.. TODO::

   * multiply reals
   * modular degree
   * anlist
   * MW group
   * symbolic det
   * poly factor
   * multivariate poly factor
"""


def suite1():
    PolySquare(10000, QQ).run()
    PolySquare(20000, ZZ).run()
    PolySquare(50000, GF(5)).run()
    PolySquare(20000, Integers(8)).run()

    SquareInts(10, 2000000).run()

    MatrixSquare(200, QQ).run()
    MatrixSquare(50, ZZ).run()

    SquareInts(10, 150000).run()

    Factorial(2 * 10**6).run(systems=['sage', 'magma'])
    Fibonacci(10**6).run()
    Fibonacci(2 * 10**7).run(systems=["sage", "magma", "mathematica"])

    MatrixKernel(150, QQ).run()

    ComplexMultiply(100000, 1000)
    ComplexMultiply(100, 100000)
    ComplexMultiply(53, 100000)

    PolyFactor(300, ZZ)
    PolyFactor(300, GF(19))
    PolyFactor(700, GF(19))

    PolyFactor(500, GF((9, 2), 'a'))
    PolyFactor(100, GF((10007, 3), 'a'))

    CharPolyTp(54, 4).run()
    CharPolyTp(389, 2).run()
    CharPolyTp(389, 2, sign=0, p=3).run()
    CharPolyTp(1000, 2, sign=1, p=2).run(systems=['sage', 'magma'])
    CharPolyTp(1, 100, sign=1, p=5).run(systems=['sage', 'magma'])   # Sage's multimodular really sucks here! (GP is way better, even)
    CharPolyTp(512, sign=1, p=3).run(systems=['sage', 'magma', 'gp'])
    CharPolyTp(512, sign=0, p=3).run(systems=['sage', 'magma', 'gp'])
    CharPolyTp(1024, sign=1, p=3).run(systems=['sage', 'magma', 'gp'])
    CharPolyTp(2006, sign=1, p=2).run(systems=['sage', 'magma', 'gp'])
    CharPolyTp(2006, sign=1, p=2).run(systems=['sage', 'magma'])    # gp takes > 1 minute.


def mpoly():
    # This includes a maxima benchmark.  Note that
    # maxima is *shockingly* slow in comparison to Singular or MAGMA.
    # It is so slow as to be useless, basically, i.e., factor
    # of 5000 slower than Singular on this example!
    MPolynomialPower(nvars=6, exp=10).run()

    main = ['sage', 'magma']  # just the main competitors
    MPolynomialPower(nvars=2, exp=200, allow_singular=False).run(main)
    MPolynomialPower(nvars=5, exp=10, allow_singular=False).run(main)
    MPolynomialPower(nvars=5, exp=30, allow_singular=True).run(main)
    MPolynomialPower(nvars=2, exp=1000, allow_singular=True).run(main)
    MPolynomialPower(nvars=10, exp=10, allow_singular=True).run(main)
    MPolynomialPower(nvars=4, exp=350, base=GF(7), allow_singular=True).run(main)
    MPolynomialMult(200, allow_singular=False).run(main)
    MPolynomialMult(400, allow_singular=True).run(main)
    MPolynomialMult(800, allow_singular=True).run(main)
    MPolynomialMult2(500, allow_singular=True).run(main)


def mpoly_all(include_maple=False):
    """
    Run benchmarks for multipoly arithmetic on all systems (except
    Maxima, since it is very very slow).  You must have mathematica,
    maple, and magma.

    .. NOTE::

       * maple is depressingly slow on these benchmarks.
       * Singular (i.e., Sage) does shockingly well.
       * mathematica is sometimes amazing.
       * macaulay2 is also quite bad (though not as bad as maple).

    EXAMPLES::

        sage: from sage.tests.benchmark import mpoly_all
        sage: mpoly_all() # not tested
        <BLANKLINE>
        ...
        ...System      min         avg         max         trials          cpu or wall
        ...
        * sage...
    """
    systems = ['sage', 'magma', 'mathematica', 'macaulay2']
    if include_maple:
        systems.append('maple')
    MPolynomialMult(200).run(systems=systems)
    MPolynomialMult(400).run(systems=systems)
    MPolynomialMult2(256).run(systems=systems)
    MPolynomialMult2(512).run(systems=systems)
    MPolynomialPower(nvars=4, exp=50).run(systems=systems)   # mathematica wins
    MPolynomialPower(nvars=10, exp=10).run(systems=systems)


def modsym_present():
    ModularSymbols1(2006, 2)
    ModularSymbols1(1, 50)
    ModularSymbols1(1, 100)
    ModularSymbols1(1, 150)
    ModularSymbols1(30, 8)
    ModularSymbols1(225, 4)
    ModularSymbols1(2, 50)
    ModularSymbols1(2, 100)


def modsym_decomp():
    ModularSymbolsDecomp1(1, 24).run()
    ModularSymbolsDecomp1(125, 2).run()
    ModularSymbolsDecomp1(389, 2).run()
    ModularSymbolsDecomp1(1, 100).run()
    ModularSymbolsDecomp1(54, 4).run()


def elliptic_curve():
    EllipticCurveTraces(100000).run()
    EllipticCurveTraces(500000).run()
    Divpoly(59).run()
    EllipticCurvePointMul(1000).run()
    EllipticCurvePointMul(2000).run()
    EllipticCurvePointMul(2500).run()      # sage is clearly using the wrong algorithm -- maybe need a balanced rep!?

    # NOTE -- Sage can also do these using Simon's program, which is
    # *way* *way* faster than MAGMA...
    EllipticCurveMW([5, 6, 7, 8, 9]).run()
    EllipticCurveMW([50, 6, 7, 8, 9]).run()
    EllipticCurveMW([1, -1, 0, -79, 289]).run(trials=1)   # rank 4
    EllipticCurveMW([0, 0, 1, -79, 342]).run(trials=1)    # rank 5  (Sage wins)
