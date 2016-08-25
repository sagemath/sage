r"""
Piecewise-defined Functions

Sage implements a very simple class of piecewise-defined functions.
Functions may be any type of symbolic expression. Infinite
intervals are not supported. The endpoints of each interval must
line up.

TODO:

- Implement max/min location and values,

- Need: parent object - ring of piecewise functions

- This class should derive from an element-type class, and should
  define ``_add_``, ``_mul_``, etc. That will automatically take care
  of left multiplication and proper coercion. The coercion mentioned
  below for scalar mult on right is bad, since it only allows ints and
  rationals. The right way is to use an element class and only define
  ``_mul_``, and have a parent, so anything gets coerced properly.

AUTHORS:

- David Joyner (2006-04): initial version

- David Joyner (2006-09): added __eq__, extend_by_zero_to, unextend,
  convolution, trapezoid, trapezoid_integral_approximation,
  riemann_sum, riemann_sum_integral_approximation, tangent_line fixed
  bugs in __mul__, __add__

- David Joyner (2007-03): adding Hann filter for FS, added general FS
  filter methods for computing and plotting, added options to plotting
  of FS (eg, specifying rgb values are now allowed). Fixed bug in
  documentation reported by Pablo De Napoli.

- David Joyner (2007-09): bug fixes due to behaviour of
  SymbolicArithmetic

- David Joyner (2008-04): fixed docstring bugs reported by J Morrow; added
  support for Laplace transform of functions with infinite support.

- David Joyner (2008-07): fixed a left multiplication bug reported by
  C. Boncelet (by defining __rmul__ = __mul__).

- Paul Butler (2009-01): added indefinite integration and default_variable

TESTS::

    sage: R.<x> = QQ[]
    sage: f = Piecewise([[(0,1),1*x^0]])
    doctest:...: DeprecationWarning: use lower-case piecewise instead
    See http://trac.sagemath.org/14801 for details.
    sage: 2*f
    doctest:...: DeprecationWarning: use lower-case piecewise instead
    See http://trac.sagemath.org/14801 for details.
    Piecewise defined function with 1 parts, [[(0, 1), 2]]
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdjoyner@gmail.com>
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

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.sage_eval import sage_eval
from sage.rings.all import QQ, RR, Integer, Rational, infinity
from sage.calculus.functional import derivative
from sage.symbolic.expression import is_Expression
from sage.symbolic.assumptions import assume, forget

from sage.calculus.calculus import SR, maxima
from sage.calculus.all import var

def Piecewise(list_of_pairs, var=None):
    """
    Deprecated spelling of :func:`sage.functions.piecewise`.

    Return a piecewise function from a list of (interval, function)
    pairs.
    
    ``list_of_pairs`` is a list of pairs (I, fcn), where
    fcn is a Sage function (such as a polynomial over RR, or functions
    using the lambda notation), and I is an interval such as I = (1,3).
    Two consecutive intervals must share a common endpoint.
    
    If the optional ``var`` is specified, then any symbolic expressions
    in the list will be converted to symbolic functions using 
    ``fcn.function(var)``.  (This says which variable is considered to
    be "piecewise".)

    We assume that these definitions are consistent (ie, no checking is
    done).
    
    EXAMPLES::
    
        sage: f1(x) = -1
        sage: f2(x) = 2
        sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
        doctest:...: DeprecationWarning: use lower-case piecewise instead
        See http://trac.sagemath.org/14801 for details.
        sage: f(1)
        -1
        sage: f(3)
        2
        sage: f = Piecewise([[(0,1),x], [(1,2),x^2]], x); f
        Piecewise defined function with 2 parts, [[(0, 1), x |--> x], [(1, 2), x |--> x^2]]
        sage: f(0.9)
        0.900000000000000
        sage: f(1.1)
        1.21000000000000
    """
    from sage.misc.superseded import deprecation
    deprecation(14801, 'use lower-case piecewise instead')
    return PiecewisePolynomial(list_of_pairs, var=var)

class PiecewisePolynomial:
    """
    Returns a piecewise function from a list of (interval, function)
    pairs.
    
    EXAMPLES::
    
        sage: f1(x) = -1
        sage: f2(x) = 2
        sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
        doctest:...: DeprecationWarning: use lower-case piecewise instead
        See http://trac.sagemath.org/14801 for details.
        sage: f(1)
        -1
        sage: f(3)
        2
    """
    def __init__(self, list_of_pairs, var=None):
        r"""
        ``list_of_pairs`` is a list of pairs (I, fcn), where
        fcn is a Sage function (such as a polynomial over RR, or functions
        using the lambda notation), and I is an interval such as I = (1,3).
        Two consecutive intervals must share a common endpoint.
        
        If the optional ``var`` is specified, then any symbolic
        expressions in the list will be converted to symbolic
        functions using ``fcn.function(var)``.  (This says which
        variable is considered to be "piecewise".)

        We assume that these definitions are consistent (ie, no checking is
        done).

        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.list()
            [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
            sage: f.length()
            2
        """
        self._length = len(list_of_pairs)
        self._intervals = [x[0] for x in list_of_pairs]
        functions = [x[1] for x in list_of_pairs]
        if var is not None:
            for i in range(len(functions)):
                if is_Expression(functions[i]):
                    functions[i] = functions[i].function(var)
        self._functions = functions
        # We regenerate self._list in case self._functions was modified
        # above.  This also protects us in case somebody mutates a list
        # after they use it as an argument to piecewise().
        self._list = [[self._intervals[i], self._functions[i]] for i in range(self._length)]
 
    def list(self):
        """
        Returns the pieces of this function as a list of functions.
        
        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.list()
            [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
        """
        return self._list
 
    def length(self):
        """
        Returns the number of pieces of this function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.length()
            2
        """
        return self._length
 
    def __repr__(self):
        """
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]]); f
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 2 parts, [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
        """
        return 'Piecewise defined function with %s parts, %s'%(
            self.length(),self.list())
 
    def _latex_(self):
        r"""
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: latex(f)
            \begin{cases}
            x \ {\mapsto}\ 1 &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ -x + 1 &\text{on $(1, 2)$}\cr
            \end{cases}
        
        ::
        
            sage: f(x) = sin(x*pi/2)
            sage: g(x) = 1-(x-1)^2
            sage: h(x) = -x
            sage: P = Piecewise([[(0,1), f], [(1,3),g], [(3,5), h]])
            sage: latex(P)
            \begin{cases}
            x \ {\mapsto}\ \sin\left(\frac{1}{2} \, \pi x\right) &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ -{\left(x - 1\right)}^{2} + 1 &\text{on $(1, 3)$}\cr
            x \ {\mapsto}\ -x &\text{on $(3, 5)$}\cr
            \end{cases}
        """
        from sage.misc.latex import latex
        tex = ['\\begin{cases}\n']
        for (left, right), f in self.list():
            tex.append('%s &\\text{on $(%s, %s)$}\\cr\n' % (latex(f), left, right))
        tex.append(r'\end{cases}')
        return ''.join(tex)

    def intervals(self):
        """
        A piecewise non-polynomial example.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.intervals()
            [(0, 1), (1, 2), (2, 3), (3, 10)]
        """
        return self._intervals
 
    def domain(self):
        """
        Returns the domain of the function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.domain()
            (0, 10)
        """
        endpoints = sum(self.intervals(), ())
        return (min(endpoints), max(endpoints))

    def functions(self):
        """
        Returns the list of functions (the "pieces").
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.functions()
            [x |--> 1, x |--> -x + 1, x |--> e^x, x |--> sin(2*x)]
        """
        return self._functions
        
    def extend_by_zero_to(self,xmin=-1000,xmax=1000):
        """
        This function simply returns the piecewise defined function which
        is extended by 0 so it is defined on all of (xmin,xmax). This is
        needed to add two piecewise functions in a reasonable way.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.extend_by_zero_to(-1, 3)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 4 parts, [[(-1, 0), 0], [(0, 1), x |--> 1], [(1, 2), x |--> -x + 1], [(2, 3), 0]]
        """
        zero = QQ['x'](0)
        list_of_pairs = self.list()
        a, b = self.domain()
        if xmin < a:
            list_of_pairs = [[(xmin, a), zero]] + list_of_pairs
        if xmax > b:
            list_of_pairs = list_of_pairs + [[(b, xmax), zero]]
        return Piecewise(list_of_pairs)

    def unextend(self):
        """
        This removes any parts in the front or back of the function which
        is zero (the inverse to extend_by_zero_to).
        
        EXAMPLES::
        
            sage: R.<x> = QQ[]
            sage: f = Piecewise([[(-3,-1),1+2+x],[(-1,1),1-x^2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: e = f.extend_by_zero_to(-10,10); e
            Piecewise defined function with 4 parts, [[(-10, -3), 0], [(-3, -1), x + 3], [(-1, 1), -x^2 + 1], [(1, 10), 0]]
            sage: d = e.unextend(); d
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 2 parts, [[(-3, -1), x + 3], [(-1, 1), -x^2 + 1]]
            sage: d==f
            True
        """
        list_of_pairs = self.list()
        funcs = self.functions()
        if funcs[0] == 0:
            list_of_pairs = list_of_pairs[1:]
        if funcs[-1] == 0:
            list_of_pairs = list_of_pairs[:-1]
        return Piecewise(list_of_pairs)

    def _riemann_sum_helper(self, N, func, initial=0):
        """
        A helper function for computing Riemann sums.
        
        INPUT:
        
        
        -  ``N`` - the number of subdivisions
        
        -  ``func`` - a function to apply to the endpoints of
           each subdivision
        
        -  ``initial`` - the starting value
        
        
        EXAMPLES::
        
            sage: f1(x) = x^2                   ## example 1
            sage: f2(x) = 5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f._riemann_sum_helper(6, lambda x0, x1: (x1-x0)*f(x1))
            19/6
        """
        a,b = self.domain()
        rsum = initial
        h = (b-a)/N
        for i in range(N):
            x0 = a+i*h
            x1 = a+(i+1)*h
            rsum += func(x0, x1)
        return rsum

    def riemann_sum_integral_approximation(self,N,mode=None):
        """
        Returns the piecewise line function defined by the Riemann sums in
        numerical integration based on a subdivision into N subintervals.
        
        Set mode="midpoint" for the height of the rectangles to be
        determined by the midpoint of the subinterval; set mode="right" for
        the height of the rectangles to be determined by the right-hand
        endpoint of the subinterval; the default is mode="left" (the height
        of the rectangles to be determined by the left-hand endpoint of
        the subinterval).
        
        EXAMPLES::
        
            sage: f1(x) = x^2                   ## example 1
            sage: f2(x) = 5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.riemann_sum_integral_approximation(6)
            17/6
            sage: f.riemann_sum_integral_approximation(6,mode="right")
            19/6
            sage: f.riemann_sum_integral_approximation(6,mode="midpoint")
            3
            sage: f.integral(definite=True)
            3
        """
        if mode is None:
            return self._riemann_sum_helper(N, lambda x0, x1: (x1-x0)*self(x0))
        elif mode == "right":
            return self._riemann_sum_helper(N, lambda x0, x1: (x1-x0)*self(x1))
        elif mode == "midpoint":
            return self._riemann_sum_helper(N, lambda x0, x1: (x1-x0)*self((x0+x1)/2))
        else:
            raise ValueError("invalid mode")

    def riemann_sum(self,N,mode=None):
        """
        Returns the piecewise line function defined by the Riemann sums in
        numerical integration based on a subdivision into N subintervals.
        Set mode="midpoint" for the height of the rectangles to be
        determined by the midpoint of the subinterval; set mode="right" for
        the height of the rectangles to be determined by the right-hand
        endpoint of the subinterval; the default is mode="left" (the height
        of the rectangles to be determined by the left-hand endpoint of
        the subinterval).
        
        EXAMPLES::
        
            sage: f1(x) = x^2
            sage: f2(x) = 5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.riemann_sum(6,mode="midpoint")
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 6 parts, [[(0, 1/3), 1/36], [(1/3, 2/3), 1/4], [(2/3, 1), 25/36], [(1, 4/3), 131/36], [(4/3, 5/3), 11/4], [(5/3, 2), 59/36]]
        
        ::
        
            sage: f = Piecewise([[(-1,1),(1-x^2).function(x)]])
            sage: rsf = f.riemann_sum(7)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = rsf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[a,0],[a,f(x=a)]],rgbcolor=(0.7,0.6,0.6)) for (a,b),f in rsf.list()])
            sage: P + Q + L
            Graphics object consisting of 15 graphics primitives

        ::
        
            sage: f = Piecewise([[(-1,1),(1/2+x-x^3)]], x) ## example 3
            sage: rsf = f.riemann_sum(8)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = rsf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[a,0],[a,f(x=a)]],rgbcolor=(0.7,0.6,0.6)) for (a,b),f in rsf.list()])
            sage: P + Q + L
            Graphics object consisting of 17 graphics primitives
            """
        if mode is None:
            rsum = self._riemann_sum_helper(N, lambda x0,x1: [[(x0,x1),SR(self(x0))]],
                                            initial=[])
        elif mode == "right":
            rsum = self._riemann_sum_helper(N, lambda x0,x1: [[(x0,x1),SR(self(x1))]],
                                            initial=[])
        elif mode == "midpoint":
            rsum = self._riemann_sum_helper(N, lambda x0,x1: [[(x0,x1),SR(self((x0+x1)/2))]],
                                            initial=[])
        else:
            raise ValueError("invalid mode")
        return Piecewise(rsum)

    def trapezoid(self,N):
        """
        Returns the piecewise line function defined by the trapezoid rule
        for numerical integration based on a subdivision into N
        subintervals.
        
        EXAMPLES::
        
            sage: R.<x> = QQ[]
            sage: f1 = x^2 
            sage: f2 = 5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.trapezoid(4)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 4 parts, [[(0, 1/2), 1/2*x], [(1/2, 1), 9/2*x - 2], [(1, 3/2), 1/2*x + 2], [(3/2, 2), -7/2*x + 8]]
        
        ::
        
            sage: R.<x> = QQ[]
            sage: f = Piecewise([[(-1,1),1-x^2]])
            sage: tf = f.trapezoid(4)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[a,0],[a,f(a)]],rgbcolor=(0.7,0.6,0.6)) for (a,b),f in tf.list()])
            sage: P+Q+L
            Graphics object consisting of 9 graphics primitives

        ::
        
            sage: R.<x> = QQ[]
            sage: f = Piecewise([[(-1,1),1/2+x-x^3]]) ## example 3
            sage: tf = f.trapezoid(6)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[a,0],[a,f(a)]],rgbcolor=(0.7,0.6,0.6)) for (a,b),f in tf.list()])
            sage: P+Q+L
            Graphics object consisting of 13 graphics primitives

        TESTS:

        Use variables other than x (:trac:`13836`)::

            sage: R.<y> = QQ[]
            sage: f1 = y^2
            sage: f2 = 5-y^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.trapezoid(4)
            Piecewise defined function with 4 parts, [[(0, 1/2), 1/2*y], [(1/2, 1), 9/2*y - 2], [(1, 3/2), 1/2*y + 2], [(3/2, 2), -7/2*y + 8]]

        """
        x = QQ[self.default_variable()].gen()
        def f(x0, x1):
            f0, f1 = self(x0), self(x1)
            return [[(x0,x1),f0+(f1-f0)*(x1-x0)**(-1)*(x-x0)]]
        rsum = self._riemann_sum_helper(N, f, initial=[])
        return Piecewise(rsum)

    def trapezoid_integral_approximation(self,N):
        """
        Returns the approximation given by the trapezoid rule for numerical
        integration based on a subdivision into N subintervals.
        
        EXAMPLES::
        
            sage: f1(x) = x^2                      ## example 1
            sage: f2(x) = 1-(1-x)^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: tf = f.trapezoid(6)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: ta = f.trapezoid_integral_approximation(6)
            sage: t = text('trapezoid approximation = %s'%ta, (1.5, 0.25))
            sage: a = f.integral(definite=True)
            sage: tt = text('area under curve = %s'%a, (1.5, -0.5))
            sage: P + Q + t + tt
            Graphics object consisting of 10 graphics primitives

        ::
        
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])  ## example 2
            sage: tf = f.trapezoid(4)
            sage: ta = f.trapezoid_integral_approximation(4)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: t = text('trapezoid approximation = %s'%ta, (1.5, 0.25))
            sage: a = f.integral(definite=True)
            sage: tt = text('area under curve = %s'%a, (1.5, -0.5))
            sage: P+Q+t+tt
            Graphics object consisting of 8 graphics primitives
        """
        def f(x0, x1):
            f0, f1 = self(x0), self(x1)
            return ((f1+f0)/2)*(x1-x0)
        return self._riemann_sum_helper(N, f)

    def critical_points(self):
        """
        Return the critical points of this piecewise function.
        
        .. warning::

           Uses maxima, which prints the warning to use results with
           caution. Only works for piecewise functions whose parts are
           polynomials with real critical not occurring on the
           interval endpoints.
        
        EXAMPLES::
        
            sage: R.<x> = QQ[]
            sage: f1 = x^0
            sage: f2 = 10*x - x^2
            sage: f3 = 3*x^4 - 156*x^3 + 3036*x^2 - 26208*x
            sage: f = Piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: expected = [5, 12, 13, 14]
            sage: all(abs(e-a) < 0.001 for e,a in zip(expected, f.critical_points()))
            True

        TESTS:

        Use variables other than x (:trac:`13836`)::

            sage: R.<y> = QQ[]
            sage: f1 = y^0
            sage: f2 = 10*y - y^2
            sage: f3 = 3*y^4 - 156*y^3 + 3036*y^2 - 26208*y
            sage: f = Piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
            sage: expected = [5, 12, 13, 14]
            sage: all(abs(e-a) < 0.001 for e,a in zip(expected, f.critical_points()))
            True
        """
        from sage.calculus.calculus import maxima
        x = QQ[self.default_variable()].gen()
        crit_pts = []
        for (a,b), f in self.list():
            for root in maxima.allroots(SR(f).diff(x)==0):
                root = float(root.rhs())
                if a < root < b:
                    crit_pts.append(root)
        return crit_pts
        
    def base_ring(self):
        """
        Returns the base ring of the function pieces.   This
        is useful when this class is extended.

        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = x^2-5
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: base_ring(f)
            Symbolic Ring

        ::

            sage: R.<x> = QQ[]
            sage: f1 = x^0
            sage: f2 = 10*x - x^2
            sage: f3 = 3*x^4 - 156*x^3 + 3036*x^2 - 26208*x
            sage: f = Piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
            sage: f.base_ring()
            Rational Field
        """
        return (self.functions()[0]).base_ring()

    def end_points(self):
        """
        Returns a list of all interval endpoints for this function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = x^2-5
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.end_points()
            [0, 1, 2, 3]
        """
        intervals = self.intervals()
        return [ intervals[0][0] ] + [b for a,b in intervals]

    def __call__(self,x0):
        """
        Evaluates self at x0. Returns the average value of the jump if x0
        is an interior endpoint of one of the intervals of self and the
        usual value otherwise.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f(0.5)
            1
            sage: f(5/2)
            e^(5/2)
            sage: f(5/2).n()
            12.1824939607035
            sage: f(1)
            1/2
        """
        #x0 = QQ(x0) ## does not allow for evaluation at pi
        n = self.length()
        endpts = self.end_points()
        for i in range(1,n):
            if x0 == endpts[i]:
                return (self.functions()[i-1](x0) + self.functions()[i](x0))/2
        if x0 == endpts[0]:
            return self.functions()[0](x0)
        if x0 == endpts[n]:
            return self.functions()[n-1](x0)
        for i in range(n):
            if endpts[i] < x0 < endpts[i+1]:
                return self.functions()[i](x0)
        raise ValueError("Value not defined outside of domain.")

    def which_function(self,x0):
        """
        Returns the function piece used to evaluate self at x0.
        
        EXAMPLES::
        
            sage: f1(z) = z
            sage: f2(x) = 1-x
            sage: f3(y) = exp(y)
            sage: f4(t) = sin(2*t)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.which_function(3/2)
            x |--> -x + 1
        """
        for (a,b), f in self.list():
            if a <= x0 <= b:
                return f
        raise ValueError("Function not defined outside of domain.")

    def default_variable(self):
        r"""
        Return the default variable. The default variable is defined as the
        first variable in the first piece that has a variable. If no pieces have
        a variable (each piece is a constant value), `x` is returned.

        The result is cached.

        AUTHOR: Paul Butler

        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 5*x
            sage: p = Piecewise([[(0,1),f1],[(1,4),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: p.default_variable()
            x

            sage: f1 = 3*var('y')
            sage: p = Piecewise([[(0,1),4],[(1,4),f1]])
            sage: p.default_variable()
            y

        """
        try:
            return self.__default_variable
        except AttributeError:
            pass
        for _, fun in self._list:
            try:
                fun = SR(fun)
                if fun.variables():
                    v = fun.variables()[0]
                    self.__default_variable = v
                    return v
            except TypeError:
                # pass if fun is lambda function
                pass
        # default to x
        v = var('x')
        self.__default_value = v
        return v

    def integral(self, x=None, a=None, b=None, definite=False):
        r"""
        By default, returns the indefinite integral of the function.
        If definite=True is given, returns the definite integral.

        AUTHOR: 

        - Paul Butler

        EXAMPLES::

            sage: f1(x) = 1-x
            sage: f = Piecewise([[(0,1),1],[(1,2),f1]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.integral(definite=True)
            1/2
        
        ::
        
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.integral(definite=True)
            1/2*pi
            
            sage: f1(x) = 2
            sage: f2(x) = 3 - x
            sage: f = Piecewise([[(-2, 0), f1], [(0, 3), f2]])
            sage: f.integral()
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 2 parts, [[(-2, 0), x |--> 2*x + 4], [(0, 3), x |--> -1/2*x^2 + 3*x + 4]]

            sage: f1(y) = -1
            sage: f2(y) = y + 3
            sage: f3(y) = -y - 1
            sage: f4(y) = y^2 - 1
            sage: f5(y) = 3
            sage: f = Piecewise([[(-4,-3),f1],[(-3,-2),f2],[(-2,0),f3],[(0,2),f4],[(2,3),f5]])
            sage: F = f.integral(y)
            sage: F
            Piecewise defined function with 5 parts, [[(-4, -3), y |--> -y - 4], [(-3, -2), y |--> 1/2*y^2 + 3*y + 7/2], [(-2, 0), y |--> -1/2*y^2 - y - 1/2], [(0, 2), y |--> 1/3*y^3 - y - 1/2], [(2, 3), y |--> 3*y - 35/6]]
            
        Ensure results are consistent with FTC::

            sage: F(-3) - F(-4)
            -1
            sage: F(-1) - F(-3)
            1
            sage: F(2) - F(0)
            2/3
            sage: f.integral(y, 0, 2)
            2/3
            sage: F(3) - F(-4)
            19/6
            sage: f.integral(y, -4, 3)
            19/6
            sage: f.integral(definite=True)
            19/6
            
        ::

            sage: f1(y) = (y+3)^2
            sage: f2(y) = y+3
            sage: f3(y) = 3
            sage: f = Piecewise([[(-infinity, -3), f1], [(-3, 0), f2], [(0, infinity), f3]])
            sage: f.integral()
            Piecewise defined function with 3 parts, [[(-Infinity, -3), y |--> 1/3*y^3 + 3*y^2 + 9*y + 9], [(-3, 0), y |--> 1/2*y^2 + 3*y + 9/2], [(0, +Infinity), y |--> 3*y + 9/2]]

        ::

            sage: f1(x) = e^(-abs(x))
            sage: f = Piecewise([[(-infinity, infinity), f1]])
            sage: f.integral(definite=True)
            2
            sage: f.integral()
            Piecewise defined function with 1 parts, [[(-Infinity, +Infinity), x |--> -1/2*((sgn(x) - 1)*e^(2*x) - 2*e^x*sgn(x) + sgn(x) + 1)*e^(-x) - 1]]

        ::

            sage: f = Piecewise([((0, 5), cos(x))])
            sage: f.integral()
            Piecewise defined function with 1 parts, [[(0, 5), x |--> sin(x)]]


        TESTS:

        Verify that piecewise integrals of zero work (:trac:`10841`)::

            sage: f0(x) = 0 
            sage: f = Piecewise([[(0,1),f0]])
            sage: f.integral(x,0,1)
            0
            sage: f = Piecewise([[(0,1), 0]])
            sage: f.integral(x,0,1)
            0
            sage: f = Piecewise([[(0,1), SR(0)]])
            sage: f.integral(x,0,1)
            0

        """
        if a != None and b != None:
            F = self.integral(x)
            return F(b) - F(a)

        if a != None or b != None:
            raise TypeError('only one endpoint given')

        area = 0 # cumulative definite integral of parts to the left of the current interval
        integrand_pieces = self.list()
        integrand_pieces.sort()
        new_pieces = []

        if x == None:
            x = self.default_variable()
        
        # The integral is computed by iterating over the pieces in order.
        # The definite integral for each piece is calculated and accumulated in `area`.
        # Thus at any time, `area` represents the definite integral of all the pieces
        # encountered so far. The indefinite integral of each piece is also calculated,
        # and the `area` before each piece is added to the piece.
        #
        # If a definite integral is requested, `area` is returned. 
        # Otherwise, a piecewise function is constructed from the indefinite integrals
        # and returned.
        #
        # An exception is made if integral is called on a piecewise function
        # that starts at -infinity. In this case, we do not try to calculate the
        # definite integral of the first piece, and the value of `area` remains 0
        # after the first piece.

        for (start, end), fun in integrand_pieces:
            fun = SR(fun)
            if start == -infinity and not definite:
                fun_integrated = fun.integral(x, end, x)
            else:
                try:
                    assume(start < x)
                except ValueError: # Assumption is redundant
                    pass
                fun_integrated = fun.integral(x, start, x) + area
                forget(start < x)
                if definite or end != infinity:
                    area += fun.integral(x, start, end)
            new_pieces.append([(start, end), SR(fun_integrated).function(x)])

        if definite:
            return SR(area)
        else:
            return Piecewise(new_pieces)

    def convolution(self, other):
        """
        Returns the convolution function,
        `f*g(t)=\int_{-\infty}^\infty f(u)g(t-u)du`, for compactly
        supported `f,g`.
        
        EXAMPLES::
        
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f = Piecewise([[(0,1),1*x^0]])  ## example 0
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g = f.convolution(f)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: h = f.convolution(g)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot(); Q = g.plot(rgbcolor=(1,1,0)); R = h.plot(rgbcolor=(0,1,1));
            sage: # Type show(P+Q+R) to view
            sage: f = Piecewise([[(0,1),1*x^0],[(1,2),2*x^0],[(2,3),1*x^0]])  ## example 1
            sage: g = f.convolution(f)
            sage: h = f.convolution(g)
            sage: P = f.plot(); Q = g.plot(rgbcolor=(1,1,0)); R = h.plot(rgbcolor=(0,1,1));
            sage: # Type show(P+Q+R) to view
            sage: f = Piecewise([[(-1,1),1]])                             ## example 2
            sage: g = Piecewise([[(0,3),x]])
            sage: f.convolution(g)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 3 parts, [[(-1, 1), 0], [(1, 2), -3/2*x], [(2, 4), -3/2*x]]
            sage: g = Piecewise([[(0,3),1*x^0],[(3,4),2*x^0]])
            sage: f.convolution(g)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 5 parts, [[(-1, 1), x + 1], [(1, 2), 3], [(2, 3), x], [(3, 4), -x + 8], [(4, 5), -2*x + 10]]
        """
        f = self
        g = other
        M = min(min(f.end_points()),min(g.end_points()))
        N = max(max(f.end_points()),max(g.end_points()))
        R2 = PolynomialRing(QQ,2,names=["tt","uu"])
        tt,uu = R2.gens()
        conv = 0
        f0 = f.functions()[0]
        g0 = g.functions()[0]
        R1 = f0.parent()
        xx = R1.gen()
        var = repr(xx)
        if len(f.intervals())==1 and len(g.intervals())==1:
            f = f.unextend()
            g = g.unextend()
            a1 = f.intervals()[0][0]
            a2 = f.intervals()[0][1]
            b1 = g.intervals()[0][0]
            b2 = g.intervals()[0][1]
            i1 = repr(f0).replace(var,repr(uu))
            i2 = repr(g0).replace(var,"("+repr(tt-uu)+")")
            cmd1 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, a1,    tt-b1)    ## if a1+b1 < tt < a2+b1
            cmd2 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, tt-b2, tt-b1)    ## if a1+b2 < tt < a2+b1
            cmd3 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, tt-b2, a2)       ## if a1+b2 < tt < a2+b2
            cmd4 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, a1, a2)          ## if a2+b1 < tt < a1+b2
            conv1 = maxima.eval(cmd1)
            conv2 = maxima.eval(cmd2)
            conv3 = maxima.eval(cmd3)
            conv4 = maxima.eval(cmd4)
            # this is a very, very, very ugly hack
            x = PolynomialRing(QQ,'x').gen()
            fg1 = sage_eval(conv1.replace("tt",var), {'x':x}) ## should be = R2(conv1)
            fg2 = sage_eval(conv2.replace("tt",var), {'x':x}) ## should be = R2(conv2)
            fg3 = sage_eval(conv3.replace("tt",var), {'x':x}) ## should be = R2(conv3)
            fg4 = sage_eval(conv4.replace("tt",var), {'x':x}) ## should be = R2(conv4)
            if a1-b1<a2-b2:
                if a2+b1!=a1+b2:
                    h = Piecewise([[(a1+b1,a1+b2),fg1],[(a1+b2,a2+b1),fg4],[(a2+b1,a2+b2),fg3]])
                else:
                    h = Piecewise([[(a1+b1,a1+b2),fg1],[(a1+b2,a2+b2),fg3]])
            else:
                if a1+b2!=a2+b1:
                    h = Piecewise([[(a1+b1,a2+b1),fg1],[(a2+b1,a1+b2),fg2],[(a1+b2,a2+b2),fg3]])
                else:
                    h = Piecewise([[(a1+b1,a2+b1),fg1],[(a2+b1,a2+b2),fg3]])
            return h
        
        if len(f.intervals())>1 or len(g.intervals())>1:
            z = Piecewise([[(-3*abs(N-M),3*abs(N-M)),0*xx**0]])
            ff = f.functions()
            gg = g.functions()
            intvlsf = f.intervals()
            intvlsg = g.intervals()
            for i in range(len(ff)):
                for j in range(len(gg)):
                    f0 = Piecewise([[intvlsf[i],ff[i]]])
                    g0 = Piecewise([[intvlsg[j],gg[j]]])
                    h = g0.convolution(f0)
                    z = z + h
            return z.unextend()
    
    def derivative(self):
        r"""
        Returns the derivative (as computed by maxima)
        Piecewise(I,`(d/dx)(self|_I)`), as I runs over the
        intervals belonging to self. self must be piecewise polynomial.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.derivative()
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 2 parts, [[(0, 1), x |--> 0], [(1, 2), x |--> -1]]
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.derivative()
            Piecewise defined function with 2 parts, [[(0, 1/2*pi), x |--> 0], [(1/2*pi, pi), x |--> 0]]
            
        ::
        
            sage: f = Piecewise([[(0,1), (x * 2)]], x)
            sage: f.derivative()
            Piecewise defined function with 1 parts, [[(0, 1), x |--> 2]]
        """
        x = self.default_variable()
        dlist = [[(a, b), derivative(f(x), x).function(x)] for (a,b),f in self.list()]
        return Piecewise(dlist)
 
    def tangent_line(self, pt):
        """
        Computes the linear function defining the tangent line of the
        piecewise function self.
        
        EXAMPLES::
        
            sage: f1(x) = x^2
            sage: f2(x) = 5-x^3+x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: tf = f.tangent_line(0.9) ## tangent line at x=0.9
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = tf.plot(rgbcolor=(0.7,0.2,0.2), plot_points=40)
            sage: P + Q
            Graphics object consisting of 4 graphics primitives
        """
        pt = QQ(pt)
        R = QQ[self.default_variable()]
        x = R.gen()
        der = self.derivative()
        tanline = (x-pt)*der(pt)+self(pt)
        dlist = [[(a, b), tanline] for (a,b),f in self.list()]
        return Piecewise(dlist)
        
    def plot(self, *args, **kwds):
        """
        Returns the plot of self.
        
        Keyword arguments are passed onto the plot command for each piece
        of the function. E.g., the plot_points keyword affects each
        segment of the plot.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot(rgbcolor=(0.7,0.1,0), plot_points=40)
            sage: P
            Graphics object consisting of 4 graphics primitives

        Remember: to view this, type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.

        TESTS:

        We should not add each piece to the legend individually, since
        this creates duplicates (:trac:`12651`). This tests that only
        one of the graphics objects in the plot has a non-``None``
        ``legend_label``::

            sage: f1 = sin(x)
            sage: f2 = cos(x)
            sage: f = piecewise([[(-1,0), f1],[(0,1), f2]])
            sage: p = f.plot(legend_label='$f(x)$')
            sage: lines = [
            ...     line
            ...     for line in p._objects
            ...     if line.options()['legend_label'] is not None ]
            sage: len(lines)
            1
        """
        from sage.plot.all import plot, Graphics

        g = Graphics()

        for i, ((a,b), f) in enumerate(self.list()):
            # If it's the first piece, pass all arguments. Otherwise,
            # filter out 'legend_label' so that we don't add each
            # piece to the legend separately (trac #12651).
            if i != 0 and 'legend_label' in kwds:
                del kwds['legend_label']

            g += plot(f, a, b, *args, **kwds)

        return g

    def fourier_series_cosine_coefficient(self,n,L):
        r"""
        Returns the n-th Fourier series coefficient of
        `\cos(n\pi x/L)`, `a_n`.
        
        INPUT:
        
        
        -  ``self`` - the function f(x), defined over -L x L
        
        -  ``n`` - an integer n=0
        
        -  ``L`` - (the period)/2
        
        
        OUTPUT:
        `a_n = \frac{1}{L}\int_{-L}^L f(x)\cos(n\pi x/L)dx`
        
        EXAMPLES::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_cosine_coefficient(2,1)
            pi^(-2)
            sage: f(x) = x^2
            sage: f = Piecewise([[(-pi,pi),f]])
            sage: f.fourier_series_cosine_coefficient(2,pi)
            1
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_cosine_coefficient(5,pi)
            -3/5/pi
        """
        from sage.all import cos, pi
        x = var('x')
        result = sum([(f(x)*cos(pi*x*n/L)/L).integrate(x, a, b)
                      for (a,b), f in self.list()])
        if is_Expression(result):
            return result.simplify_trig()
        return result
    
    def fourier_series_sine_coefficient(self,n,L):
        r"""
        Returns the n-th Fourier series coefficient of
        `\sin(n\pi x/L)`, `b_n`.
        
        INPUT:
        
        
        -  ``self`` - the function f(x), defined over -L x L
        
        -  ``n`` - an integer n0
        
        -  ``L`` - (the period)/2
        
        
        OUTPUT:
        `b_n = \frac{1}{L}\int_{-L}^L f(x)\sin(n\pi x/L)dx`
        
        EXAMPLES::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_sine_coefficient(2,1)  # L=1, n=2
            0
        """
        from sage.all import sin, pi
        x = var('x')
        result = sum([(f(x)*sin(pi*x*n/L)/L).integrate(x, a, b)
                      for (a,b), f in self.list()])
        if is_Expression(result):
            return result.simplify_trig()
        return result

    def _fourier_series_helper(self, N, L, scale_function):
        r"""
        A helper function for the construction of Fourier series. The
        argument scale_function is a function which takes in n,
        representing the `n^{th}` coefficient, and return an
        expression to scale the sine and cosine coefficients by.
        
        EXAMPLES::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f._fourier_series_helper(3, 1, lambda n: 1)
            cos(2*pi*x)/pi^2 - 4*cos(pi*x)/pi^2 + 1/3
        """
        from sage.all import pi, sin, cos, srange
        x = self.default_variable()
        a0 = self.fourier_series_cosine_coefficient(0,L)
        result = a0/2 + sum([(self.fourier_series_cosine_coefficient(n,L)*cos(n*pi*x/L) +
                              self.fourier_series_sine_coefficient(n,L)*sin(n*pi*x/L))*
                             scale_function(n)
                             for n in srange(1,N)])
        return result.expand()
        

    def fourier_series_partial_sum(self,N,L):
        r"""
        Returns the partial sum
        
        .. math::
        
           f(x) \sim \frac{a_0}{2} + \sum_{n=1}^N [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        as a string.
        
        EXAMPLE::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_partial_sum(3,1)
            cos(2*pi*x)/pi^2 - 4*cos(pi*x)/pi^2 + 1/3
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum(3,pi)
            -3*cos(x)/pi - 3*sin(2*x)/pi + 3*sin(x)/pi - 1/4
        """
        return self._fourier_series_helper(N, L, lambda n: 1)
  
    def fourier_series_partial_sum_cesaro(self,N,L):
        r"""
        Returns the Cesaro partial sum
        
        .. math::
        
           f(x) \sim \frac{a_0}{2} + \sum_{n=1}^N (1-n/N)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        
        as a string. This is a "smoother" partial sum - the Gibbs
        phenomenon is mollified.
        
        EXAMPLE::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_partial_sum_cesaro(3,1)
            1/3*cos(2*pi*x)/pi^2 - 8/3*cos(pi*x)/pi^2 + 1/3
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum_cesaro(3,pi)
            -2*cos(x)/pi - sin(2*x)/pi + 2*sin(x)/pi - 1/4
        """
        return self._fourier_series_helper(N, L, lambda n: 1-n/N)

    def fourier_series_partial_sum_hann(self,N,L):
        r"""
        Returns the Hann-filtered partial sum (named after von Hann, not
        Hamming)
        
        .. math::
        
           f(x) \sim \frac{a_0}{2} + \sum_{n=1}^N H_N(n)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        as a string, where `H_N(x) = (1+\cos(\pi x/N))/2`. This is
        a "smoother" partial sum - the Gibbs phenomenon is mollified.
        
        EXAMPLE::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_partial_sum_hann(3,1)
            1/4*cos(2*pi*x)/pi^2 - 3*cos(pi*x)/pi^2 + 1/3
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum_hann(3,pi)
            -9/4*cos(x)/pi - 3/4*sin(2*x)/pi + 9/4*sin(x)/pi - 1/4
        """
        from sage.all import cos, pi
        return self._fourier_series_helper(N, L, lambda n: (1+cos(pi*n/N))/2)

    def fourier_series_partial_sum_filtered(self,N,L,F):
        r"""
        Returns the "filtered" partial sum
        
        .. math::
        
           f(x) \sim \frac{a_0}{2} + \sum_{n=1}^N F_n*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        as a string, where `F = [F_1,F_2, ..., F_{N}]` is a list
        of length `N` consisting of real numbers. This can be used
        to plot FS solutions to the heat and wave PDEs.
        
        EXAMPLE::
        
            sage: f(x) = x^2
            sage: f = Piecewise([[(-1,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_partial_sum_filtered(3,1,[1,1,1])
            cos(2*pi*x)/pi^2 - 4*cos(pi*x)/pi^2 + 1/3
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum_filtered(3,pi,[1,1,1])
            -3*cos(x)/pi - 3*sin(2*x)/pi + 3*sin(x)/pi - 1/4
        """
        return self._fourier_series_helper(N, L, lambda n: F[n])
        
    def plot_fourier_series_partial_sum(self,N,L,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        
        .. math::
        
           f(x) \sim \frac{a_0}{2} +  sum_{n=1}^N [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        over xmin x xmin.
        
        EXAMPLE::
        
            sage: f1(x) = -2
            sage: f2(x) = 1
            sage: f3(x) = -1
            sage: f4(x) = 2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot_fourier_series_partial_sum(3,pi,-5,5)    # long time
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum(15,pi,-5,5)   # long time
        
        Remember, to view this type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.
        """
        from sage.plot.all import plot
        return plot(self.fourier_series_partial_sum(N,L), xmin, xmax, **kwds)

    def plot_fourier_series_partial_sum_cesaro(self,N,L,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        
        .. math::
        
                     f(x) \sim \frac{a_0}{2} +                     \sum_{n=1}^N (1-n/N)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        
        over xmin x xmin. This is a "smoother" partial sum - the Gibbs
        phenomenon is mollified.
        
        EXAMPLE::
        
            sage: f1(x) = -2
            sage: f2(x) = 1
            sage: f3(x) = -1
            sage: f4(x) = 2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot_fourier_series_partial_sum_cesaro(3,pi,-5,5)    # long time
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum_cesaro(15,pi,-5,5)   # long time
        
        Remember, to view this type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.
        """     
        from sage.plot.all import plot
        return plot(self.fourier_series_partial_sum_cesaro(N,L), xmin, xmax, **kwds)
    
    def plot_fourier_series_partial_sum_hann(self,N,L,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        
        .. math::
        
           f(x) \sim \frac{a_0}{2} + \sum_{n=1}^N H_N(n)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        
        over xmin x xmin, where H_N(x) = (0.5)+(0.5)\*cos(x\*pi/N) is the
        N-th Hann filter.
        
        EXAMPLE::
        
            sage: f1(x) = -2
            sage: f2(x) = 1
            sage: f3(x) = -1
            sage: f4(x) = 2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot_fourier_series_partial_sum_hann(3,pi,-5,5)    # long time
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum_hann(15,pi,-5,5)   # long time
        
        Remember, to view this type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.
        """     
        from sage.plot.all import plot
        return plot(self.fourier_series_partial_sum_hann(N,L), xmin, xmax, **kwds)
        
    def plot_fourier_series_partial_sum_filtered(self,N,L,F,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        
        .. math::
        
                     f(x) \sim \frac{a_0}{2} +                     \sum_{n=1}^N F_n*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         
        
        
        over xmin x xmin, where `F = [F_1,F_2, ..., F_{N}]` is a
        list of length `N` consisting of real numbers. This can be
        used to plot FS solutions to the heat and wave PDEs.
        
        EXAMPLE::
        
            sage: f1(x) = -2
            sage: f2(x) = 1
            sage: f3(x) = -1
            sage: f4(x) = 2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot_fourier_series_partial_sum_filtered(3,pi,[1]*3,-5,5)    # long time
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum_filtered(15,pi,[1]*15,-5,5)   # long time
        
        Remember, to view this type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.
        """
        from sage.plot.all import plot
        return plot(self.fourier_series_partial_sum_filtered(N,L,F), xmin, xmax, **kwds)
                
    def fourier_series_value(self,x,L):
        r"""
        Returns the value of the Fourier series coefficient of self at
        `x`,
        
        
        .. math::
        
                     f(x) \sim \frac{a_0}{2} +                     \sum_{n=1}^\infty [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],         \ \ \ -L<x<L.         
        
        
        This method applies to piecewise non-polynomial functions as well.
        
        INPUT:
        
        
        -  ``self`` - the function f(x), defined over -L x L
        
        -  ``x`` - a real number
        
        -  ``L`` - (the period)/2
        
        
        OUTPUT: `(f^*(x+)+f^*(x-)/2`, where `f^*` denotes
        the function `f` extended to `\RR` with period
        `2L` (Dirichlet's Theorem for Fourier series).
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(-10,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.fourier_series_value(101,10)  
            1/2
            sage: f.fourier_series_value(100,10)
            1
            sage: f.fourier_series_value(10,10)
            1/2*sin(20) + 1/2
            sage: f.fourier_series_value(20,10)
            1
            sage: f.fourier_series_value(30,10)
            1/2*sin(20) + 1/2
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(-pi,0),lambda x:0],[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_value(-1,pi)
            0
            sage: f.fourier_series_value(20,pi)
            -1
            sage: f.fourier_series_value(pi/2,pi)
            1/2
        """
        xnew = x - int(RR(x/(2*L)))*2*L
        endpts = self.end_points()
        if xnew == endpts[0] or xnew == endpts[-1]:
            return (self.functions()[0](endpts[0]) + self.functions()[-1](endpts[-1]))/2
        else:
            return self(xnew)

    def cosine_series_coefficient(self,n,L):
        r"""
        Returns the n-th cosine series coefficient of
        `\cos(n\pi x/L)`, `a_n`.
        
        INPUT:
        
        
        -  ``self`` - the function f(x), defined over 0 x L (no
           checking is done to insure this)
        
        -  ``n`` - an integer n=0
        
        -  ``L`` - (the period)/2
        
        
        OUTPUT:
        `a_n = \frac{2}{L}\int_{-L}^L f(x)\cos(n\pi x/L)dx` such
        that
        
        .. math::
        
                     f(x) \sim \frac{a_0}{2} +                     \sum_{n=1}^\infty a_n\cos(\frac{n\pi x}{L}),\ \ 0<x<L.         
        
        
        
        EXAMPLES::
        
            sage: f(x) = x
            sage: f = Piecewise([[(0,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.cosine_series_coefficient(2,1)  
            0
            sage: f.cosine_series_coefficient(3,1)
            -4/9/pi^2
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.cosine_series_coefficient(2,pi)
            0
            sage: f.cosine_series_coefficient(3,pi)
            2/pi
            sage: f.cosine_series_coefficient(111,pi)
            2/37/pi
            sage: f1 = lambda x: x*(pi-x)
            sage: f = Piecewise([[(0,pi),f1]])
            sage: f.cosine_series_coefficient(0,pi)
            1/3*pi^2

        """
        from sage.all import cos, pi
        x = var('x')
        result = sum([(2*f(x)*cos(pi*x*n/L)/L).integrate(x, a, b)
                      for (a,b), f in self.list()])
        if is_Expression(result):
            return result.simplify_trig()
        return result


    def sine_series_coefficient(self,n,L):
        r"""
        Returns the n-th sine series coefficient of
        `\sin(n\pi x/L)`, `b_n`.
        
        INPUT:
        
        -  ``self`` - the function f(x), defined over 0 x L (no
           checking is done to insure this)
        
        -  ``n`` - an integer n0
        
        -  ``L`` - (the period)/2
        
        OUTPUT:

        `b_n = \frac{2}{L}\int_{-L}^L f(x)\sin(n\pi x/L)dx` such
        that
        
        .. math::
        
           f(x) \sim \sum_{n=1}^\infty b_n\sin(\frac{n\pi x}{L}),\ \ 0<x<L.         
        
        EXAMPLES::
        
            sage: f(x) = 1
            sage: f = Piecewise([[(0,1),f]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.sine_series_coefficient(2,1)  
            0
            sage: f.sine_series_coefficient(3,1)
            4/3/pi
        """
        from sage.all import sin, pi
        x = var('x')
        result = sum([(2*f(x)*sin(pi*x*n/L)/L).integrate(x, a, b)
                      for (a,b), f in self.list()])
        if is_Expression(result):
            return result.simplify_trig()
        return result

    def laplace(self, x='x', s='t'):
        r"""
        Returns the Laplace transform of self with respect to the variable
        var.
        
        INPUT:
        
        
        -  ``x`` - variable of self
        
        -  ``s`` - variable of Laplace transform.
        
        
        We assume that a piecewise function is 0 outside of its domain and
        that the left-most endpoint of the domain is 0.
        
        EXAMPLES::
        
            sage: x, s, w = var('x, s, w')
            sage: f = Piecewise([[(0,1),1],[(1,2), 1-x]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.laplace(x, s)
            -e^(-s)/s + (s + 1)*e^(-2*s)/s^2 + 1/s - e^(-s)/s^2
            sage: f.laplace(x, w)
            -e^(-w)/w + (w + 1)*e^(-2*w)/w^2 + 1/w - e^(-w)/w^2
            
        ::
        
            sage: y, t = var('y, t')
            sage: f = Piecewise([[(1,2), 1-y]]) 
            sage: f.laplace(y, t)
            (t + 1)*e^(-2*t)/t^2 - e^(-t)/t^2
            
        ::
        
            sage: s = var('s')
            sage: t = var('t')
            sage: f1(t) = -t
            sage: f2(t) = 2
            sage: f = Piecewise([[(0,1),f1],[(1,infinity),f2]])
            sage: f.laplace(t,s)
            (s + 1)*e^(-s)/s^2 + 2*e^(-s)/s - 1/s^2
        """
        from sage.all import assume, exp, forget
        x = var(x)
        s = var(s)
        assume(s>0)
        result =  sum([(SR(f)*exp(-s*x)).integral(x,a,b)
                       for (a,b),f in self.list()])
        forget(s>0)
        return result

    def _make_compatible(self, other):
        """
        Returns self and other extended to be defined on the same domain as
        well as a refinement of their intervals. This is used for adding
        and multiplying piecewise functions.
        
        EXAMPLES::
        
            sage: R.<x> = QQ[]
            sage: f1 = Piecewise([[(0, 2), x]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f2 = Piecewise([[(1, 3), x^2]])
            sage: f1._make_compatible(f2)
            (Piecewise defined function with 2 parts, [[(0, 2), x], [(2, 3), 0]],
            Piecewise defined function with 2 parts, [[(0, 1), 0], [(1, 3), x^2]],
            [(0, 1), (1, 2), (2, 3)])
        """
        a1, b1 = self.domain()
        a2, b2 = other.domain()
        a = min(a1, a2)
        b = max(b1, b2)
        F = self.extend_by_zero_to(a,b)
        G = other.extend_by_zero_to(a,b)
        endpts = list(set(F.end_points()).union(set(G.end_points())))
        endpts.sort()
        return F, G, zip(endpts, endpts[1:])
   
    def __add__(self,other):
        """
        Returns the piecewise defined function which is the sum of self and
        other. Does not require both domains be the same.
        
        EXAMPLES::
        
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g1 = x-2
            sage: g2 = x-5
            sage: g = Piecewise([[(0,5),g1],[(5,10),g2]])
            sage: h = f+g
            sage: h
            Piecewise defined function with 5 parts, [[(0, 1), x - 1], [(1, 2), -1], [(2, 3), 3*x - 2], [(3, 5), 8], [(5, 10), 5]]
        
        Note that in this case the functions must be defined using
        polynomial expressions *not* using the lambda notation.
        """
        F, G, intervals = self._make_compatible(other)
        fcn = []
        for a,b in intervals:
            fcn.append([(a,b), F.which_function(b)+G.which_function(b)])        
        return Piecewise(fcn)
        
    def __mul__(self,other):
        r"""
        Returns the piecewise defined function which is the product of one
        piecewise function (self) with another one (other).
        
        EXAMPLES::
        
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g1 = x-2
            sage: g2 = x-5
            sage: g = Piecewise([[(0,5),g1],[(5,10),g2]])
            sage: h = f*g
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: h
            Piecewise defined function with 5 parts, [[(0, 1), x - 2], [(1, 2), -x^2 + 3*x - 2], [(2, 3), 2*x^2 - 4*x], [(3, 5), -x^2 + 12*x - 20], [(5, 10), -x^2 + 15*x - 50]]
            sage: g*(11/2)
            Piecewise defined function with 2 parts, [[(0, 5), 11/2*x - 11], [(5, 10), 11/2*x - 55/2]]
            
        Note that in this method the functions must be defined using
        polynomial expressions *not* using the lambda notation.
        """
        ## needed for scalar multiplication
        if isinstance(other,Rational) or isinstance(other,Integer):
            return Piecewise([[(a,b), other*f] for (a,b),f in self.list()])
        else:
            F, G, intervals = self._make_compatible(other)
            fcn = []
            for a,b in intervals:
                fcn.append([(a,b),F.which_function(b)*G.which_function(b)])     
            return Piecewise(fcn)

    __rmul__ = __mul__

    def __eq__(self,other):
        """
        Implements Boolean == operator.

        EXAMPLES::
        
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g = Piecewise([[(0,1),1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f==g
            True
        """
        return self.list()==other.list()
