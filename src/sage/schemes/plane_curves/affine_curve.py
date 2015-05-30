"""
Affine plane curves over a general ring

AUTHORS:

- William Stein (2005-11-13)

- David Joyner (2005-11-13)

- David Kohel (2006-01)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.all import singular

from sage.misc.all import add

from sage.rings.all import degree_lowest_rational_function

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_affine


from curve import Curve_generic

class AffineSpaceCurve_generic(Curve_generic, AlgebraicScheme_subscheme_affine):
    def _repr_type(self):
        return "Affine Space"

    def __init__(self, A, X):
        if not is_AffineSpace(A):
            raise TypeError("A (=%s) must be an affine space"%A)
        Curve_generic.__init__(self, A, X)
        d = self.dimension()
        if d != 1:
            raise ValueError("defining equations (=%s) define a scheme of dimension %s != 1"%(X,d))

class AffineCurve_generic(Curve_generic):
    def __init__(self, A, f):
        P = f.parent()
        if not (is_AffineSpace(A) and A.dimension != 2):
            raise TypeError("Argument A (= %s) must be an affine plane."%A)
        Curve_generic.__init__(self, A, [f])

    def _repr_type(self):
        return "Affine"

    def divisor_of_function(self, r):
        """
        Return the divisor of a function on a curve.

        INPUT: r is a rational function on X

        OUTPUT:


        -  ``list`` - The divisor of r represented as a list of
           coefficients and points. (TODO: This will change to a more
           structural output in the future.)


        EXAMPLES::

            sage: F = GF(5)
            sage: P2 = AffineSpace(2, F, names = 'xy')
            sage: R = P2.coordinate_ring()
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: K = FractionField(R)
            sage: r = 1/x
            sage: C.divisor_of_function(r)     # todo: not implemented (broken)
                  [[-1, (0, 0, 1)]]
            sage: r = 1/x^3
            sage: C.divisor_of_function(r)     # todo: not implemented (broken)
                  [[-3, (0, 0, 1)]]
        """
        F = self.base_ring()
        f = self.defining_polynomial()
        pts = self.places_on_curve()
        numpts = len(pts)
        R = f.parent()
        x,y = R.gens()
        R0 = PolynomialRing(F,3,names = [str(x),str(y),"t"])
        vars0 = R0.gens()
        t = vars0[2]
        divf = []
        for pt0 in pts:
            if pt0[2] != F(0):
                lcs = self.local_coordinates(pt0,5)
                yt = lcs[1]
                xt = lcs[0]
                ldg = degree_lowest_rational_function(r(xt,yt),t)
                if ldg[0] != 0:
                    divf.append([ldg[0],pt0])
        return divf

    def local_coordinates(self, pt, n):
        r"""
        Return local coordinates to precision n at the given point.

            Behaviour is flaky - some choices of `n` are worst that
            others.


        INPUT:


        -  ``pt`` - an F-rational point on X which is not a
           point of ramification for the projection (x,y) - x.

        -  ``n`` - the number of terms desired


        OUTPUT: x = x0 + t y = y0 + power series in t

        EXAMPLES::

            sage: F = GF(5)
            sage: pt = (2,3)
            sage: R = PolynomialRing(F,2, names = ['x','y'])
            sage: x,y = R.gens()
            sage: f = y^2-x^9-x
            sage: C = Curve(f)
            sage: C.local_coordinates(pt, 9)
            [t + 2, -2*t^12 - 2*t^11 + 2*t^9 + t^8 - 2*t^7 - 2*t^6 - 2*t^4 + t^3 - 2*t^2 - 2]
        """
        f = self.defining_polynomial()
        R = f.parent()
        F = self.base_ring()
        p = F.characteristic()
        x0 = F(pt[0])
        y0 = F(pt[1])
        astr = ["a"+str(i) for i in range(1,2*n)]
        x,y = R.gens()
        R0 = PolynomialRing(F,2*n+2,names = [str(x),str(y),"t"]+astr)
        vars0 = R0.gens()
        t = vars0[2]
        yt = y0*t**0+add([vars0[i]*t**(i-2) for i in range(3,2*n+2)])
        xt = x0+t
        ft = f(xt,yt)
        S = singular
        S.eval('ring s = '+str(p)+','+str(R0.gens())+',lp;')
        S.eval('poly f = '+str(ft) + ';')
        c = S('coeffs(%s, t)'%ft)
        N = int(c.size())
        b = ["%s[%s,1],"%(c.name(), i) for i in range(2,N//2-4)]
        b = ''.join(b)
        b = b[:len(b)-1] # to cut off the trailing comma
        cmd = 'ideal I = '+b
        S.eval(cmd)
        S.eval('short=0')    # print using *'s and ^'s.
        c = S.eval('slimgb(I)')
        d = c.split("=")
        d = d[1:]
        d[len(d)-1] += "\n"
        e = [x[:x.index("\n")] for x in d]
        vals = []
        for x in e:
            for y in vars0:
                if str(y) in x:
                    if len(x.replace(str(y),"")) != 0:
                        i = x.find("-")
                        if i>0:
                            vals.append([eval(x[1:i]),x[:i],F(eval(x[i+1:]))])
                        i = x.find("+")
                        if i>0:
                            vals.append([eval(x[1:i]),x[:i],-F(eval(x[i+1:]))])
                    else:
                        vals.append([eval(str(y)[1:]),str(y),F(0)])
        vals.sort()
        k = len(vals)
        v = [x0+t,y0+add([vals[i][2]*t**(i+1) for i in range(k)])]
        return v

    def plot(self, *args, **kwds):
        """
        Plot the real points on this affine plane curve.

        INPUT:


        -  ``self`` - an affine plane curve

        -  ``*args`` - optional tuples (variable, minimum, maximum) for
           plotting dimensions

        -  ``**kwds`` - optional keyword arguments passed on to
           ``implicit_plot``


        EXAMPLES:

        A cuspidal curve::

            sage: R.<x, y> = QQ[]
            sage: C = Curve(x^3 - y^2)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive

        A 5-nodal curve of degree 11.  This example also illustrates
        some of the optional arguments::

            sage: R.<x, y> = ZZ[]
            sage: C = Curve(32*x^2 - 2097152*y^11 + 1441792*y^9 - 360448*y^7 + 39424*y^5 - 1760*y^3 + 22*y - 1)
            sage: C.plot((x, -1, 1), (y, -1, 1), plot_points=400)
            Graphics object consisting of 1 graphics primitive

        A line over `\mathbf{RR}`::

            sage: R.<x, y> = RR[]
            sage: C = Curve(R(y - sqrt(2)*x))
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
        """
        I = self.defining_ideal()
        return I.plot(*args, **kwds)

class AffineCurve_finite_field(AffineCurve_generic):
    def rational_points(self, algorithm="enum"):
        r"""
        Return sorted list of all rational points on this curve.

        Use *very* naive point enumeration to find all rational points on
        this curve over a finite field.

        EXAMPLE::

            sage: A.<x,y> = AffineSpace(2,GF(9,'a'))
            sage: C = Curve(x^2 + y^2 - 1)
            sage: C
            Affine Curve over Finite Field in a of size 3^2 defined by x^2 + y^2 - 1
            sage: C.rational_points()
            [(0, 1), (0, 2), (1, 0), (2, 0), (a + 1, a + 1), (a + 1, 2*a + 2), (2*a + 2, a + 1), (2*a + 2, 2*a + 2)]
        """
        f = self.defining_polynomial()
        R = f.parent()
        K = R.base_ring()
        points = []
        for x in K:
            for y in K:
                if f(x,y) == 0:
                    points.append(self((x,y)))
        points.sort()
        return points


class AffineCurve_prime_finite_field(AffineCurve_finite_field):
    # CHECK WHAT ASSUMPTIONS ARE MADE REGARDING AFFINE VS. PROJECTIVE MODELS!!!
    # THIS IS VERY DIRTY STILL -- NO DATASTRUCTURES FOR DIVISORS.

    def riemann_roch_basis(self,D):
        """
        Interfaces with Singular's BrillNoether command.

        INPUT:


        -  ``self`` - a plane curve defined by a polynomial eqn f(x,y)
           = 0 over a prime finite field F = GF(p) in 2 variables x,y
           representing a curve X: f(x,y) = 0 having n F-rational
           points (see the Sage function places_on_curve)

        -  ``D`` - an n-tuple of integers
           `(d1, ..., dn)` representing the divisor
           `Div = d1*P1+...+dn*Pn`, where
           `X(F) = \{P1,...,Pn\}`.
           *The ordering is that dictated by places_on_curve.*


        OUTPUT: basis of L(Div)

        EXAMPLE::

            sage: R = PolynomialRing(GF(5),2,names = ["x","y"])
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: D = [6,0,0,0,0,0]
            sage: C.riemann_roch_basis(D)
            [1, (y^2*z^4 - x*z^5)/x^6, (y^2*z^5 - x*z^6)/x^7, (y^2*z^6 - x*z^7)/x^8]
        """
        f = self.defining_polynomial()
        R = f.parent()
        F = self.base_ring()
        p = F.characteristic()
        Dstr = str(tuple(D))
        G = singular(','.join([str(x) for x in D]), type='intvec')

        singular.LIB('brnoeth.lib')

        S = singular.ring(p, R.gens(), 'lp')
        fsing = singular(str(f))

        X = fsing.Adj_div()
        P = singular.NSplaces(1, X)
        T = P[1][2]
        T.set_ring()

        LG = G.BrillNoether(P)

        dim = len(LG)
        basis = [(LG[i][1], LG[i][2]) for i in range(1,dim+1)]
        x, y, z = PolynomialRing(F, 3, names = ["x","y","z"]).gens()
        V = []
        for g in basis:
            T.set_ring()  # necessary...
            V.append(eval(g[0].sage_polystring())/eval(g[1].sage_polystring()))
        return V


    def rational_points(self, algorithm="enum"):
        r"""
        Return sorted list of all rational points on this curve.

        INPUT:


        -  ``algorithm`` - string:

           +  ``'enum'`` - straightforward enumeration

           +  ``'bn'`` - via Singular's Brill-Noether package.

           +  ``'all'`` - use all implemented algorithms and
              verify that they give the same answer, then return it


        .. note::

           The Brill-Noether package does not always work. When it
           fails a RuntimeError exception is raised.

        EXAMPLE::

            sage: x, y = (GF(5)['x,y']).gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f); C
            Affine Curve over Finite Field of size 5 defined by -x^9 + y^2 - x
            sage: C.rational_points(algorithm='bn')
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: C = Curve(x - y + 1)
            sage: C.rational_points()
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]

        We compare Brill-Noether and enumeration::

            sage: x, y = (GF(17)['x,y']).gens()
            sage: C = Curve(x^2 + y^5 + x*y - 19)
            sage: v = C.rational_points(algorithm='bn')
            sage: w = C.rational_points(algorithm='enum')
            sage: len(v)
            20
            sage: v == w
            True
        """
        if algorithm == "enum":

            return AffineCurve_finite_field.rational_points(self, algorithm="enum")

        elif algorithm == "bn":
            f = self.defining_polynomial()._singular_()
            singular = f.parent()
            singular.lib('brnoeth')
            try:
                X1 = f.Adj_div()
            except (TypeError, RuntimeError) as s:
                raise RuntimeError(str(s) + "\n\n ** Unable to use the Brill-Noether Singular package to compute all points (see above).")

            X2 = singular.NSplaces(1, X1)
            R = X2[5][1][1]
            singular.set_ring(R)

            # We use sage_flattened_str_list since iterating through
            # the entire list through the sage/singular interface directly
            # would involve hundreds of calls to singular, and timing issues
            # with the expect interface could crop up.  Also, this is vastly
            # faster (and more robust).
            v = singular('POINTS').sage_flattened_str_list()
            pnts = [self(int(v[3*i]), int(v[3*i+1])) for i in range(len(v)//3) if int(v[3*i+2])!=0]
            # remove multiple points
            pnts = sorted(set(pnts))
            return pnts

        elif algorithm == "all":

            S_enum = self.rational_points(algorithm = "enum")
            S_bn = self.rational_points(algorithm = "bn")
            if S_enum != S_bn:
                raise RuntimeError("Bug in rational_points -- different algorithms give different answers for curve %s!"%self)
            return S_enum

        else:
            raise ValueError("No algorithm '%s' known"%algorithm)
