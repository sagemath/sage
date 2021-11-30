r"""
Projective plane conics over a field

AUTHORS:

- Marco Streng (2010-07-20)

- Nick Alexander (2008-01-08)

"""
#*****************************************************************************
#       Copyright (C) 2008 Nick Alexander <ncalexander@gmail.com>
#       Copyright (C) 2009/2010 Marco Streng <marco.streng@gmail.com>
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

from sage.rings.all import PolynomialRing

import sage.rings.abc

from sage.modules.free_module_element import vector
from sage.structure.sequence import Sequence
from sage.structure.element import is_Vector
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.matrix.constructor import Matrix
from sage.structure.element import is_Matrix

from sage.schemes.curves.projective_curve import ProjectivePlaneCurve

from sage.categories.fields import Fields
_Fields = Fields()

class ProjectiveConic_field(ProjectivePlaneCurve):
    r"""
    Create a projective plane conic curve over a field.
    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - Z^2)
        Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Rational Field defined by X^2 + Y^2 - Z^2

    TESTS::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: Conic([K(1), 1, -1])._test_pickling()
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES:

        ::

            sage: c = Conic([1, 1, 1]); c
            Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2
        """
        ProjectivePlaneCurve.__init__(self, A, f)
        self._coefficients = [f[(2,0,0)], f[(1,1,0)], f[(1,0,1)],
                                f[(0,2,0)], f[(0,1,1)], f[(0,0,2)]]
        self._parametrization = None
        self._diagonal_matrix = None

        self._rational_point = None




    def _repr_type(self):
        r"""
        Returns ``'Projective Conic'``, which is the first part of the
        plain text representation of this object as output by
        the function ``_repr_`` of the class ``Curve_generic``.

        EXAMPLES::

            sage: c = Conic([1, 1, 1]); c
            Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2
            sage: c._repr_()
            'Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2'
            sage: c._repr_type()
            'Projective Conic'
        """
        return "Projective Conic"

    def base_extend(self, S):
        r"""
        Returns the conic over ``S`` given by the same equation as ``self``.

        EXAMPLES::

            sage: c = Conic([1, 1, 1]); c
            Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2
            sage: c.has_rational_point()
            False
            sage: d = c.base_extend(QuadraticField(-1, 'i')); d
            Projective Conic Curve over Number Field in i with defining polynomial x^2 + 1 with i = 1*I defined by x^2 + y^2 + z^2
            sage: d.rational_point(algorithm = 'rnfisnorm')
            (i : 1 : 0)
        """
        if S in _Fields:
            B = self.base_ring()
            if B == S:
                return self
            if not S.has_coerce_map_from(B):
                raise ValueError("No natural map from the base ring of self " \
                                  "(= %s) to S (= %s)" % (self, S))
            from .constructor import Conic
            con = Conic([S(c) for c in self.coefficients()], \
                        self.variable_names())
            if self._rational_point is not None:
                pt = [S(c) for c in Sequence(self._rational_point)]
                if not pt == [0,0,0]:
                    # The following line stores the point in the cache
                    # if (and only if) there is no point in the cache.
                    pt = con.point(pt)
            return con
        return ProjectivePlaneCurve.base_extend(self, S)

    def cache_point(self, p):
        r"""
        Replace the point in the cache of ``self`` by ``p`` for use
        by ``self.rational_point()`` and ``self.parametrization()``.

        EXAMPLES::

            sage: c = Conic([1, -1, 1])
            sage: c.point([15, 17, 8])
            (15/8 : 17/8 : 1)
            sage: c.rational_point()
            (15/8 : 17/8 : 1)
            sage: c.cache_point(c.rational_point(read_cache = False))
            sage: c.rational_point()
            (-1 : 1 : 0)
        """
        if isinstance(p, (tuple, list)):
            p = self.point(p)
        self._rational_point = p

    def coefficients(self):
        r"""
        Gives a the `6` coefficients of the conic ``self``
        in lexicographic order.

        EXAMPLES::

            sage: Conic(QQ, [1,2,3,4,5,6]).coefficients()
            [1, 2, 3, 4, 5, 6]

            sage: P.<x,y,z> = GF(13)[]
            sage: a = Conic(x^2+5*x*y+y^2+z^2).coefficients(); a
            [1, 5, 0, 1, 0, 1]
            sage: Conic(a)
            Projective Conic Curve over Finite Field of size 13 defined by x^2 + 5*x*y + y^2 + z^2
        """
        return self._coefficients


    def derivative_matrix(self):
        r"""
        Gives the derivative of the defining polynomial of
        the conic ``self``, which is a linear map,
        as a `3 \times 3` matrix.

        EXAMPLES:

        In characteristic different from `2`, the
        derivative matrix is twice the symmetric matrix:

        ::

            sage: c = Conic(QQ, [1,1,1,1,1,0])
            sage: c.symmetric_matrix()
            [  1 1/2 1/2]
            [1/2   1 1/2]
            [1/2 1/2   0]
            sage: c.derivative_matrix()
            [2 1 1]
            [1 2 1]
            [1 1 0]

        An example in characteristic `2`:

        ::

            sage: P.<t> = GF(2)[]
            sage: c = Conic([t, 1, t^2, 1, 1, 0]); c
            Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 2 (using GF2X) defined by t*x^2 + x*y + y^2 + (t^2)*x*z + y*z
            sage: c.is_smooth()
            True
            sage: c.derivative_matrix()
            [  0   1 t^2]
            [  1   0   1]
            [t^2   1   0]
        """
        a, b, c, d, e, f = self.coefficients()
        return Matrix([[ 2*a ,   b ,   c ],
                       [   b , 2*d ,   e ],
                       [   c ,   e , 2*f ]])

    def determinant(self):
        r"""
        Returns the determinant of the symmetric matrix that defines
        the conic ``self``.

        This is defined only if the base field has characteristic
        different from `2`.

        EXAMPLES:

        ::

            sage: C = Conic([1,2,3,4,5,6])
            sage: C.determinant()
            41/4
            sage: C.symmetric_matrix().determinant()
            41/4

        Determinants are only defined in characteristic different from `2`::

            sage: C = Conic(GF(2), [1, 1, 1, 1, 1, 0])
            sage: C.is_smooth()
            True
            sage: C.determinant()
            Traceback (most recent call last):
            ...
            ValueError: The conic self (= Projective Conic Curve over Finite Field of size 2 defined by x^2 + x*y + y^2 + x*z + y*z) has no symmetric matrix because the base field has characteristic 2
        """
        return self.symmetric_matrix().determinant()

    def diagonal_matrix(self):
        r"""
        Returns a diagonal matrix `D` and a matrix `T` such that `T^t A T = D`
        holds, where `(x, y, z) A (x, y, z)^t` is the defining polynomial
        of the conic ``self``.

        EXAMPLES:

        ::

            sage: c = Conic(QQ, [1,2,3,4,5,6])
            sage: d, t = c.diagonal_matrix(); d, t
            (
            [    1     0     0]  [   1   -1 -7/6]
            [    0     3     0]  [   0    1 -1/3]
            [    0     0 41/12], [   0    0    1]
            )
            sage: t.transpose()*c.symmetric_matrix()*t
            [    1     0     0]
            [    0     3     0]
            [    0     0 41/12]

        Diagonal matrices are only defined in characteristic different
        from `2`:

        ::

            sage: c = Conic(GF(4, 'a'), [0, 1, 1, 1, 1, 1])
            sage: c.is_smooth()
            True
            sage: c.diagonal_matrix()
            Traceback (most recent call last):
            ...
            ValueError: The conic self (= Projective Conic Curve over Finite Field in a of size 2^2 defined by x*y + y^2 + x*z + y*z + z^2) has no symmetric matrix because the base field has characteristic 2
        """
        A = self.symmetric_matrix()
        B = self.base_ring()
        basis = [vector(B,{2:0,i:1}) for i in range(3)]
        for i in range(3):
            zerovalue = (basis[i]*A*basis[i].column()== 0)
            if zerovalue:
                for j in range(i+1,3):
                    if basis[j]*A*basis[j].column() != 0:
                        b = basis[i]
                        basis[i] = basis[j]
                        basis[j] = b
                        zerovalue = False
            if zerovalue:
                for j in range(i+1,3):
                    if basis[i]*A*basis[j].column() != 0:
                        basis[i] = basis[i]+basis[j]
                        zerovalue = False
            if not zerovalue:
                l = (basis[i]*A*basis[i].column())
                for j in range(i+1,3):
                    basis[j] = basis[j] - \
                               (basis[i]*A*basis[j].column())/l * basis[i]
        T = Matrix(basis).transpose()
        return T.transpose()*A*T, T

    def diagonalization(self, names=None):
        r"""
        Returns a diagonal conic `C`, an isomorphism of schemes `M: C` -> ``self``
        and the inverse `N` of `M`.

        EXAMPLES::

            sage: Conic(GF(5), [1,0,1,1,0,1]).diagonalization()
            (Projective Conic Curve over Finite Field of size 5 defined by x^2 + y^2 + 2*z^2,
             Scheme morphism:
              From: Projective Conic Curve over Finite Field of size 5 defined by x^2 + y^2 + 2*z^2
              To:   Projective Conic Curve over Finite Field of size 5 defined by x^2 + y^2 + x*z + z^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x + 2*z : y : z),
             Scheme morphism:
              From: Projective Conic Curve over Finite Field of size 5 defined by x^2 + y^2 + x*z + z^2
              To:   Projective Conic Curve over Finite Field of size 5 defined by x^2 + y^2 + 2*z^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x - 2*z : y : z))

        The diagonalization is only defined in characteristic different
        from 2:

        ::

            sage: Conic(GF(2), [1,1,1,1,1,0]).diagonalization()
            Traceback (most recent call last):
            ...
            ValueError: The conic self (= Projective Conic Curve over Finite Field of size 2 defined by x^2 + x*y + y^2 + x*z + y*z) has no symmetric matrix because the base field has characteristic 2

        An example over a global function field:

        ::

            sage: K = FractionField(PolynomialRing(GF(7), 't'))
            sage: (t,) = K.gens()
            sage: C = Conic(K, [t/2,0, 1, 2, 0, 3])
            sage: C.diagonalization()
            (Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7 defined by (-3*t)*x^2 + 2*y^2 + (3*t + 3)/t*z^2,
             Scheme morphism:
               From: Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7 defined by (-3*t)*x^2 + 2*y^2 + (3*t + 3)/t*z^2
               To:   Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7 defined by (-3*t)*x^2 + 2*y^2 + x*z + 3*z^2
               Defn: Defined on coordinates by sending (x : y : z) to
                     (x - 1/t*z : y : z),
             Scheme morphism:
               From: Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7 defined by (-3*t)*x^2 + 2*y^2 + x*z + 3*z^2
               To:   Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7 defined by (-3*t)*x^2 + 2*y^2 + (3*t + 3)/t*z^2
               Defn: Defined on coordinates by sending (x : y : z) to
                     (x + 1/t*z : y : z))


        """
        if names is None:
            names = self.defining_polynomial().parent().variable_names()
        from .constructor import Conic
        D, T = self.diagonal_matrix()
        con = Conic(D, names = names)
        return con, con.hom(T, self), self.hom(T.inverse(), con)

    def gens(self):
        r"""
        Returns the generators of the coordinate ring of ``self``.

        EXAMPLES:

        ::

            sage: P.<x,y,z> = QQ[]
            sage: c = Conic(x^2+y^2+z^2)
            sage: c.gens()
            (xbar, ybar, zbar)
            sage: c.defining_polynomial()(c.gens())
            0

        The function ``gens()`` is required for the following construction:

        ::

            sage: C.<a,b,c> = Conic(GF(3), [1, 1, 1])
            sage: C
            Projective Conic Curve over Finite Field of size 3 defined by a^2 + b^2 + c^2

        """
        return self.coordinate_ring().gens()

    def has_rational_point(self, point = False,
                           algorithm = 'default', read_cache = True):
        r"""
        Returns True if and only if the conic ``self``
        has a point over its base field `B`.

        If ``point`` is True, then returns a second output, which is
        a rational point if one exists.

        Points are cached whenever they are found. Cached information
        is used if and only if ``read_cache`` is True.

        ALGORITHM:

        The parameter ``algorithm`` specifies the algorithm
        to be used:

         - ``'default'`` -- If the base field is real or complex,
           use an elementary native Sage implementation.

         - ``'magma'`` (requires Magma to be installed) --
           delegates the task to the Magma computer algebra
           system.

        EXAMPLES::

            sage: Conic(RR, [1, 1, 1]).has_rational_point()
            False
            sage: Conic(CC, [1, 1, 1]).has_rational_point()
            True

            sage: Conic(RR, [1, 2, -3]).has_rational_point(point = True)
            (True, (1.73205080756888 : 0.000000000000000 : 1.00000000000000))

        Conics over polynomial rings can be solved internally::

            sage: R.<t> = QQ[]
            sage: C = Conic([-2,t^2+1,t^2-1])
            sage: C.has_rational_point()
            True

        And they can also be solved with Magma::

            sage: C.has_rational_point(algorithm='magma') # optional - magma
            True
            sage: C.has_rational_point(algorithm='magma', point=True) # optional - magma
            (True, (-t : 1 : 1))

            sage: D = Conic([t,1,t^2])
            sage: D.has_rational_point(algorithm='magma') # optional - magma
            False

        TESTS:

        One of the following fields comes with an embedding into the complex
        numbers, one does not. Check that they are both handled correctly by
        the Magma interface. ::

            sage: K.<i> = QuadraticField(-1)
            sage: K.coerce_embedding()
            Generic morphism:
              From: Number Field in i with defining polynomial x^2 + 1 with i = 1*I
              To:   Complex Lazy Field
              Defn: i -> 1*I
            sage: Conic(K, [1,1,1]).rational_point(algorithm='magma') # optional - magma
            (-i : 1 : 0)

            sage: x = QQ['x'].gen()
            sage: L.<i> = NumberField(x^2+1, embedding=None)
            sage: Conic(L, [1,1,1]).rational_point(algorithm='magma') # optional - magma
            (-i : 1 : 0)
            sage: L == K
            False
        """
        if read_cache:
            if self._rational_point is not None:
                if point:
                    return True, self._rational_point
                else:
                    return True

        B = self.base_ring()

        if algorithm == 'magma':
            from sage.interfaces.magma import magma
            M = magma(self)
            b = M.HasRationalPoint().sage()
            if not point:
                return b
            if not b:
                return False, None
            M_pt = M.HasRationalPoint(nvals=2)[1]

            # Various attempts will be made to convert `pt` to
            # a Sage object. The end result will always be checked
            # by self.point().

            pt = [M_pt[1], M_pt[2], M_pt[3]]

            # The first attempt is to use sequences. This is efficient and
            # succeeds in cases where the Magma interface fails to convert
            # number field elements, because embeddings between number fields
            # may be lost on conversion to and from Magma.
            # This should deal with all absolute number fields.
            try:
                return True, self.point([B(c.Eltseq().sage()) for c in pt])
            except TypeError:
                pass

            # The second attempt tries to split Magma elements into
            # numerators and denominators first. This is necessary
            # for the field of rational functions, because (at the moment of
            # writing) fraction field elements are not converted automatically
            # from Magma to Sage.
            try:
                return True, self.point( \
                  [B(c.Numerator().sage()/c.Denominator().sage()) for c in pt])
            except (TypeError, NameError):
                pass

            # Finally, let the Magma interface handle conversion.
            try:
                return True, self.point([B(c.sage()) for c in pt])
            except (TypeError, NameError):
                pass

            raise NotImplementedError("No correct conversion implemented for converting the Magma point %s on %s to a correct Sage point on self (=%s)" % (M_pt, M, self))

        if algorithm != 'default':
            raise ValueError("Unknown algorithm: %s" % algorithm)

        if isinstance(B, sage.rings.abc.ComplexField):
            if point:
                [_,_,_,d,e,f] = self._coefficients
                if d == 0:
                    return True, self.point([0,1,0])
                return True, self.point([0, ((e**2-4*d*f).sqrt()-e)/(2*d), 1],
                                        check = False)
            return True
        if isinstance(B, sage.rings.abc.RealField):
            D, T = self.diagonal_matrix()
            [a, b, c] = [D[0,0], D[1,1], D[2,2]]
            if a == 0:
                ret = True, self.point(T*vector([1,0,0]), check = False)
            elif a*c <= 0:
                ret = True, self.point(T*vector([(-c/a).sqrt(),0,1]),
                                       check = False)
            elif b == 0:
                ret = True, self.point(T*vector([0,1,0]), check = False)
            elif b*c <= 0:
                ret = True, self.point(T*vector([0,(-c/b).sqrt(),0,1]),
                                       check = False)
            else:
                ret = False, None
            if point:
                return ret
            return ret[0]
        raise NotImplementedError("has_rational_point not implemented for " \
                                   "conics over base field %s" % B)

    def has_singular_point(self, point = False):
        r"""
        Return True if and only if the conic ``self`` has a rational
        singular point.

        If ``point`` is True, then also return a rational singular
        point (or ``None`` if no such point exists).

        EXAMPLES:

        ::

            sage: c = Conic(QQ, [1,0,1]); c
            Projective Conic Curve over Rational Field defined by x^2 + z^2
            sage: c.has_singular_point(point = True)
            (True, (0 : 1 : 0))

            sage: P.<x,y,z> = GF(7)[]
            sage: e = Conic((x+y+z)*(x-y+2*z)); e
            Projective Conic Curve over Finite Field of size 7 defined by x^2 - y^2 + 3*x*z + y*z + 2*z^2
            sage: e.has_singular_point(point = True)
            (True, (2 : 4 : 1))

            sage: Conic([1, 1, -1]).has_singular_point()
            False
            sage: Conic([1, 1, -1]).has_singular_point(point = True)
            (False, None)

        ``has_singular_point`` is not implemented over all fields
        of characteristic `2`. It is implemented over finite fields.

        ::

            sage: F.<a> = FiniteField(8)
            sage: Conic([a, a+1, 1]).has_singular_point(point = True)
            (True, (a + 1 : 0 : 1))

            sage: P.<t> = GF(2)[]
            sage: C = Conic(P, [t,t,1]); C
            Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 2 (using GF2X) defined by t*x^2 + t*y^2 + z^2
            sage: C.has_singular_point(point = False)
            Traceback (most recent call last):
            ...
            NotImplementedError: Sorry, find singular point on conics not implemented over all fields of characteristic 2.
        """
        if not point:
           ret = self.has_singular_point(point = True)
           return ret[0]
        B = self.base_ring()
        if B.characteristic() == 2:
            [a,b,c,d,e,f] = self.coefficients()
            if b == 0 and c == 0 and e == 0:
                for i in range(3):
                    if [a, d, f][i] == 0:
                        return True, self.point(vector(B, {2:0, i:1}))
                if hasattr(a/f, 'is_square') and hasattr(a/f, 'sqrt'):
                    if (a/f).is_square():
                        return True, self.point([1,0,(a/f).sqrt()])
                    if (d/f).is_square():
                        return True, self.point([0,1,(d/f).sqrt()])
                raise NotImplementedError("Sorry, find singular point on conics not implemented over all fields of characteristic 2.")
            pt = [e, c, b]
            if self.defining_polynomial()(pt) == 0:
                return True, self.point(pt)
            return False, None
        D = self.symmetric_matrix()
        if D.determinant() == 0:
            return True, self.point(Sequence(D.right_kernel().gen()))
        return False, None

    def hom(self, x, Y=None):
        r"""
        Return the scheme morphism from ``self`` to ``Y`` defined by ``x``.
        Here ``x`` can be a matrix or a sequence of polynomials.
        If ``Y`` is omitted, then a natural image is found if possible.

        EXAMPLES:

        Here are a few Morphisms given by matrices. In the first
        example, ``Y`` is omitted, in the second example, ``Y`` is specified.

        ::

            sage: c = Conic([-1, 1, 1])
            sage: h = c.hom(Matrix([[1,1,0],[0,1,0],[0,0,1]])); h
            Scheme morphism:
              From: Projective Conic Curve over Rational Field defined by -x^2 + y^2 + z^2
              To:   Projective Conic Curve over Rational Field defined by -x^2 + 2*x*y + z^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x + y : y : z)
            sage: h([-1, 1, 0])
            (0 : 1 : 0)

            sage: c = Conic([-1, 1, 1])
            sage: d = Conic([4, 1, -1])
            sage: c.hom(Matrix([[0, 0, 1/2], [0, 1, 0], [1, 0, 0]]), d)
            Scheme morphism:
              From: Projective Conic Curve over Rational Field defined by -x^2 + y^2 + z^2
              To:   Projective Conic Curve over Rational Field defined by 4*x^2 + y^2 - z^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (1/2*z : y : x)

        ``ValueError`` is raised if the wrong codomain ``Y`` is specified:

        ::

            sage: c = Conic([-1, 1, 1])
            sage: c.hom(Matrix([[0, 0, 1/2], [0, 1, 0], [1, 0, 0]]), c)
            Traceback (most recent call last):
            ...
            ValueError: The matrix x (= [  0   0 1/2]
            [  0   1   0]
            [  1   0   0]) does not define a map from self (= Projective Conic Curve over Rational Field defined by -x^2 + y^2 + z^2) to Y (= Projective Conic Curve over Rational Field defined by -x^2 + y^2 + z^2)

        The identity map between two representations of the same conic:

        ::

            sage: C = Conic([1,2,3,4,5,6])
            sage: D = Conic([2,4,6,8,10,12])
            sage: C.hom(identity_matrix(3), D)
            Scheme morphism:
              From: Projective Conic Curve over Rational Field defined by x^2 + 2*x*y + 4*y^2 + 3*x*z + 5*y*z + 6*z^2
              To:   Projective Conic Curve over Rational Field defined by 2*x^2 + 4*x*y + 8*y^2 + 6*x*z + 10*y*z + 12*z^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x : y : z)

        An example not over the rational numbers:

        ::

            sage: P.<t> = QQ[]
            sage: C = Conic([1,0,0,t,0,1/t])
            sage: D = Conic([1/t^2, 0, -2/t^2, t, 0, (t + 1)/t^2])
            sage: T = Matrix([[t,0,1],[0,1,0],[0,0,1]])
            sage: C.hom(T, D)
            Scheme morphism:
              From: Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Rational Field defined by x^2 + t*y^2 + 1/t*z^2
              To:   Projective Conic Curve over Fraction Field of Univariate Polynomial Ring in t over Rational Field defined by 1/(t^2)*x^2 + t*y^2 - 2/(t^2)*x*z + (t + 1)/(t^2)*z^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (t*x + z : y : z)

        """
        if is_Matrix(x):
            from .constructor import Conic
            y = x.inverse()
            A = y.transpose()*self.matrix()*y
            im = Conic(A)
            if Y is None:
                Y = im
            elif not Y == im:
                raise ValueError("The matrix x (= %s) does not define a " \
                                 "map from self (= %s) to Y (= %s)" % \
                                 (x, self, Y))
            x = Sequence(x*vector(self.ambient_space().gens()))
            return self.Hom(Y)(x, check = False)
        return ProjectivePlaneCurve.hom(self, x, Y)

    def is_diagonal(self):
        r"""
        Return True if and only if the conic has the form
        `a*x^2 + b*y^2 + c*z^2`.

        EXAMPLES:

        ::

            sage: c=Conic([1,1,0,1,0,1]); c
            Projective Conic Curve over Rational Field defined by x^2 + x*y + y^2 + z^2
            sage: d,t = c.diagonal_matrix()
            sage: c.is_diagonal()
            False
            sage: c.diagonalization()[0].is_diagonal()
            True
        """
        return all(self.coefficients()[i] == 0 for i in [1, 2, 4])

    def is_smooth(self):
        r"""
        Returns True if and only if ``self`` is smooth.

        EXAMPLES:

        ::

            sage: Conic([1,-1,0]).is_smooth()
            False
            sage: Conic(GF(2),[1,1,1,1,1,0]).is_smooth()
            True
        """
        if self.base_ring().characteristic() == 2:
            [a,b,c,d,e,f] = self.coefficients()
            if b == 0 and c == 0 and e == 0:
                return False
            return self.defining_polynomial()([e, c, b]) != 0
        return self.determinant() != 0


    def _magma_init_(self, magma):
        """
        Internal function. Returns a string to initialize this
        conic in the Magma subsystem.

        EXAMPLES::

            sage: C = Conic(QQ, [1,2,3])
            sage: C._magma_init_(magma)                     # optional - magma
            'Conic([_sage_ref...|1/1,2/1,3/1,0/1,0/1,0/1])'
            sage: C = Conic(GF(41), [-1,2,5])               # optional - magma
            sage: C._magma_init_(magma)                     # optional - magma
            'Conic([_sage_ref...|GF(41)!40,GF(41)!2,GF(41)!5,GF(41)!0,GF(41)!0,GF(41)!0])'
            sage: F.<a> = GF(25)
            sage: C = Conic([3,0,1,4,a,2])
            sage: C
            Projective Conic Curve over Finite Field in a of size 5^2 defined by -2*x^2 - y^2 + x*z + a*y*z + 2*z^2
            sage: magma(C)                                  # optional - magma
            Conic over GF(5^2) defined by
            3*X^2 + 4*Y^2 + X*Z + a*Y*Z + 2*Z^2
            sage: magma(Conic([1/2,2/3,-4/5,6/7,8/9,-10/11])) # optional - magma
            Conic over Rational Field defined by
            1/2*X^2 + 2/3*X*Y + 6/7*Y^2 - 4/5*X*Z + 8/9*Y*Z - 10/11*Z^2
            sage: R.<x> = Frac(QQ['x'])
            sage: magma(Conic([x,1+x,1-x]))                 # optional - magma
            Conic over Univariate rational function field over Rational Field defined by
            x*X^2 + (x + 1)*Y^2 + (-x + 1)*Z^2
            sage: P.<x> = QQ[]
            sage: K.<b> = NumberField(x^3+x+1)
            sage: magma(Conic([b,1,2]))                     # optional - magma
            Conic over Number Field with defining polynomial x^3 + x + 1 over the Rational Field defined by
            b*X^2 + Y^2 + 2*Z^2
        """
        kmn = magma(self.base_ring())._ref()
        coeffs = self.coefficients()
        magma_coeffs = [coeffs[i]._magma_init_(magma) for i in [0, 3, 5, 1, 4, 2]]
        return 'Conic([%s|%s])' % (kmn,','.join(magma_coeffs))


    def matrix(self):
        r"""
        Returns a matrix `M` such that `(x, y, z) M (x, y, z)^t`
        is the defining equation of ``self``.

        The matrix `M` is upper triangular if the base field has
        characteristic `2` and symmetric otherwise.

        EXAMPLES::

            sage: R.<x, y, z> = QQ[]
            sage: C = Conic(x^2 + x*y + y^2 + z^2)
            sage: C.matrix()
            [  1 1/2   0]
            [1/2   1   0]
            [  0   0   1]

            sage: R.<x, y, z> = GF(2)[]
            sage: C = Conic(x^2 + x*y + y^2 + x*z + z^2)
            sage: C.matrix()
            [1 1 1]
            [0 1 0]
            [0 0 1]
        """
        if self.base_ring().characteristic() == 2:
            return self.upper_triangular_matrix()
        return self.symmetric_matrix()

    _matrix_ = matrix

    def parametrization(self, point=None, morphism=True):
        r"""
        Return a parametrization `f` of ``self`` together with the
        inverse of `f`.

        .. warning::

           The second map is currently broken and neither the inverse nor
           well-defined.

        If ``point`` is specified, then that point is used
        for the parametrization. Otherwise, use ``self.rational_point()``
        to find a point.

        If ``morphism`` is True, then `f` is returned in the form
        of a Scheme morphism. Otherwise, it is a tuple of polynomials
        that gives the parametrization.

        EXAMPLES:

        An example over a finite field ::

            sage: c = Conic(GF(2), [1,1,1,1,1,0])
            sage: f, g = c.parametrization(); f, g
            (Scheme morphism:
              From: Projective Space of dimension 1 over Finite Field of size 2
              To:   Projective Conic Curve over Finite Field of size 2 defined by x^2 + x*y
            + y^2 + x*z + y*z
              Defn: Defined on coordinates by sending (x : y) to ...,
             Scheme morphism:
              From: Projective Conic Curve over Finite Field of size 2 defined by x^2 + x*y
            + y^2 + x*z + y*z
              To:   Projective Space of dimension 1 over Finite Field of size 2
              Defn: Defined on coordinates by sending (x : y : z) to ...)
            sage: set(f(p) for p in f.domain())
            {(0 : 0 : 1), (0 : 1 : 1), (1 : 0 : 1)}
            sage: (g*f).is_one()  # known bug  (see :trac:`31892`)
            True

        An example with ``morphism = False`` ::

            sage: R.<x,y,z> = QQ[]
            sage: C = Curve(7*x^2 + 2*y*z + z^2)
            sage: (p, i) = C.parametrization(morphism = False); (p, i)
            ([-2*x*y, x^2 + 7*y^2, -2*x^2], [-1/2*x, 1/7*y + 1/14*z])
            sage: C.defining_polynomial()(p)
            0
            sage: i[0](p) / i[1](p)
            x/y

        A ``ValueError`` is raised if ``self`` has no rational point ::

            sage: C = Conic(x^2 + y^2 + 7*z^2)
            sage: C.parametrization()
            Traceback (most recent call last):
            ...
            ValueError: Conic Projective Conic Curve over Rational Field defined by x^2 + y^2 + 7*z^2 has no rational points over Rational Field!

        A ``ValueError`` is raised if ``self`` is not smooth ::

            sage: C = Conic(x^2 + y^2)
            sage: C.parametrization()
            Traceback (most recent call last):
            ...
            ValueError: The conic self (=Projective Conic Curve over Rational Field defined by x^2 + y^2) is not smooth, hence does not have a parametrization.
        """
        if (self._parametrization is not None) and not point:
            par = self._parametrization
        else:
            if not self.is_smooth():
                raise ValueError("The conic self (=%s) is not smooth, hence does not have a parametrization." % self)
            if point is None:
                point = self.rational_point()
            point = Sequence(point)
            B = self.base_ring()
            Q = PolynomialRing(B, 'x,y')
            [x, y] = Q.gens()
            gens = self.ambient_space().gens()
            P = PolynomialRing(B, 4, ['X', 'Y', 'T0', 'T1'])
            [X, Y, T0, T1] = P.gens()
            c3 = [j for j in range(2,-1,-1) if point[j] != 0][0]
            c1 = [j for j in range(3) if j != c3][0]
            c2 = [j for j in range(3) if j != c3 and j != c1][0]
            L = [0,0,0]
            L[c1] = Y*T1*point[c1] + Y*T0
            L[c2] = Y*T1*point[c2] + X*T0
            L[c3] = Y*T1*point[c3]
            bezout = P(self.defining_polynomial()(L) / T0)
            t = [bezout([x,y,0,-1]),bezout([x,y,1,0])]
            par = (tuple([Q(p([x,y,t[0],t[1]])/y) for  p in L]),
                   tuple([gens[m]*point[c3]-gens[c3]*point[m]
                       for m in [c2,c1]]))
            if self._parametrization is None:
                self._parametrization = par
        if not morphism:
            return par
        P1 = ProjectiveSpace(self.base_ring(), 1, 'x,y')
        return P1.hom(par[0],self), self.Hom(P1)(par[1], check = False)

    def point(self, v, check=True):
        r"""
        Constructs a point on ``self`` corresponding to the input ``v``.

        If ``check`` is True, then checks if ``v`` defines a valid
        point on ``self``.

        If no rational point on ``self`` is known yet, then also caches the point
        for use by ``self.rational_point()`` and ``self.parametrization()``.

        EXAMPLES::

            sage: c = Conic([1, -1, 1])
            sage: c.point([15, 17, 8])
            (15/8 : 17/8 : 1)
            sage: c.rational_point()
            (15/8 : 17/8 : 1)
            sage: d = Conic([1, -1, 1])
            sage: d.rational_point()
            (-1 : 1 : 0)
        """
        if is_Vector(v):
            v = Sequence(v)
        p = ProjectivePlaneCurve.point(self, v, check=check)
        if self._rational_point is None:
            self._rational_point = p
        return p


    def random_rational_point(self, *args1, **args2):
        r"""
        Return a random rational point of the conic ``self``.

        ALGORITHM:

            1. Compute a parametrization `f` of ``self`` using
               ``self.parametrization()``.
            2. Computes a random point `(x:y)` on the projective
               line.
            3. Output `f(x:y)`.

        The coordinates x and y are computed using
        ``B.random_element``, where ``B`` is the base field of
        ``self`` and additional arguments to ``random_rational_point``
        are passed to ``random_element``.

        If the base field is a finite field, then the
        output is uniformly distributed over the points of self.

        EXAMPLES::

            sage: c = Conic(GF(2), [1,1,1,1,1,0])
            sage: [c.random_rational_point() for i in range(10)] # output is random
            [(1 : 0 : 1), (1 : 0 : 1), (1 : 0 : 1), (0 : 1 : 1), (1 : 0 : 1), (0 : 0 : 1), (1 : 0 : 1), (1 : 0 : 1), (0 : 0 : 1), (1 : 0 : 1)]

            sage: d = Conic(QQ, [1, 1, -1])
            sage: d.random_rational_point(den_bound = 1, num_bound = 5) # output is random
            (-24/25 : 7/25 : 1)

            sage: Conic(QQ, [1, 1, 1]).random_rational_point()
            Traceback (most recent call last):
            ...
            ValueError: Conic Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2 has no rational points over Rational Field!

        """
        if not self.is_smooth():
            raise NotImplementedError("Sorry, random points not implemented " \
                                       "for non-smooth conics")
        par = self.parametrization()
        x = 0
        y = 0
        B = self.base_ring()
        while x == 0 and y == 0:
            x = B.random_element(*args1, **args2)
            y = B.random_element(*args1, **args2)
        return par[0]([x,y])


    def rational_point(self, algorithm = 'default', read_cache = True):
        r"""
        Return a point on ``self`` defined over the base field.

        Raises ``ValueError`` if no rational point exists.

        See ``self.has_rational_point`` for the algorithm used
        and for the use of the parameters ``algorithm`` and ``read_cache``.

        EXAMPLES:

        Examples over `\QQ` ::

            sage: R.<x,y,z> = QQ[]
            sage: C = Conic(7*x^2 + 2*y*z + z^2)
            sage: C.rational_point()
            (0 : 1 : 0)

            sage: C = Conic(x^2 + 2*y^2 + z^2)
            sage: C.rational_point()
            Traceback (most recent call last):
            ...
            ValueError: Conic Projective Conic Curve over Rational Field defined by x^2 + 2*y^2 + z^2 has no rational points over Rational Field!

            sage: C = Conic(x^2 + y^2 + 7*z^2)
            sage: C.rational_point(algorithm = 'rnfisnorm')
            Traceback (most recent call last):
            ...
            ValueError: Conic Projective Conic Curve over Rational Field defined by x^2 + y^2 + 7*z^2 has no rational points over Rational Field!

        Examples over number fields ::

            sage: P.<x> = QQ[]
            sage: L.<b> = NumberField(x^3-5)
            sage: C = Conic(L, [3, 2, -b])
            sage: p = C.rational_point(algorithm = 'rnfisnorm')
            sage: p                                         # output is random
            (1/3*b^2 - 4/3*b + 4/3 : b^2 - 2 : 1)
            sage: C.defining_polynomial()(list(p))
            0

            sage: K.<i> = QuadraticField(-1)
            sage: D = Conic(K, [3, 2, 5])
            sage: D.rational_point(algorithm = 'rnfisnorm') # output is random
            (-3 : 4*i : 1)

            sage: L.<s> = QuadraticField(2)
            sage: Conic(QQ, [1, 1, -3]).has_rational_point()
            False
            sage: E = Conic(L, [1, 1, -3])
            sage: E.rational_point()                        # output is random
            (-1 : -s : 1)

        Currently Magma is better at solving conics over number fields than
        Sage, so it helps to use the algorithm 'magma' if Magma is installed::

            sage: q = C.rational_point(algorithm = 'magma', read_cache=False) # optional - magma
            sage: q                       # output is random, optional - magma
            (1/5*b^2 : 1/5*b^2 : 1)
            sage: C.defining_polynomial()(list(p))          # optional - magma
            0
            sage: len(str(p)) > 1.5*len(str(q))             # optional - magma
            True

            sage: D.rational_point(algorithm = 'magma', read_cache=False) # random, optional - magma
            (1 : 2*i : 1)

            sage: E.rational_point(algorithm='magma', read_cache=False) # random, optional - magma
            (-s : 1 : 1)

            sage: F = Conic([L.gen(), 30, -20])
            sage: q = F.rational_point(algorithm='magma')   # optional - magma
            sage: q                       # output is random, optional - magma
            (-10/7*s + 40/7 : 5/7*s - 6/7 : 1)
            sage: p = F.rational_point(read_cache=False)
            sage: p                       # output is random
            (788210*s - 1114700 : -171135*s + 242022 : 1)
            sage: len(str(p)) > len(str(q))                 # optional - magma
            True

            sage: Conic([L.gen(), 30, -21]).has_rational_point(algorithm='magma') # optional - magma
            False

        Examples over finite fields ::

            sage: F.<a> = FiniteField(7^20)
            sage: C = Conic([1, a, -5]); C
            Projective Conic Curve over Finite Field in a of size 7^20 defined by x^2 + a*y^2 + 2*z^2
            sage: C.rational_point()  # output is random
            (4*a^19 + 5*a^18 + 4*a^17 + a^16 + 6*a^15 + 3*a^13 + 6*a^11 + a^9 + 3*a^8 + 2*a^7 + 4*a^6 + 3*a^5 + 3*a^4 + a^3 + a + 6 : 5*a^18 + a^17 + a^16 + 6*a^15 + 4*a^14 + a^13 + 5*a^12 + 5*a^10 + 2*a^9 + 6*a^8 + 6*a^7 + 6*a^6 + 2*a^4 + 3 : 1)

        Examples over `\RR` and `\CC` ::

            sage: Conic(CC, [1, 2, 3]).rational_point()
            (0 : 1.22474487139159*I : 1)

            sage: Conic(RR, [1, 1, 1]).rational_point()
            Traceback (most recent call last):
            ...
            ValueError: Conic Projective Conic Curve over Real Field with 53 bits of precision defined by x^2 + y^2 + z^2 has no rational points over Real Field with 53 bits of precision!
        """
        bl,pt = self.has_rational_point(point = True, algorithm = algorithm,
                                        read_cache = read_cache)
        if bl:
            return pt
        raise ValueError("Conic %s has no rational points over %s!" % \
                          (self, self.ambient_space().base_ring()))


    def singular_point(self):
        r"""
        Returns a singular rational point of ``self``

        EXAMPLES:

        ::

            sage: Conic(GF(2), [1,1,1,1,1,1]).singular_point()
            (1 : 1 : 1)

        ``ValueError`` is raised if the conic has no rational singular point

        ::

            sage: Conic(QQ, [1,1,1,1,1,1]).singular_point()
            Traceback (most recent call last):
            ...
            ValueError: The conic self (= Projective Conic Curve over Rational Field defined by x^2 + x*y + y^2 + x*z + y*z + z^2) has no rational singular point
        """
        b = self.has_singular_point(point = True)
        if not b[0]:
            raise ValueError("The conic self (= %s) has no rational " \
                              "singular point" % self)
        return b[1]

    def symmetric_matrix(self):
        r"""
        The symmetric matrix `M` such that `(x y z) M (x y z)^t`
        is the defining equation of ``self``.

        EXAMPLES::

            sage: R.<x, y, z> = QQ[]
            sage: C = Conic(x^2 + x*y/2 + y^2 + z^2)
            sage: C.symmetric_matrix()
            [  1 1/4   0]
            [1/4   1   0]
            [  0   0   1]

            sage: C = Conic(x^2 + 2*x*y + y^2 + 3*x*z + z^2)
            sage: v = vector([x, y, z])
            sage: v * C.symmetric_matrix() * v
            x^2 + 2*x*y + y^2 + 3*x*z + z^2
        """
        a, b, c, d, e, f = self.coefficients()
        if self.base_ring().characteristic() == 2:
            if b == 0 and c == 0 and e == 0:
                return Matrix([[a,0,0],[0,d,0],[0,0,f]])
            raise ValueError("The conic self (= %s) has no symmetric matrix " \
                              "because the base field has characteristic 2" % \
                              self)
        return Matrix([[  a , b/2, c/2 ],
                       [ b/2,  d , e/2 ],
                       [ c/2, e/2,  f  ]])


    def upper_triangular_matrix(self):
        r"""
        The upper-triangular matrix `M` such that `(x y z) M (x y z)^t`
        is the defining equation of ``self``.

        EXAMPLES::

            sage: R.<x, y, z> = QQ[]
            sage: C = Conic(x^2 + x*y + y^2 + z^2)
            sage: C.upper_triangular_matrix()
            [1 1 0]
            [0 1 0]
            [0 0 1]

            sage: C = Conic(x^2 + 2*x*y + y^2 + 3*x*z + z^2)
            sage: v = vector([x, y, z])
            sage: v * C.upper_triangular_matrix() * v
            x^2 + 2*x*y + y^2 + 3*x*z + z^2
        """
        from sage.matrix.constructor import matrix
        [a,b,c,d,e,f] = self.coefficients()
        return matrix([[ a, b, c ],
                       [ 0, d, e ],
                       [ 0, 0, f ]])

    def variable_names(self):
        r"""
        Returns the variable names of the defining polynomial
        of ``self``.

        EXAMPLES:

        ::

            sage: c=Conic([1,1,0,1,0,1], 'x,y,z')
            sage: c.variable_names()
            ('x', 'y', 'z')
            sage: c.variable_name()
            'x'

        The function ``variable_names()`` is required
        for the following construction:

        ::

            sage: C.<p,q,r> = Conic(QQ, [1, 1, 1])
            sage: C
            Projective Conic Curve over Rational Field defined by p^2 + q^2 + r^2

        """
        return self.defining_polynomial().parent().variable_names()

