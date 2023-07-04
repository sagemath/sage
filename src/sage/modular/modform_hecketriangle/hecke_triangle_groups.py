r"""
Hecke triangle groups

AUTHORS:

- Jonas Jermann (2013): initial version

"""

# ****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import divisors
from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_generic
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.unique_representation import UniqueRepresentation

lazy_import("sage.functions.log", "exp")
lazy_import("sage.functions.gamma", "psi1")
lazy_import("sage.functions.trig", "sec")
lazy_import("sage.rings.number_field.number_field", "NumberField")
lazy_import("sage.rings.qqbar", ["AA", "AlgebraicField"])
lazy_import("sage.rings.universal_cyclotomic_field", "E")
lazy_import("sage.symbolic.constants", "pi")

from .hecke_triangle_group_element import HeckeTriangleGroupElement, cyclic_representative, coerce_AA


class HeckeTriangleGroup(FinitelyGeneratedMatrixGroup_generic,
                         UniqueRepresentation):
    r"""
    Hecke triangle group `(2, n, \infty)`.
    """

    Element = HeckeTriangleGroupElement

    @staticmethod
    def __classcall__(cls, n=3):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(QQ(3)) == HeckeTriangleGroup(int(3))
            True
        """
        if n == infinity:
            n = infinity
        else:
            n = ZZ(n)
            if n < 3:
                raise AttributeError("n has to be infinity or an Integer >= 3.")

        return super().__classcall__(cls, n)

    def __init__(self, n):
        r"""
        Hecke triangle group (2, n, infinity).

        Namely the von Dyck group corresponding to the triangle group
        with angles (pi/2, pi/n, 0).

        INPUT:

        - ``n``   - ``infinity`` or an integer greater or equal to ``3``.

        OUTPUT:

        The Hecke triangle group for the given parameter ``n``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(12)
            sage: G
            Hecke triangle group for n = 12
            sage: G.category()
            Category of groups
        """

        self._n = n
        self.element_repr_method("default")

        if n in [3, infinity]:
            self._base_ring = ZZ
            self._lam = ZZ.one() if n == 3 else ZZ(2)
        else:
            lam_symbolic = coerce_AA(E(2 * n) + ~E(2 * n))
            K = NumberField(self.lam_minpoly(), 'lam',
                            embedding=lam_symbolic)
            # self._base_ring = K.order(K.gens())
            self._base_ring = K.maximal_order()
            self._lam = self._base_ring.gen(1)

        T = matrix(self._base_ring, [[1, self._lam], [0, 1]])
        S = matrix(self._base_ring, [[0, -1], [1, 0]])

        FinitelyGeneratedMatrixGroup_generic.__init__(self, ZZ(2), self._base_ring, [S, T])

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10)
            Hecke triangle group for n = 10
        """
        return "Hecke triangle group for n = {}".format(self._n)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: a = HeckeTriangleGroup(5)
            sage: latex(a)
            \Gamma^{(5)}
        """
        return '\\Gamma^{(%s)}' % latex(self._n)

    def element_repr_method(self, method=None):
        r"""
        Either return or set the representation method for elements of ``self``.

        INPUT:

        - ``method``  -- If ``method=None`` (default) the current default representation
                         method is returned. Otherwise the default method is set to ``method``.
                         If ``method`` is not available a ValueError is raised. Possible methods are:

                         ``default``: Use the usual representation method for matrix group elements.

                         ``basic``:   The representation is given as a word in ``S`` and powers of ``T``.

                         ``conj``:    The conjugacy representative of the element is represented
                                      as a word in powers of the basic blocks, together with
                                      an unspecified conjugation matrix.

                         ``block``:   Same as ``conj`` but the conjugation matrix is specified as well.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(5)
            sage: G.element_repr_method()
            'default'
            sage: G.element_repr_method("basic")
            sage: G.element_repr_method()
            'basic'
        """
        if method is None:
            return self._element_repr_method
        elif method in ["default", "basic", "block", "conj"]:
            self._element_repr_method = method
        else:
            raise ValueError(f"the specified method {method} is not supported")

    def one(self):
        r"""
        Return the identity element/matrix for ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(10)
            sage: G(1) == G.one()
            True
            sage: G(1)
            [1 0]
            [0 1]

            sage: G(1).parent()
            Hecke triangle group for n = 10
        """
        return self.I()

    def lam_minpoly(self):
        r"""
        Return the minimal polynomial of the corresponding lambda parameter of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).lam_minpoly()
            x^4 - 5*x^2 + 5
            sage: HeckeTriangleGroup(17).lam_minpoly()
            x^8 - x^7 - 7*x^6 + 6*x^5 + 15*x^4 - 10*x^3 - 10*x^2 + 4*x + 1
            sage: HeckeTriangleGroup(infinity).lam_minpoly()
            x - 2
        """
        # TODO: Write an explicit (faster) implementation
        n = self._n
        if n == infinity:
            lam_symbolic = QQ(2)
        else:
            lam_symbolic = E(2 * n) + ~E(2 * n)
        return lam_symbolic.minpoly()

    def base_ring(self):
        r"""
        Return the base ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(n=infinity).base_ring()
            Integer Ring
            sage: HeckeTriangleGroup(n=7).base_ring()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1 with lam = 1.801937735804839?
        """
        return self._base_ring

    def base_field(self):
        r"""
        Return the base field of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(n=infinity).base_field()
            Rational Field
            sage: HeckeTriangleGroup(n=7).base_field()
            Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1 with lam = 1.801937735804839?
        """
        if self._n in [3, infinity]:
            return QQ
        else:
            return self._base_ring.number_field()

    def n(self):
        r"""
        Return the parameter ``n`` of ``self``, where
        ``pi/n`` is the angle at ``rho`` of the corresponding
        basic hyperbolic triangle with vertices ``i``, ``rho``
        and ``infinity``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).n()
            10
            sage: HeckeTriangleGroup(infinity).n()
            +Infinity
        """
        return self._n

    def lam(self):
        r"""
        Return the parameter ``lambda`` of ``self``,
        where ``lambda`` is twice the real part of ``rho``,
        lying between ``1`` (when ``n=3``) and ``2`` (when ``n=infinity``).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).lam()
            1
            sage: HeckeTriangleGroup(4).lam()
            lam
            sage: HeckeTriangleGroup(4).lam()^2
            2
            sage: HeckeTriangleGroup(6).lam()^2
            3
            sage: AA(HeckeTriangleGroup(10).lam())
            1.9021130325903...?
            sage: HeckeTriangleGroup(infinity).lam()
            2
        """
        return self._lam

    def rho(self):
        r"""
        Return the vertex ``rho`` of the basic hyperbolic
        triangle which describes ``self``. ``rho`` has
        absolute value 1 and angle ``pi/n``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).rho() == QQbar(1/2 + sqrt(3)/2*i)
            True
            sage: HeckeTriangleGroup(4).rho() == QQbar(sqrt(2)/2*(1 + i))
            True
            sage: HeckeTriangleGroup(6).rho() == QQbar(sqrt(3)/2 + 1/2*i)
            True
            sage: HeckeTriangleGroup(10).rho()
            0.95105651629515...? + 0.30901699437494...?*I
            sage: HeckeTriangleGroup(infinity).rho()
            1
        """
        # TODO: maybe rho should be replaced by -rhobar
        # Also we could use NumberFields...
        if self._n == infinity:
            return AA.one()
        else:
            rho = AlgebraicField()(E(2 * self._n))
            rho.simplify()
            return rho

    def alpha(self):
        r"""
        Return the parameter ``alpha`` of ``self``.

        This is the first parameter of the hypergeometric series used
        in the calculation of the Hauptmodul of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).alpha()
            1/12
            sage: HeckeTriangleGroup(4).alpha()
            1/8
            sage: HeckeTriangleGroup(5).alpha()
            3/20
            sage: HeckeTriangleGroup(6).alpha()
            1/6
            sage: HeckeTriangleGroup(10).alpha()
            1/5
            sage: HeckeTriangleGroup(infinity).alpha()
            1/4
        """
        return QQ((1, 2)) * (QQ((1, 2)) - ZZ.one() / self._n)

    def beta(self):
        r"""
        Return the parameter ``beta`` of ``self``.

        This is the second parameter of the hypergeometric series used
        in the calculation of the Hauptmodul of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).beta()
            5/12
            sage: HeckeTriangleGroup(4).beta()
            3/8
            sage: HeckeTriangleGroup(5).beta()
            7/20
            sage: HeckeTriangleGroup(6).beta()
            1/3
            sage: HeckeTriangleGroup(10).beta()
            3/10
            sage: HeckeTriangleGroup(infinity).beta()
            1/4
        """
        return QQ((1, 2)) * (QQ((1, 2)) + ZZ.one() / self._n)

    # We use cached method here to create unique instances of basic matrices
    # (major performance gain)
    @cached_method
    def I(self):
        r"""
        Return the identity element/matrix for ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).I()
            [1 0]
            [0 1]

            sage: HeckeTriangleGroup(10).I().parent()
            Hecke triangle group for n = 10
        """
        return self(matrix(self._base_ring, [[1, 0], [0, 1]]), check=False)

    # We use cached method here to create unique instances of basic matrices
    # (major performance gain)
    @cached_method
    def T(self, m=1):
        r"""
        Return the element in ``self`` corresponding to the translation by ``m*self.lam()``.

        INPUT:

        - ``m`` -- An integer, default: ``1``, namely the second generator of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).T()
            [1 1]
            [0 1]
            sage: HeckeTriangleGroup(10).T(-4)
            [     1 -4*lam]
            [     0      1]
            sage: HeckeTriangleGroup(10).T().parent()
            Hecke triangle group for n = 10
        """
        return self(matrix(self._base_ring, [[1, self._lam * m], [0, 1]]),
                    check=False)

    # We use cached method here to create unique instances of basic matrices
    # (major performance gain)
    @cached_method
    def S(self):
        r"""
        Return the generator of ``self`` corresponding to the
        conformal circle inversion.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).S()
            [ 0 -1]
            [ 1  0]
            sage: HeckeTriangleGroup(10).S()
            [ 0 -1]
            [ 1  0]
            sage: HeckeTriangleGroup(10).S()^2 == -HeckeTriangleGroup(10).I()
            True
            sage: HeckeTriangleGroup(10).S()^4 == HeckeTriangleGroup(10).I()
            True

            sage: HeckeTriangleGroup(10).S().parent()
            Hecke triangle group for n = 10
        """
        return self.gen(0)

    # We use cached method here to create unique instances of basic matrices
    # (major performance gain)
    @cached_method
    def U(self):
        r"""
        Return an alternative generator of ``self`` instead of ``T``.
        ``U`` stabilizes ``rho`` and has order ``2*self.n()``.

        If ``n=infinity`` then ``U`` is parabolic and has infinite order,
        it then fixes the cusp ``[-1]``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).U()
            [ 1 -1]
            [ 1  0]
            sage: HeckeTriangleGroup(3).U()^3 == -HeckeTriangleGroup(3).I()
            True
            sage: HeckeTriangleGroup(3).U()^6 == HeckeTriangleGroup(3).I()
            True
            sage: HeckeTriangleGroup(10).U()
            [lam  -1]
            [  1   0]
            sage: HeckeTriangleGroup(10).U()^10 == -HeckeTriangleGroup(10).I()
            True
            sage: HeckeTriangleGroup(10).U()^20 == HeckeTriangleGroup(10).I()
            True

            sage: HeckeTriangleGroup(10).U().parent()
            Hecke triangle group for n = 10
        """
        return self.T() * self.S()

    def V(self, j):
        r"""
        Return the j'th generator for the usual representatives of
        conjugacy classes of ``self``. It is given by ``V=U^(j-1)*T``.

        INPUT:

        - ``j``  -- Any integer. To get the usual representatives
          ``j`` should range from ``1`` to ``self.n()-1``.

        OUTPUT:

        The corresponding matrix/element.

        The matrix is parabolic if ``j`` is congruent to `\pm 1` modulo ``self.n()``.
        It is elliptic if ``j`` is congruent to 0 modulo ``self.n()``.
        It is hyperbolic otherwise.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(3)
            sage: G.V(0) == -G.S()
            True
            sage: G.V(1) == G.T()
            True
            sage: G.V(2)
            [1 0]
            [1 1]
            sage: G.V(3) == G.S()
            True

            sage: G = HeckeTriangleGroup(5)
            sage: G.element_repr_method("default")
            sage: G.V(1)
            [  1 lam]
            [  0   1]
            sage: G.V(2)
            [lam lam]
            [  1 lam]
            sage: G.V(3)
            [lam   1]
            [lam lam]
            sage: G.V(4)
            [  1   0]
            [lam   1]
            sage: G.V(5) == G.S()
            True
        """
        return self.U()**(j-1) * self.T()

    def dvalue(self):
        r"""
        Return a symbolic expression (or an exact value in case n=3, 4, 6)
        for the transfinite diameter (or capacity) of ``self``.

        This is the first nontrivial Fourier coefficient of the Hauptmodul
        for the Hecke triangle group in case it is normalized to ``J_inv(i)=1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).dvalue()
            1/1728
            sage: HeckeTriangleGroup(4).dvalue()
            1/256
            sage: HeckeTriangleGroup(5).dvalue()
            e^(2*euler_gamma - 4*pi/(sqrt(5) + 1) + psi(17/20) + psi(13/20))
            sage: HeckeTriangleGroup(6).dvalue()
            1/108
            sage: HeckeTriangleGroup(10).dvalue()
            e^(2*euler_gamma - 4*pi/sqrt(2*sqrt(5) + 10) + psi(4/5) + psi(7/10))
            sage: HeckeTriangleGroup(infinity).dvalue()
            1/64
        """
        n = self._n
        if n == 3:
            return 1 / ZZ(2**6*3**3)
        if n == 4:
            return 1 / ZZ(2**8)
        if n == 6:
            return 1 / ZZ(2**2*3**3)
        if n == infinity:
            return 1 / ZZ(2**6)
        # reference for this formula ?
        return exp(-ZZ(2)*psi1(ZZ.one()) + psi1(ZZ.one()-self.alpha())
                   + psi1(ZZ.one()-self.beta()) - pi*sec(pi/self._n))

    def is_arithmetic(self) -> bool:
        r"""
        Return ``True`` if ``self`` is an arithmetic subgroup.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).is_arithmetic()
            True
            sage: HeckeTriangleGroup(4).is_arithmetic()
            True
            sage: HeckeTriangleGroup(5).is_arithmetic()
            False
            sage: HeckeTriangleGroup(6).is_arithmetic()
            True
            sage: HeckeTriangleGroup(10).is_arithmetic()
            False
            sage: HeckeTriangleGroup(infinity).is_arithmetic()
            True
        """
        return self._n in [ZZ(3), ZZ(4), ZZ(6), infinity]

    def get_FD(self, z):
        r"""
        Return a tuple (A,w) which determines how to map ``z``
        to the usual (strict) fundamental domain of ``self``.

        INPUT:

        - ``z`` -- a complex number or an element of AlgebraicField().

        OUTPUT:

        A tuple ``(A, w)``.

        - ``A`` -- a matrix in ``self`` such that ``A.acton(w)==z``
          (if ``z`` is exact at least).

        - ``w`` -- a complex number or an element of AlgebraicField()
          (depending on the type ``z``) which lies inside the (strict)
          fundamental domain of ``self`` (``self.in_FD(w)==True``) and
          which is equivalent to ``z`` (by the above property).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(8)
            sage: z = AlgebraicField()(1+i/2)
            sage: (A, w) = G.get_FD(z)
            sage: A
            [-lam    1]
            [  -1    0]
            sage: A.acton(w) == z
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: z = (134.12 + 0.22*i).n()
            sage: (A, w) = G.get_FD(z)
            sage: A
            [-73*lam^3 + 74*lam       73*lam^2 - 1]
            [        -lam^2 + 1                lam]
            sage: w
            0.769070776942... + 0.779859114103...*I
            sage: z
            134.120000000... + 0.220000000000...*I
            sage: A.acton(w)
            134.1200000... + 0.2200000000...*I
        """
        ID = self.I()
        T = self.T()
        S = self.S()
        TI = self.T(-1)

        A = ID
        w = z
        while abs(w) < ZZ.one() or abs(w.real()) > self.lam() / ZZ(2):
            if abs(w) < ZZ.one():
                w = self.S().acton(w)
                A = S*A
            while w.real() >= self.lam() / ZZ(2):
                w = TI.acton(w)
                A = TI*A
            while w.real() < -self.lam() / ZZ(2):
                w = T.acton(w)
                A = T*A
        if w.real() == self.lam() / ZZ(2):
            w = TI.acton(w)
            A = TI*A
        if abs(w) == ZZ.one() and w.real() > ZZ.zero():
            w = S.acton(w)
            A = S*A

        AI = A.inverse()

        return (AI, A.acton(z))

    def in_FD(self, z) -> bool:
        r"""
        Return ``True`` if ``z`` lies in the (strict) fundamental
        domain of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(5).in_FD(CC(1.5/2 + 0.9*i))
            True
            sage: HeckeTriangleGroup(4).in_FD(CC(1.5/2 + 0.9*i))
            False
        """
        return self.get_FD(z)[0] == self.I()

    def root_extension_field(self, D):
        r"""
        Return the quadratic extension field of the base field by
        the square root of the given discriminant ``D``.

        INPUT:

        - ``D`` -- An element of the base ring of ``self``
          corresponding to a discriminant.

        OUTPUT:

        A relative (at most quadratic) extension to the base field
        of self in the variable ``e`` which corresponds to ``sqrt(D)``.
        If the extension degree is ``1`` then the base field is returned.

        The correct embedding is the positive resp. positive imaginary one.
        Unfortunately no default embedding can be specified for relative
        number fields yet.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)
            sage: G.root_extension_field(32)
            Number Field in e with defining polynomial x^2 - 32
            sage: G.root_extension_field(-4)
            Number Field in e with defining polynomial x^2 + 4
            sage: G.root_extension_field(4) == G.base_field()
            True

            sage: G = HeckeTriangleGroup(n=7)
            sage: lam = G.lam()
            sage: D = 4*lam^2 + 4*lam - 4
            sage: G.root_extension_field(D)
            Number Field in e with defining polynomial x^2 - 4*lam^2 - 4*lam + 4 over its base field
            sage: G.root_extension_field(4) == G.base_field()
            True
            sage: D = lam^2 - 4
            sage: G.root_extension_field(D)
            Number Field in e with defining polynomial x^2 - lam^2 + 4 over its base field
        """
        K = self.base_field()
        x = PolynomialRing(K, 'x').gen()
        D = self.base_ring()(D)

        if D.is_square():
            return K
        else:
            # unfortunately we can't set embeddings for relative extensions :-(
            # return K.extension(x**2 - D, 'e', embedding=AA(D).sqrt())

            L = K.extension(x**2 - D, 'e')

            # e = AA(D).sqrt()
            # emb = L.hom([e])
            # L._unset_embedding()
            # L.register_embedding(emb)

            # return NumberField(L.absolute_polynomial(), 'e', structure=AbsoluteFromRelative(L), embedding=(???))
            return L

    # We cache this method for performance reasons (it is repeatedly reused)
    @cached_method
    def root_extension_embedding(self, D, K=None):
        r"""
        Return the correct embedding from the root extension field
        of the given discriminant ``D``  to the field ``K``.

        Also see the method ``root_extension_embedding(K)`` of
        ``HeckeTriangleGroupElement`` for more examples.

        INPUT:

        - ``D`` -- An element of the base ring of ``self``
          corresponding to a discriminant.

        - ``K`` -- A field to which we want the (correct) embedding.
          If ``K=None`` (default) then ``AlgebraicField()`` is
          used for positive ``D`` and ``AlgebraicRealField()``
          otherwise.

        OUTPUT:

        The corresponding embedding if it was found.
        Otherwise a ``ValueError`` is raised.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)
            sage: G.root_extension_embedding(32)
            Ring morphism:
              From: Number Field in e with defining polynomial x^2 - 32
              To:   Algebraic Real Field
              Defn: e |--> 5.656854249492...?
            sage: G.root_extension_embedding(-4)
            Ring morphism:
              From: Number Field in e with defining polynomial x^2 + 4
              To:   Algebraic Field
              Defn: e |--> 2*I
            sage: G.root_extension_embedding(4)
            Coercion map:
              From: Rational Field
              To:   Algebraic Real Field

            sage: G = HeckeTriangleGroup(n=7)
            sage: lam = G.lam()
            sage: D = 4*lam^2 + 4*lam - 4
            sage: G.root_extension_embedding(D, CC)
            Relative number field morphism:
              From: Number Field in e with defining polynomial x^2 - 4*lam^2 - 4*lam + 4 over its base field
              To:   Complex Field with 53 bits of precision
              Defn: e |--> 4.02438434522...
                    lam |--> 1.80193773580...
            sage: D = lam^2 - 4
            sage: G.root_extension_embedding(D)
            Relative number field morphism:
              From: Number Field in e with defining polynomial x^2 - lam^2 + 4 over its base field
              To:   Algebraic Field
              Defn: e |--> 0.?... + 0.867767478235...?*I
                    lam |--> 1.801937735804...?
        """
        D = self.base_ring()(D)
        F = self.root_extension_field(D)
        if K is None:
            if coerce_AA(D) > 0:
                K = AA
            else:
                K = AlgebraicField()

        L = [emb for emb in F.embeddings(K)]

        # Three possibilities up to numerical artefacts:
        # (1) emb = e, purely imaginary
        # (2) emb = e or lam (can't distinguish), purely real
        # (3) emb = (e,lam), e purely imaginary, lam purely real
        # (4) emb = (e,lam), e purely real, lam purely real
        # There always exists one emb with "e" positive resp. positive imaginary
        # and if there is a lam there exists a positive one...
        #
        # Criteria to pick the correct "maximum":
        # 1. First figure out if e resp. lam is purely real or imaginary
        #    (using "abs(e.imag()) > abs(e.real())")
        # 2. In the purely imaginary case we don't want anything negative imaginary
        #    and we know the positive case is unique after sorting lam
        # 3. For the remaining cases we want the biggest real part
        #    (and lam should get comparison priority)
        def emb_key(emb):
            L = []
            gens_len = len(emb.im_gens())
            for k in range(gens_len):
                a = emb.im_gens()[k]
                try:
                    a.simplify()
                    a.exactify()
                except AttributeError:
                    pass
                # If a is purely imaginary:
                if abs(a.imag()) > abs(a.real()):
                    if a.imag() < 0:
                        a = -infinity
                    else:
                        a = ZZ(0)
                else:
                    a = a.real()

                L.append(a)

            L.reverse()
            return L

        if len(L) > 1:
            L.sort(key=emb_key)
        return L[-1]

    def _elliptic_conj_reps(self):
        r"""
        Store all elliptic conjugacy representatives in the
        internal dictionary. This is a helper function
        used by :meth:`_conjugacy_representatives`.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)
            sage: G.element_repr_method("conj")
            sage: G._elliptic_conj_reps()
            sage: sorted(G._conj_prim.items())
            [(-4, [[S], [S]]), (lam - 3, [[U], [U]]), (0, [[V(4)]])]
            sage: sorted(G._conj_nonprim.items())
            [(-lam - 2, [[U^(-2)], [U^2], [U^(-2)], [U^2]]), (lam - 3, [[U^(-1)], [U^(-1)]])]
        """
        if not hasattr(self, "_max_block_length"):
            self._conjugacy_representatives()
        elif ZZ(-4) in self._conj_prim:
            return

        D = self.U().discriminant()
        if D not in self._conj_prim:
            self._conj_prim[D] = []
        self._conj_prim[D].append(self.U())

        D = self.S().discriminant()
        if D not in self._conj_prim:
            self._conj_prim[D] = []
        self._conj_prim[D].append(self.S())

        other_reps = [self.U()**k
                      for k in range(-((self.n() - 1) / 2).floor(),
                                     self.n() // 2 + 1)
                      if k not in [0, 1]]

        for v in other_reps:
            D = v.discriminant()
            if D not in self._conj_nonprim:
                self._conj_nonprim[D] = []
            self._conj_nonprim[D].append(v)

    def _conjugacy_representatives(self, max_block_length=0, D=None):
        r"""
        Store conjugacy representatives up to block length
        ``max_block_length`` (a non-negative integer, default: 0)
        in the internal dictionary. Previously calculated data is reused.
        This is a helper function for e.g. :meth:`class_number`.

        The set of all (hyperbolic) conjugacy types of block length
        ``t`` is stored in ``self._conj_block[t]``.
        The set of all primitive representatives (so far) with
        discriminant ``D`` is stored in ``self._conj_prim[D]``.
        The set of all non-primitive representatives (so far) with
        discriminant ``D`` is stored in ``self._conj_nonprim[D]``.

        The case of non-positive discriminants is done manually.

        INPUT:

        - ``max_block_length`` -- A non-negative integer (default: ``0``),
                                  the maximal block length.

        - ``D``                -- An element/discriminant of the base ring or
                                  more generally an upper bound for the
                                  involved discriminants. If ``D != None``
                                  then an upper bound for ``max_block_length``
                                  is deduced from ``D`` (default: ``None``).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)
            sage: G.element_repr_method("conj")
            sage: G._conjugacy_representatives(2)

            sage: sorted(G._conj_block[2])
            [((2, 1), (1, 1)),
             ((2, 2),),
             ((3, 1), (1, 1)),
             ((3, 1), (2, 1)),
             ((3, 2),),
             ((4, 1), (1, 1)),
             ((4, 1), (2, 1)),
             ((4, 1), (3, 1))]

            sage: [key for key in sorted(G._conj_prim)]
            [-4, lam - 3, 0, 4*lam, 7*lam + 6, 9*lam + 5, 15*lam + 6, 33*lam + 21]
            sage: for key in sorted(G._conj_prim):
            ....:     print(sorted(G._conj_prim[key]))
            [[S], [S]]
            [[U], [U]]
            [[V(4)]]
            [[V(3)], [V(2)]]
            [[V(1)*V(4)]]
            [[V(3)*V(4)], [V(1)*V(2)]]
            [[V(2)*V(4)], [V(1)*V(3)]]
            [[V(2)*V(3)]]
            sage: [key for key in sorted(G._conj_nonprim)]
            [-lam - 2, lam - 3, 32*lam + 16]

            sage: for key in sorted(G._conj_nonprim):
            ....:     print(sorted(G._conj_nonprim[key]))
            [[U^(-2)], [U^(-2)], [U^2], [U^2]]
            [[U^(-1)], [U^(-1)]]
            [[V(3)^2], [V(2)^2]]

            sage: G.element_repr_method("default")
        """
        from sage.combinat.partition import OrderedPartitions
        from sage.combinat.combinat import tuples

        if D is not None:
            max_block_length = max(AA.zero(),
                                   coerce_AA((D + 4)/(self.lam()**2))).sqrt().floor()
        else:
            try:
                max_block_length = ZZ(max_block_length)
                if max_block_length < 0:
                    raise TypeError
            except TypeError:
                raise ValueError("max_block_length must be a non-negative integer!")

        if not hasattr(self, "_max_block_length"):
            self._max_block_length = ZZ.zero()
            self._conj_block = {}
            self._conj_nonprim = {}
            self._conj_prim = {}

            # It is not clear how to define the class number for D=0:
            # Conjugacy classes are V(n-1)^(+-k) for arbitrary k
            # and the trivial class (what about self_conj_block[0]?).
            #
            # One way is to define it using the fixed points and in
            # that case V(n-1) would be a good representative.
            # The non-primitive case is unclear however...
            #
            # We set it here to ensure that 0 is enlisted as a discriminant...
            #
            self._conj_prim[ZZ.zero()] = []
            self._conj_prim[ZZ.zero()].append(self.V(self.n()-1))

            self._elliptic_conj_reps()

        if max_block_length <= self._max_block_length:
            return

        def is_cycle(seq):
            length = len(seq)
            for n in divisors(length):
                if n < length and is_cycle_of_length(seq, n):
                    return True
            return False

        def is_cycle_of_length(seq, n):
            for j in range(n, len(seq)):
                if seq[j] != seq[j % n]:
                    return False
            return True

        j_list = range(1, self.n())

        for t in range(self._max_block_length + 1, max_block_length + 1):
            t_ZZ = ZZ(t)
            if t_ZZ not in self._conj_block:
                self._conj_block[t_ZZ] = set()

            partitions = OrderedPartitions(t).list()
            for par in partitions:
                len_par = len(par)
                exp_list = tuples(j_list, len_par)
                for ex in exp_list:
                    keep = True
                    if len_par > 1:
                        for k in range(-1, len_par - 1):
                            if ex[k] == ex[k + 1]:
                                keep = False
                                break
                    # We don't want powers of V(1)
                    elif ex[0] == 1:
                        keep = False
                    # But: Do we maybe want powers of V(n-1)??
                    # For now we exclude the parabolic cases...
                    elif ex[0] == self.n() - 1:
                        keep = False

                    if keep:
                        conj_type = cyclic_representative(tuple((ZZ(ex[k]), ZZ(par[k])) for k in range(len_par)))
                        self._conj_block[t_ZZ].add(conj_type)

            for el in self._conj_block[t_ZZ]:
                group_el = prod([self.V(el[k][0])**el[k][1] for k in range(len(el))])

                D = group_el.discriminant()
                assert coerce_AA(D) > 0

                # The primitive cases
                # if group_el.is_primitive():
                if not ((len(el) == 1 and el[0][1] > 1) or is_cycle(el)):
                    if D not in self._conj_prim:
                        self._conj_prim[D] = []
                    self._conj_prim[D].append(group_el)
                # The remaining cases
                else:
                    if D not in self._conj_nonprim:
                        self._conj_nonprim[D] = []
                    self._conj_nonprim[D].append(group_el)

        self._max_block_length = max_block_length

    def class_representatives(self, D, primitive=True):
        r"""
        Return a representative for each conjugacy class for the
        discriminant ``D`` (ignoring the sign).

        If ``primitive=True`` only one representative for each
        fixed point is returned (ignoring sign).

        INPUT:

        - ``D``          -- An element of the base ring corresponding
          to a valid discriminant.

        - ``primitive``  -- If ``True`` (default) then only primitive
          representatives are considered.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)
            sage: G.element_repr_method("conj")

            sage: R = G.class_representatives(4)
            sage: R
            [[V(2)]]
            sage: [v.continued_fraction()[1] for v in R]
            [(2,)]

            sage: R = G.class_representatives(0)
            sage: R
            [[V(3)]]
            sage: [v.continued_fraction()[1] for v in R]
            [(1, 2)]

            sage: R = G.class_representatives(-4)
            sage: R
            [[S]]
            sage: R = G.class_representatives(-4, primitive=False)
            sage: R
            [[S], [U^2]]

            sage: R = G.class_representatives(G.lam()^2 - 4)
            sage: R
            [[U]]
            sage: R = G.class_representatives(G.lam()^2 - 4, primitive=False)
            sage: R
            [[U], [U^(-1)]]

            sage: R = G.class_representatives(14)
            sage: sorted(R)
            [[V(2)*V(3)], [V(1)*V(2)]]
            sage: sorted(v.continued_fraction()[1] for v in R)
            [(1, 2, 2), (3,)]

            sage: R = G.class_representatives(32)
            sage: sorted(R)
            [[V(3)^2*V(1)], [V(1)^2*V(3)]]
            sage: [v.continued_fraction()[1] for v in sorted(R)]
            [(1, 2, 1, 3), (1, 4)]

            sage: R = G.class_representatives(32, primitive=False)
            sage: sorted(R)
            [[V(3)^2*V(1)], [V(1)^2*V(3)], [V(2)^2]]

            sage: G.element_repr_method("default")
        """
        if coerce_AA(D) == 0 and not primitive:
            raise ValueError("There are infinitely many non-primitive conjugacy classes of discriminant 0.")

        self._conjugacy_representatives(D=D)

        L = []
        if D in self._conj_prim:
            L += self._conj_prim[D]
        if not primitive and D in self._conj_nonprim:
            L += self._conj_nonprim[D]

        if not L:
            raise ValueError("D = {} is not a{} discriminant for {}".format(D, " primitive" if primitive else "", self))
        else:
            return L

    def class_number(self, D, primitive=True):
        r"""
        Return the class number of ``self`` for the discriminant ``D``.

        This is the number of conjugacy classes of (primitive) elements
        of discriminant ``D``.

        Note: Due to the 1-1 correspondence with hyperbolic fixed
        points resp. hyperbolic binary quadratic forms
        this also gives the class number in those cases.

        INPUT:

        - ``D``          -- An element of the base ring corresponding
          to a valid discriminant.

        - ``primitive``  -- If ``True`` (default) then only primitive
          elements are considered.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)

            sage: G.class_number(4)
            1
            sage: G.class_number(4, primitive=False)
            1
            sage: G.class_number(14)
            2
            sage: G.class_number(32)
            2
            sage: G.class_number(32, primitive=False)
            3
            sage: G.class_number(68)
            4
        """

        if coerce_AA(D) <= 0:
            raise NotImplementedError

        self._conjugacy_representatives(D=D)

        num = ZZ(0)
        if D in self._conj_prim:
            num = len(self._conj_prim[D])
        if not primitive and D in self._conj_nonprim:
            num += len(self._conj_nonprim[D])

        if num == 0:
            raise ValueError("D = {} is not a{} discriminant for {}".format(D, " primitive" if primitive else "", self))
        else:
            return num

    def is_discriminant(self, D, primitive=True) -> bool:
        r"""
        Return whether ``D`` is a discriminant of an element of ``self``.

        Note: Checking that something isn't a discriminant takes much
        longer than checking for valid discriminants.

        INPUT:

        - ``D``          -- An element of the base ring.

        - ``primitive``  -- If ``True`` (default) then only primitive
          elements are considered.

        OUTPUT:

        ``True`` if ``D`` is a primitive discriminant (a discriminant of
        a primitive element) and ``False`` otherwise.
        If ``primitive=False`` then also non-primitive elements are considered.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)

            sage: G.is_discriminant(68)
            True
            sage: G.is_discriminant(196, primitive=False)    # long time
            True
            sage: G.is_discriminant(2)
            False
        """
        self._conjugacy_representatives(0)
        t_bound = max(AA.zero(),
                      coerce_AA((D + 4) / (self.lam()**2))).sqrt().floor()
        for t in range(self._max_block_length + 1, t_bound + 1):
            self._conjugacy_representatives(t)

            if D in self._conj_prim:
                return True
            if not primitive and D in self._conj_nonprim:
                return True

        if D in self._conj_prim:
            return True
        elif not primitive and D in self._conj_nonprim:
            return True
        else:
            return False

    def list_discriminants(self, D, primitive=True, hyperbolic=True, incomplete=False):
        r"""
        Return a list of all discriminants up to some upper bound ``D``.

        INPUT:

        - ``D``          -- An element/discriminant of the base ring or
          more generally an upper bound for the discriminant.

        - ``primitive``  -- If ``True`` (default) then only primitive
          discriminants are listed.

        - ``hyperbolic`` -- If ``True`` (default) then only positive
          discriminants are listed.

        - ``incomplete`` -- If ``True`` (default: ``False``) then all
          (also higher) discriminants which were gathered so far are listed
          (however there might be missing discriminants inbetween).

        OUTPUT:

        A list of discriminants less than or equal to ``D``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)
            sage: G.list_discriminants(D=68)
            [4, 12, 14, 28, 32, 46, 60, 68]
            sage: G.list_discriminants(D=0, hyperbolic=False, primitive=False)
            [-4, -2, 0]

            sage: G = HeckeTriangleGroup(n=5)
            sage: G.list_discriminants(D=20)
            [4*lam, 7*lam + 6, 9*lam + 5]
            sage: G.list_discriminants(D=0, hyperbolic=False, primitive=False)
            [-4, -lam - 2, lam - 3, 0]
        """

        self._conjugacy_representatives(D=D)

        if incomplete:
            max_D = infinity
        else:
            max_D = coerce_AA(D)

        if hyperbolic:
            L = [key for key in self._conj_prim if 0 < coerce_AA(key) <= max_D]
        else:
            L = [key for key in self._conj_prim if coerce_AA(key) <= max_D]

        if not primitive:
            if hyperbolic:
                L += [key for key in self._conj_nonprim
                      if 0 < coerce_AA(key) <= max_D and key not in L]
            else:
                L += [key for key in self._conj_nonprim
                      if coerce_AA(key) <= max_D and key not in L]

        return sorted(L, key=coerce_AA)

    # TODO: non-primitive ones?
    def reduced_elements(self, D):
        r"""
        Return all reduced (primitive) elements of discriminant ``D``.

        Also see the element method :meth:`is_reduced` for more information.

        - ``D`` -- An element of the base ring corresponding
          to a valid discriminant.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)
            sage: R = G.reduced_elements(D=12)
            sage: R
            [
            [    5  -lam]  [     5 -3*lam]
            [3*lam    -1], [   lam     -1]
            ]
            sage: [v.continued_fraction() for v in R]
            [((), (1, 3)), ((), (3, 1))]
            sage: R = G.reduced_elements(D=14)
            sage: sorted(R)
            [
            [3*lam    -1]  [4*lam    -3]  [ 5*lam     -7]  [ 5*lam     -3]
            [    1     0], [    3  -lam], [     3 -2*lam], [     7 -2*lam]
            ]
            sage: sorted(v.continued_fraction() for v in R)
            [((), (1, 2, 2)), ((), (2, 1, 2)), ((), (2, 2, 1)), ((), (3,))]
        """
        L = self.class_representatives(D=D, primitive=True)
        R = []
        for v in L:
            R += v.reduced_elements()

        return R

    def simple_elements(self, D):
        r"""
        Return all simple elements of discriminant ``D``.

        Also see the element method ``is_simple()`` for more information.

        - ``D`` -- An element of the base ring corresponding
                   to a valid discriminant.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)
            sage: sorted(G.simple_elements(D=12))
            [
            [  1 lam]  [  3 lam]
            [lam   3], [lam   1]
            ]
            sage: sorted(G.simple_elements(D=14))
            [
            [  lam     1]  [  lam     3]  [2*lam     1]  [2*lam     3]
            [    3 2*lam], [    1 2*lam], [    3   lam], [    1   lam]
            ]
        """
        L = self.class_representatives(D=D, primitive=True)
        R = []
        for v in L:
            R += v.simple_elements()

        return R

    def rational_period_functions(self, k, D):
        r"""
        Return a list of basic rational period functions of weight ``k`` for discriminant ``D``.

        The list is expected to be a generating set for all rational
        period functions of the given weight and discriminant (unknown).

        The method assumes that ``D > 0``.

        Also see the element method `rational_period_function` for more information.

        - ``k`` -- An even integer, the desired weight
          of the rational period functions.

        - ``D`` -- An element of the base ring corresponding
          to a valid discriminant.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=4)
            sage: sorted(G.rational_period_functions(k=4, D=12))
            [(z^4 - 1)/z^4]
            sage: sorted(G.rational_period_functions(k=-2, D=12))
            [-z^2 + 1, 4*lam*z^2 - 4*lam]
            sage: sorted(G.rational_period_functions(k=2, D=14))
            [(24*z^6 - 120*z^4 + 120*z^2 - 24)/(9*z^8 - 80*z^6 + 146*z^4 - 80*z^2 + 9),
             (24*z^6 - 120*z^4 + 120*z^2 - 24)/(9*z^8 - 80*z^6 + 146*z^4 - 80*z^2 + 9),
             1/z,
             (z^2 - 1)/z^2]
            sage: sorted(G.rational_period_functions(k=-4, D=14))
            [-16*z^4 + 16, -z^4 + 1, 16*z^4 - 16]
        """
        try:
            k = ZZ(k)
            if not ZZ(2).divides(k):
                raise TypeError
        except TypeError:
            raise ValueError("k={} has to be an even integer!".format(k))

        z = PolynomialRing(self.base_ring(), 'z').gen()

        R = []
        if k != 0:
            R.append(ZZ(1) - z**(-k))
        if k == 2:
            R.append(z**(-1))

        L = self.class_representatives(D=D, primitive=True)
        for v in L:
            rat = v.rational_period_function(k)
            if rat != 0:
                R.append(rat)

        return R
