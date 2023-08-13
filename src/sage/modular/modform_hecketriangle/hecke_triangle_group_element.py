r"""
Hecke triangle group elements

AUTHORS:

- Jonas Jermann (2014): initial version

"""

# ****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.rings.cc import CC

lazy_import('sage.rings.qqbar', 'AA')

from sage.groups.matrix_gps.group_element import MatrixGroupElement_generic
from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane


# We want to simplify p after the coercion (pari bug for AA)
def coerce_AA(p):
    r"""
    Return the argument first coerced into ``AA`` and then simplified.

    This leads to a major performance gain with some operations.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.hecke_triangle_group_element import coerce_AA
        sage: p = (791264*AA(2*cos(pi/8))^2 - 463492).sqrt()
        sage: AA(p)._exact_field()
        Number Field in a with defining polynomial y^8 ... with a in ...
        sage: coerce_AA(p)._exact_field()
        Number Field in a with defining polynomial y^4 - 1910*y^2 - 3924*y + 681058 with a in ...?
    """
    el = AA(p)
    el.simplify()
    #el.exactify()

    return el


def cyclic_representative(L):
    r"""
    Return a unique representative among all cyclic permutations
    of the given list/tuple.

    INPUT:

    - ``L`` -- A list or tuple.

    OUTPUT:

    The maximal element among all cyclic permutations with respect
    to lexicographical ordering.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.hecke_triangle_group_element import cyclic_representative
        sage: cyclic_representative((1,))
        (1,)
        sage: cyclic_representative((2,2))
        (2, 2)
        sage: cyclic_representative((1,2,1,2))
        (2, 1, 2, 1)
        sage: cyclic_representative((1,2,3,2,3,1))
        (3, 2, 3, 1, 1, 2)
    """
    L = list(L)
    n = len(L)
    Lmax = L[:]
    for _ in range(n - 1):
        L.append(L.pop(0))
        if L > Lmax:
            Lmax = L[:]

    return tuple(Lmax)


class HeckeTriangleGroupElement(MatrixGroupElement_generic):
    r"""
    Elements of HeckeTriangleGroup.
    """
    def __init__(self, parent, M, check=True, **kwargs):
        r"""
        An element of HeckeTriangle group given by a matrix ``M``.

        INPUT:

        - ``parent`` -- A ``HeckeTriangleGroup``.

        - ``M``      -- A matrix which coerces into the matrix space
                        of ``parent``. For example with entries in a
                        polynomial ring over ``ZZ`` with parameter ``lam``.

        - ``check``  -- ``True`` (default) or ``False``. If ``True``
                        then a (possibly long) check is performed
                        to see whether ``M`` really corresponds to a
                        group element of ``parent``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup, HeckeTriangleGroupElement
            sage: lam = PolynomialRing(ZZ, 'lam').gen()
            sage: M = matrix([[-1, 0], [-lam^4 + 5*lam^2 + lam - 5, -1]])
            sage: G = HeckeTriangleGroup(4)
            sage: G(M)
            Traceback (most recent call last):
            ...
            TypeError: The matrix is not an element of Hecke triangle group for n = 4, up to equivalence it identifies two nonequivalent points.

            sage: G = HeckeTriangleGroup(10)
            sage: el = G(M)
            sage: el == HeckeTriangleGroupElement(G, M)
            True
            sage: type(el)
            <class 'sage.modular.modform_hecketriangle.hecke_triangle_groups.HeckeTriangleGroup_with_category.element_class'>
            sage: el.category()
            Category of elements of Hecke triangle group for n = 10
            sage: type(HeckeTriangleGroupElement(G, M))
            <class 'sage.modular.modform_hecketriangle.hecke_triangle_group_element.HeckeTriangleGroupElement'>
            sage: HeckeTriangleGroupElement(G, M).category()
            Category of elements of Hecke triangle group for n = 10
            sage: el
            [ -1   0]
            [lam  -1]
            sage: el.matrix().parent()
            Full MatrixSpace of 2 by 2 dense matrices over Maximal Order in Number Field in lam with defining polynomial x^4 - 5*x^2 + 5 with lam = 1.902113032590308?

            sage: M = matrix([[-1, lam], [0, 1]])
            sage: G(M)
            Traceback (most recent call last):
            ...
            TypeError: The matrix is not an element of Hecke triangle group for n = 10, it has determinant -1 != 1.

            sage: G.T().inverse()
            [   1 -lam]
            [   0    1]
            sage: G.U() == G.T()*G.S()
            True
            sage: G.U()^(-10) == -G.I()
            True
        """
        MatrixGroupElement_generic.__init__(self, parent, M, check=check, convert=True)

        # The matrix check involves a lengthy element method (_word_S_T_data)
        # whose result is also used for other purposes. For performance reason the
        # results are stored/cached in the element. Moreover this avoids code duplication.
        # In particular this means we cannot call the method from _matrix_check().
        # Instead it is called here in the __init__ method of the element
        # (after the preliminary checks).
        if check:
            if self._matrix.determinant() != 1:
                raise TypeError("The matrix is not an element of {}, it has determinant {} != 1.".format(parent, self._matrix.determinant()))
            self._word_S_T_data()

    @cached_method
    def _word_S_T_data(self):
        r"""
        Return a tuple ``(L, sgn)`` which describes the decomposition
        of ``self`` as a product of the generators ``S`` and ``T``
        together with a sign correction ``sgn``.

        If this decomposition is not possible a ``TypeError``
        is raised. In particular this function can be used to
        check the membership in ``parent`` of an arbitrary matrix
        over the base ring.

        OUTPUT:

        The tuple entries of ``L`` are either of the form ``(0, 1)``,
        corresponding to ``S`` or ``(1, m)`` corresponding to
        ``T^m``, where ``m`` is a non-trivial integer. ``sgn`` is +-1.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=17)
            sage: G.I()._word_S_T_data()
            ((), 1)
            sage: (-G.V(2))._word_S_T_data()
            (((1, 1), (0, 1), (1, 1)), -1)
            sage: G.U()._word_S_T_data()
            (((1, 1), (0, 1)), 1)

            sage: G = HeckeTriangleGroup(n=infinity)
            sage: (-G.V(2)*G.V(3))._word_S_T_data()
            (((1, 1), (0, 1), (1, 2), (0, 1), (1, 1), (0, 1), (1, 1)), -1)
            sage: G.U()._word_S_T_data()
            (((1, 1), (0, 1)), 1)
        """
        res = []
        ID = self.parent().I()._matrix
        T = self.parent().T()._matrix
        S = self.parent().S()._matrix
        M = self._matrix
        lam = self.parent().lam()
        zero = ZZ.zero()
        one = ZZ.one()
        half = one / ZZ(2)

        while True:
            a,b,c,d = M.list()
            mshift = coerce_AA((4*a*c + b*d) / (4*c*c + d*d))
            m = (mshift / lam + half).floor()
            if m != zero:
                res.append((one, m),)
                M = T**(-m) * M
                a,b,c,d = M.list()

            abs_t = coerce_AA((4*a*a + b*b) / (4*c*c + d*d))
            if coerce_AA(abs_t) < 1:
                M = (-S) * M
                res.append((zero, one),)
            elif M == ID:
                return (tuple(res), one)
            elif M == -ID:
                return (tuple(res), -one)
            else:
                raise TypeError("The matrix is not an element of {}, up to equivalence it identifies two nonequivalent points.".format(self.parent()))

    def word_S_T(self):
        r"""
        Decompose ``self`` into a product of the generators
        ``S`` and ``T`` of its parent, together with a sign
        correction matrix, namely: ``self = sgn * prod(L)``.

        Warning:
        If ``self`` is +- the identity ``prod(L)`` is an empty product
        which produces ``1`` instead of the identity matrix.

        OUTPUT:

        The function returns a tuple ``(L, sgn)`` where the entries
        of ``L`` are either the generator ``S`` or a non-trivial
        integer power of the generator ``T``. ``sgn`` is +- the identity.

        If this decomposition is not possible a ``TypeError`` is raised.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=17)
            sage: (-G.I()).word_S_T()[0]
            ()
            sage: (-G.I()).word_S_T()[1]
            [-1  0]
            [ 0 -1]
            sage: (L, sgn) = (-G.V(2)).word_S_T()
            sage: L
            (
            [  1 lam]  [ 0 -1]  [  1 lam]
            [  0   1], [ 1  0], [  0   1]
            )
            sage: sgn == -G.I()
            True
            sage: -G.V(2) == sgn * prod(L)
            True
            sage: (L, sgn) = G.U().word_S_T()
            sage: L
            (
            [  1 lam]  [ 0 -1]
            [  0   1], [ 1  0]
            )
            sage: sgn == G.I()
            True
            sage: G.U() == sgn * prod(L)
            True

            sage: G = HeckeTriangleGroup(n=infinity)
            sage: (L, sgn) = (-G.V(2)*G.V(3)).word_S_T()
            sage: L
            (
            [1 2]  [ 0 -1]  [1 4]  [ 0 -1]  [1 2]  [ 0 -1]  [1 2]
            [0 1], [ 1  0], [0 1], [ 1  0], [0 1], [ 1  0], [0 1]
            )
            sage: -G.V(2)*G.V(3) == sgn * prod(L)
            True
        """
        Tf = self.parent().T
        S = self.parent().S()
        (L, sgn) = self._word_S_T_data()

        M = [S if v[0]==0 else Tf(v[1]) for v in L]
        if sgn > 0:
            sgn = self.parent().I()
        else:
            sgn = -self.parent().I()

        return (tuple(M), sgn)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.
        The result depends on the default element representation
        method of the parent: ``self.parent().element_repr_method()``.

        See :meth:`string_repr` for a list of possible methods
        and for more examples.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)
            sage: el = G.S()*G.T(3)*G.S()*G.T(-2)
            sage: el
            [        -1      2*lam]
            [     3*lam -6*lam - 7]
            sage: el.string_repr()
            '[        -1      2*lam]\n[     3*lam -6*lam - 7]'
            sage: G.element_repr_method("basic")
            sage: el
            S*T^3*S*T^(-2)
            sage: el.string_repr("basic")
            'S*T^3*S*T^(-2)'
        """
        return self.string_repr(self.parent().element_repr_method())

    def string_repr(self, method="default"):
        r"""
        Return a string representation of ``self`` using the specified ``method``.
        This method is used to represent ``self``.
        The default representation method can be set for the parent with
        ``self.parent().element_repr_method(method)``.

        INPUT:

        - ``method``  -- ``default``: Use the usual representation method for matrix group elements.

                         ``basic``:   The representation is given as a word in ``S`` and powers of ``T``.
                                      Note: If ``S, T`` are defined accordingly the output can
                                      be used/evaluated directly to recover ``self``.

                         ``conj``:    The conjugacy representative of the element is represented
                                      as a word in powers of the basic blocks, together with
                                      an unspecified conjugation matrix.

                         ``block``:   Same as ``conj`` but the conjugation matrix is specified as well.
                                      Note: Assuming ``S, T, U, V`` are defined accordingly the output
                                      can directly be used/evaluated to recover ``self``.

        Warning: For ``n=infinity`` the methods ``conj`` and ``block`` are not verified at all
        and are probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)
            sage: el1 = -G.I()
            sage: el2 = G.S()*G.T(3)*G.S()*G.T(-2)
            sage: el3 = G.V(2)*G.V(3)^2*G.V(4)^3
            sage: el4 = G.U()^4
            sage: el5 = (G.V(2)*G.T()).acton(-G.S())

            sage: el4.string_repr(method="basic")
            'S*T^(-1)'

            sage: G.element_repr_method("default")
            sage: el1
            [-1  0]
            [ 0 -1]
            sage: el2
            [        -1      2*lam]
            [     3*lam -6*lam - 7]
            sage: el3
            [34*lam + 19   5*lam + 4]
            [27*lam + 18   5*lam + 2]
            sage: el4
            [   0   -1]
            [   1 -lam]
            sage: el5
            [-7*lam - 4  9*lam + 6]
            [-4*lam - 5  7*lam + 4]

            sage: G.element_repr_method("basic")
            sage: el1
            -1
            sage: el2
            S*T^3*S*T^(-2)
            sage: el3
            -T*S*T*S*T^(-1)*S*T^(-2)*S*T^(-4)*S
            sage: el4
            S*T^(-1)
            sage: el5
            T*S*T^2*S*T^(-2)*S*T^(-1)

            sage: G.element_repr_method("conj")
            sage: el1
            [-1]
            sage: el2
            [-V(4)^2*V(1)^3]
            sage: el3
            [V(3)^2*V(4)^3*V(2)]
            sage: el4
            [-U^(-1)]
            sage: el5
            [-S]

            sage: G.element_repr_method("block")
            sage: el1
            -1
            sage: el2
            -(S*T^3) * (V(4)^2*V(1)^3) * (S*T^3)^(-1)
            sage: el3
            (T*S*T) * (V(3)^2*V(4)^3*V(2)) * (T*S*T)^(-1)
            sage: el4
            -U^(-1)
            sage: el5
            -(T*S*T^2) * (S) * (T*S*T^2)^(-1)

            sage: G.element_repr_method("default")

            sage: G = HeckeTriangleGroup(n=infinity)
            sage: el = G.S()*G.T(3)*G.S()*G.T(-2)
            sage: print(el.string_repr())
            [ -1   4]
            [  6 -25]
            sage: print(el.string_repr(method="basic"))
            S*T^3*S*T^(-2)
        """
        if method == "default":
            return MatrixGroupElement_generic._repr_(self)
        elif method == "basic":
            (L, sgn) = self._word_S_T_data()

            if not L:
                return "-1" if sgn < 0 else "1"

            Lstr = list(L)
            for i,(v0,v1) in enumerate(Lstr):
                if v0 == 0:
                    Lstr[i] = "S"
                elif v1 == 1:
                    Lstr[i] = "T"
                else:
                    if v1 < 0:
                        exp = "(" + str(v1) + ")"
                    else:
                        exp = str(v1)
                    Lstr[i] = "T^" + exp
            Lstr = "*".join(Lstr)

            return "-" + Lstr if sgn < 0 else Lstr

        elif method == "block":
            if self.parent().n() == infinity:
                from warnings import warn
                warn("The case n=infinity here is not verified at all and probably wrong!")

            (L, R, sgn) = self._block_decomposition_data()

            repr_str = self.string_repr(method="conj")
            repr_str = repr_str[1:-1]
            if sgn < 0:
                repr_str = repr_str[1:]

            #if self != R.inverse().acton(self):
            if R.is_identity():
                repr_str = "{}{}".format("-" if sgn < 0 else "", repr_str)
            else:
                R_str = "({})".format(R.string_repr(method="basic"))
                repr_str = "{}{} * ({}) * {}^(-1)".format("-" if sgn < 0 else "", R_str, repr_str, R_str)

            return repr_str

        elif method == "conj":
            if self.parent().n() == infinity:
                from warnings import warn
                warn("The case n=infinity here is not verified at all and probably wrong!")

            (L, R, sgn) = self._block_decomposition_data()

            if self.is_elliptic():
                L = [ L ]

            repr_str = ""
            begin = True
            for v in L:
                if self.is_identity():
                    pass
                elif self.is_elliptic():
                    if v[0] == 0:
                        repr_str += "S"
                    elif v[1] == 1:
                        repr_str += "U"
                    else:
                        exp = "{}".format(v[1])
                        if v[1] < 0:
                            exp = "({})".format(exp)
                        repr_str += "U^{}".format(exp)
                    begin = False
                elif v[1] == 0:
                    pass
                else:
                    if not begin:
                        repr_str += "*"
                    factor = "V({})".format(v[0])
                    if v[1] == 1:
                        repr_str += factor
                    else:
                        exp = "{}".format(v[1])
                        if v[1] < 0:
                            exp = "({})".format(exp)
                        repr_str += "{}^{}".format(factor, exp)
                    begin = False

            if begin:
                repr_str += "1"

            if sgn < 0:
                repr_str = "-{}".format(repr_str)

            repr_str = "[{}]".format(repr_str)

            return repr_str
        else:
            raise NotImplementedError

    # We cache this method since the calculation is rather long and the
    # result is being reused:
    # - For block decompositions
    # - For calculating reduced (and simple) elements (and all dependent methods)
    @cached_method
    def continued_fraction(self):
        r"""
        For hyperbolic and parabolic elements: Return the (negative)
        lambda-continued fraction expansion (lambda-CF) of the (attracting)
        hyperbolic fixed point of ``self``.

        Let ``r_j in Z`` for ``j >= 0``. A finite lambda-CF is defined as:
        ``[r_0; r_1, ..., r_k] := (T^(r_0)*S* ... *T^(r_k)*S)(infinity)``,
        where ``S`` and ``T`` are the generators of ``self``. An infinite
        lambda-CF is defined as a corresponding limit value (k->infinity)
        if it exists.

        In this case the lambda-CF of parabolic and hyperbolic fixed points
        are returned which have an eventually periodic lambda-CF.
        The parabolic elements are exactly those with a cyclic permutation
        of the period ``[2, 1, ..., 1]`` with ``n-3`` ones.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        OUTPUT:

        A tuple ``(preperiod, period)`` with the preperiod and period
        tuples of the lambda-CF.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.T().continued_fraction()
            ((0, 1), (1, 1, 1, 1, 2))
            sage: G.V(2).acton(G.T(-3)).continued_fraction()
            ((), (2, 1, 1, 1, 1))
            sage: (-G.V(2)).continued_fraction()
            ((1,), (2,))
            sage: (-G.V(2)^3*G.V(6)^2*G.V(3)).continued_fraction()
            ((1,), (2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2))
            sage: (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3)).continued_fraction()
            ((1, 1, 1, 2), (2, 2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1))
            sage: (G.V(1)^5*G.V(2)*G.V(3)^3).continued_fraction()
            ((6,), (2, 1, 2, 1, 2, 1, 7))

            sage: G = HeckeTriangleGroup(n=8)
            sage: G.T().continued_fraction()
            ((0, 1), (1, 1, 1, 1, 1, 2))
            sage: G.V(2).acton(G.T(-3)).continued_fraction()
            ((), (2, 1, 1, 1, 1, 1))
            sage: (-G.V(2)).continued_fraction()
            ((1,), (2,))
            sage: (-G.V(2)^3*G.V(6)^2*G.V(3)).continued_fraction()
            ((1,), (2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2))
            sage: (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3)).continued_fraction()
            ((1, 1, 1, 2), (2, 2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1))
            sage: (G.V(1)^5*G.V(2)*G.V(3)^3).continued_fraction()
            ((6,), (2, 1, 2, 1, 2, 1, 7))
            sage: (G.V(2)^3*G.V(5)*G.V(1)*G.V(6)^2*G.V(4)).continued_fraction()
            ((1,), (2, 2, 2, 1, 1, 1, 3, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2))
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if self.is_identity():
            raise NotImplementedError
        if self.is_elliptic():
            # Note: The algorithm still produces "something"
            # emb = self.root_extension_embedding(QQbar)
            raise NotImplementedError

        emb = self.root_extension_embedding(AA)
        G = self.parent()
        S = G.S()
        TI = G.T().inverse()
        lam = G.lam()

        p = self.fixed_points()[0]

        cf_dict = {}
        L = []
        cf_index = ZZ.zero()
        one = ZZ.one()

        while p not in cf_dict:
            cf_dict[p] = cf_index
            if p == infinity:
                # TODO: The choice of r doesn't matter?
                r = ZZ.zero()
            #elif self.is_elliptic():
            #    r = ZZ(emb(p/lam).real().floor() + 1)
            else:
                emb_res = emb(p/lam)
                emb_res.simplify()
                emb_res.exactify()
                r = emb_res.floor() + one
            L.append(r)
            p = (S*TI**r).acton(p)
            cf_index += one

        preperiod_len = cf_dict[p]
        #period_len = cf_index - preperiod_len

        return (tuple(L[:preperiod_len]), tuple(L[preperiod_len:]))

    # TODO: allow output as word?
    # We cache this method since the calculation is rather long and the
    # data is being reused when working with primitive representatives
    # and conjugacy classes.
    @cached_method
    def _primitive_block_decomposition_data(self):
        r"""
        Return a tuple ``(L, R)`` which describes the
        decomposition of ``self`` into a very specific
        primitive conjugacy representative whose
        decomposition is further described by the tuple ``L``,
        and the corresponding  conjugation matrix ``R``.

        Together they describe the primitive part of self.
        I.e. an element which is equal to ``self`` up
        to a sign after taking the appropriate power
        and which itself cannot be written as a non-trivial
        power (at least for non-elliptic elements).

        To construct the representative see
        :meth:`primitive_representative`. To construct
        the primitive part see :meth:`primitive_part`.
        To get a corresponding decomposition of ``self``
        see :meth:`block_decomposition`.

        In the hyperbolic and parabolic case the
        representative is given as a product of powers of
        ``V(j)`` (more precisely ``self.parent().V(j)``),
        where ``1 <= j <= n-1``.

        The number of such factors is called ``block length``
        (see :meth:`block_length`). Each block (and also
        their product) has a positive sign and
        non-negative entries.

        In the elliptic case the primitive representative
        is either ``S`` or ``U``.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        OUTPUT:

        A tuple ``(L, R)``, where ``R`` is an element of
        the Hecke triangle group that conjugates the
        described primitive representative to the primitive
        part of ``self``.

        In the hyperbolic and parabolic case ``L`` is an
        ordered tuple of (tuple) data ``(j, k)``, corresponding
        to a factor ``V(j)^k``.

        If the representative is the identity then ``((1,0),)``
        is returned (consistent with the previous notation).

        In the elliptic case ``L=(a, 1)``, with either ``a=0``
        corresponding to the representative ``S`` or ``a=1``
        corresponding to the representative ``U``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.element_repr_method("basic")

            sage: (L, R) = G.T()._primitive_block_decomposition_data()
            sage: L
            ((1, 1),)
            sage: R
            T^(-1)
            sage: (L, R) = G.V(2).acton(G.T(-3))._primitive_block_decomposition_data()
            sage: L
            ((6, 1),)
            sage: R
            T
            sage: (L, R) = (-G.V(2))._primitive_block_decomposition_data()
            sage: L
            ((2, 1),)
            sage: R
            T*S*T
            sage: (L, R) = (-G.V(2)^3*G.V(6)^2*G.V(3))._primitive_block_decomposition_data()
            sage: L
            ((2, 3), (6, 2), (3, 1))
            sage: R
            1
            sage: (L, R) = (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3))._primitive_block_decomposition_data()
            sage: L
            ((2, 3), (6, 2), (3, 1))
            sage: R
            T*S*T*S*T*S*T^2*S*T
            sage: (L, R) = (G.V(1)^5*G.V(2)*G.V(3)^3)._primitive_block_decomposition_data()
            sage: L
            ((3, 3), (1, 5), (2, 1))
            sage: R
            T^6*S*T

            sage: G.element_repr_method("default")
            sage: (L, R) = G.I()._primitive_block_decomposition_data()
            sage: L
            ((6, 0),)
            sage: R
            [1 0]
            [0 1]

            sage: (L, R) = G.U()._primitive_block_decomposition_data()
            sage: L
            (1, 1)
            sage: R
            [1 0]
            [0 1]
            sage: (L, R) = (-G.S())._primitive_block_decomposition_data()
            sage: L
            (0, 1)
            sage: R
            [-1  0]
            [ 0 -1]
            sage: (L, R) = (G.V(2)*G.V(3)).acton(G.U()^6)._primitive_block_decomposition_data()
            sage: L
            (1, 1)
            sage: R
            [-2*lam^2 - 2*lam + 2 -2*lam^2 - 2*lam + 1]
            [        -2*lam^2 + 1   -2*lam^2 - lam + 2]
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        G = self.parent()
        zero = ZZ.zero()
        one = ZZ.one()

        # The elliptic case (for this case we use a special notation):
        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            from sage.rings.qqbar import QQbar

            emb = self.root_extension_embedding(QQbar)
            p = self.fixed_points()[0]
            embp = emb(p)
            embp.simplify()
            embp.exactify()
            (R, embw) = G.get_FD(embp)
            w = R.inverse().acton(p)
            # we should have: embw == emb(w)
            embw = emb(w)
            embw.simplify()
            embw.exactify()

            if embw == QQbar.gen():
                R = -R
                L = (zero, one)
            elif (embw == -one/G.rho()):
                R = R*G.T().inverse()
                L = (one, one)
            else:
                raise RuntimeError("There is something wrong in the method "
                         "_primitive_block_decomposition_data. Please contact sage-devel@googlegroups.com")

            return (L, R)

        # The identity case (consistent with the notation in the parabolic case):
        if self.is_identity():
            return (((ZZ(self.parent().n()-one), zero),), G.I())

        # The hyperbolic and parabolic case:
        # The parabolic case is much simpler but the same algorithm
        # as in the hyperbolic case still works

        (preperiod, period) = self.continued_fraction()

        number_of_ones = []
        list_larger = []
        ones = 0
        for l in period:
            if l > 1:
                number_of_ones.append(ones)
                ones = 0
                list_larger.append(l)
            else:
                ones += 1
        number_of_ones.append(ones)

        initial_ones = number_of_ones.pop(0)
        if not list_larger:
            list_v1 = [-ZZ(1)]
            list_vlarger = [ initial_ones + 2 ]
        else:
            list_v1 = [ v-2 for v in list_larger ]
            list_vlarger = [ v+2 for v in number_of_ones ]
            list_vlarger[-1] += initial_ones

        L = []
        for k in range(len(list_vlarger)):
            if list_v1[k] != 0:
                L.append([ZZ(1), list_v1[k]])
            L.append([ZZ(list_vlarger[k]), ZZ(1)])

        L_len = len(L)
        k = 0
        while k < L_len - 1:
            if L[k][0] == L[k+1][0]:
                k_entry = L.pop(k+1)
                L[k][1] += k_entry[1]
                L_len -= 1
            else:
                k += 1
        if L_len > 1 and L[-1][0] == L[0][0]:
            k_entry = L.pop(-1)
            L[0][1] += k_entry[1]
            R = G.V(L[0][0])**(-k_entry[1])
        else:
            R = G.I()

        # This should determine whether self is conjugate to a positive power of V(1) or V(n-1)
        # sign((a+d)*(b-c)) is actually a conjugacy invariant for the parabolic subspace
        # and distinguishes the two (three) cases (sign(0):=0)
        if self.is_parabolic() and coerce_AA(self.trace() * (self.b() - self.c())).sign() > 0:
            # In this case self should be conjugate to a positive power of V(1)
            # in either case L is / should be (at the moment) always equal to [n-1, 1]
            L[0][0] = 1
            R = R * (-G.S())

        R = G.V(initial_ones + 1) * R
        R = prod((G.T(r) * G.S() for r in preperiod), G.I()) * R

        L = tuple(tuple(v) for v in L)

        return (L, R)

    def primitive_representative(self, method="block"):
        r"""
        Return a tuple ``(P, R)`` which gives the
        decomposition of the primitive part of ``self``,
        namely ``R*P*R.inverse()`` into a specific
        representative ``P`` and the corresponding
        conjugation matrix ``R`` (the result depends on
        the method used).

        Together they describe the primitive part of self.
        I.e. an element which is equal to ``self`` up
        to a sign after taking the appropriate power.

        See :meth:`_primitive_block_decomposition_data` for a description
        about the representative in case the default method
        ``block`` is used. Also see :meth:`primitive_part`
        to construct the primitive part of ``self``.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``method`` -- ``block`` (default) or ``cf``. The method
                        used to determine ``P`` and ``R``. If
                        ``self`` is elliptic this parameter is
                        ignored and if ``self`` is +- the identity
                        then the ``block`` method is used.

                        With ``block`` the decomposition described
                        in :meth:`_primitive_block_decomposition_data` is used.

                        With ``cf`` a reduced representative from
                        the lambda-CF of ``self`` is used (see
                        :meth:`continued_fraction`). In that case
                        ``P`` corresponds to the period and ``R``
                        to the preperiod.

        OUTPUT:

        A tuple ``(P, R)`` of group elements such that
        ``R*P*R.inverse()`` is a/the primitive part of ``self``

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.element_repr_method("basic")
            sage: el = G.T().primitive_representative(method="cf")
            sage: el
            (S*T^(-1)*S*T^(-1)*S*T*S, S*T*S)
            sage: (el[0]).is_primitive()
            True
            sage: el = G.V(2).acton(G.T(-3)).primitive_representative(method="cf")
            sage: el
            (-T*S*T^(-1)*S*T^(-1), 1)
            sage: (el[0]).is_primitive()
            True
            sage: el = (-G.V(2)).primitive_representative(method="cf")
            sage: el
            (T^2*S, T*S)
            sage: (el[0]).is_primitive()
            True
            sage: el = (-G.V(2)^3*G.V(6)^2*G.V(3)).primitive_representative(method="cf")
            sage: el
            (-T^2*S*T^2*S*T*S*T^(-2)*S*T*S*T*S*T^2*S, T*S)
            sage: (el[0]).is_primitive()
            True
            sage: el = (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3)).primitive_representative(method="cf")
            sage: el
            (-T^2*S*T^2*S*T^2*S*T*S*T^(-2)*S*T*S*T*S, T*S*T*S*T*S*T^2*S)
            sage: (el[0]).is_primitive()
            True
            sage: el = (G.V(1)^5*G.V(2)*G.V(3)^3).primitive_representative(method="cf")
            sage: el
            (T^2*S*T*S*T^2*S*T*S*T^2*S*T*S*T^7*S, T^6*S)
            sage: (el[0]).is_primitive()
            True
            sage: el = (G.V(2)*G.V(3)).acton(G.U()^6).primitive_representative(method="cf")
            sage: el
            (T*S, -T*S*T^2*S*T*S*T)
            sage: (el[0]).is_primitive()
            True

            sage: G.element_repr_method("block")
            sage: el = G.T().primitive_representative()
            sage: (el[0]).is_primitive()
            True
            sage: el = G.V(2).acton(G.T(-3)).primitive_representative()
            sage: el
            ((-S*T^(-1)*S) * (V(6)) * (-S*T^(-1)*S)^(-1), (T^(-1)) * (V(1)) * (T^(-1))^(-1))
            sage: (el[0]).is_primitive()
            True
            sage: el = (-G.V(2)).primitive_representative()
            sage: el
            ((T*S*T) * (V(2)) * (T*S*T)^(-1), (T*S*T) * (V(2)) * (T*S*T)^(-1))
            sage: (el[0]).is_primitive()
            True
            sage: el = (-G.V(2)^3*G.V(6)^2*G.V(3)).primitive_representative()
            sage: el
            (V(2)^3*V(6)^2*V(3), 1)
            sage: (el[0]).is_primitive()
            True
            sage: el = (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3)).primitive_representative()
            sage: el
            (V(2)^3*V(6)^2*V(3), (T*S*T*S*T*S*T) * (V(2)*V(4)) * (T*S*T*S*T*S*T)^(-1))
            sage: (el[0]).is_primitive()
            True
            sage: el = (G.V(1)^5*G.V(2)*G.V(3)^3).primitive_representative()
            sage: el
            (V(3)^3*V(1)^5*V(2), (T^6*S*T) * (V(1)^5*V(2)) * (T^6*S*T)^(-1))
            sage: (el[0]).is_primitive()
            True

            sage: G.element_repr_method("default")
            sage: el = G.I().primitive_representative()
            sage: el
            (
            [1 0]  [1 0]
            [0 1], [0 1]
            )
            sage: (el[0]).is_primitive()
            True

            sage: el = G.U().primitive_representative()
            sage: el
            (
            [lam  -1]  [1 0]
            [  1   0], [0 1]
            )
            sage: (el[0]).is_primitive()
            True
            sage: el = (-G.S()).primitive_representative()
            sage: el
            (
            [ 0 -1]  [-1  0]
            [ 1  0], [ 0 -1]
            )
            sage: (el[0]).is_primitive()
            True
            sage: el = (G.V(2)*G.V(3)).acton(G.U()^6).primitive_representative()
            sage: el
            (
            [lam  -1]  [-2*lam^2 - 2*lam + 2 -2*lam^2 - 2*lam + 1]
            [  1   0], [        -2*lam^2 + 1   -2*lam^2 - lam + 2]
            )
            sage: (el[0]).is_primitive()
            True
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        G = self.parent()

        if self.is_identity():
            method="block"

        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            (data, R) = self._primitive_block_decomposition_data()
            if data[0] == 0:
                P = G.S()
            else:
                P = G.U()

            return (P, R)

        if method == "cf":
            (preperiod, period) = self.continued_fraction()

            P = prod((G.T()**r*G.S() for r in period), G.I())
            R = prod((G.T()**r*G.S() for r in preperiod), G.I())

            return (P, R)

        elif method == "block":
            (data_list, R) = self._primitive_block_decomposition_data()
            P = prod((G.V(v[0])**v[1] for v in data_list), G.I())

            return (P, R)

        else:
            raise ValueError("if the element is not elliptic, then method must "
                             "be either be 'cf' or 'block'")

    def primitive_part(self, method="cf"):
        r"""
        Return the primitive part of ``self``. I.e. a group element
        ``A`` with non-negative trace such that
        ``self = sign * A^power``, where ``sign = self.sign()``
        is +- the identity (to correct the sign) and
        ``power = self.primitive_power()``.

        The primitive part itself is choosen such that it cannot be
        written as a non-trivial power of another element.
        It is a generator of the stabilizer of the corresponding
        (attracting) fixed point.

        If ``self`` is elliptic then the primitive part is
        chosen as a conjugate of ``S`` or ``U``.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``method``  -- The method used to determine the primitive
                         part (see :meth:`primitive_representative`),
                         default: "cf". The parameter is ignored
                         for elliptic elements or +- the identity.

                         The result should not depend on the method.

        OUTPUT:

        The primitive part as a group element of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.element_repr_method("block")
            sage: G.T().primitive_part()
            (T^(-1)*S) * (V(6)) * (T^(-1)*S)^(-1)
            sage: G.V(2).acton(G.T(-3)).primitive_part()
            (T) * (V(6)) * (T)^(-1)
            sage: (-G.V(2)).primitive_part()
            (T*S*T) * (V(2)) * (T*S*T)^(-1)
            sage: (-G.V(2)^3*G.V(6)^2*G.V(3)).primitive_part()
            V(2)^3*V(6)^2*V(3)
            sage: (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3)).primitive_part()
            (T*S*T*S*T*S*T^2*S*T) * (V(2)^3*V(6)^2*V(3)) * (T*S*T*S*T*S*T^2*S*T)^(-1)
            sage: (G.V(1)^5*G.V(2)*G.V(3)^3).primitive_part()
            (T^6*S*T) * (V(3)^3*V(1)^5*V(2)) * (T^6*S*T)^(-1)
            sage: (G.V(2)*G.V(3)).acton(G.U()^6).primitive_part()
            (-T*S*T^2*S*T*S*T) * (U) * (-T*S*T^2*S*T*S*T)^(-1)

            sage: (-G.I()).primitive_part()
            1

            sage: G.U().primitive_part()
            U
            sage: (-G.S()).primitive_part()
            S
            sage: el = (G.V(2)*G.V(3)).acton(G.U()^6)
            sage: el.primitive_part()
            (-T*S*T^2*S*T*S*T) * (U) * (-T*S*T^2*S*T*S*T)^(-1)
            sage: el.primitive_part() == el.primitive_part(method="block")
            True

            sage: G.T().primitive_part()
            (T^(-1)*S) * (V(6)) * (T^(-1)*S)^(-1)
            sage: G.T().primitive_part(method="block")
            (T^(-1)) * (V(1)) * (T^(-1))^(-1)
            sage: G.V(2).acton(G.T(-3)).primitive_part() == G.V(2).acton(G.T(-3)).primitive_part(method="block")
            True
            sage: (-G.V(2)).primitive_part() == (-G.V(2)).primitive_part(method="block")
            True
            sage: el = -G.V(2)^3*G.V(6)^2*G.V(3)
            sage: el.primitive_part() == el.primitive_part(method="block")
            True
            sage: el = (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3))
            sage: el.primitive_part() == el.primitive_part(method="block")
            True
            sage: el=G.V(1)^5*G.V(2)*G.V(3)^3
            sage: el.primitive_part() == el.primitive_part(method="block")
            True

            sage: G.element_repr_method("default")
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        (P, R) = self.primitive_representative(method=method)

        return R*P*R.inverse()

    def reduce(self, primitive=True):
        r"""
        Return a reduced version of ``self`` (with the same
        the same fixed points). Also see :meth:`is_reduced`.

        If ``self`` is elliptic (or +- the identity) the result
        is never reduced (by definition). Instead a more canonical
        conjugation representative of ``self`` (resp. it's
        primitive part) is choosen.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``primitive`` -- If ``True`` (default) then a primitive
                           representative for ``self`` is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: print(G.T().reduce().string_repr("basic"))
            S*T^(-1)*S*T^(-1)*S*T*S
            sage: G.T().reduce().is_reduced(require_hyperbolic=False)
            True
            sage: print(G.V(2).acton(-G.T(-3)).reduce().string_repr("basic"))
            -T*S*T^(-1)*S*T^(-1)
            sage: print(G.V(2).acton(-G.T(-3)).reduce(primitive=False).string_repr("basic"))
            T*S*T^(-3)*S*T^(-1)
            sage: print((-G.V(2)).reduce().string_repr("basic"))
            T^2*S
            sage: (-G.V(2)).reduce().is_reduced()
            True
            sage: print((-G.V(2)^3*G.V(6)^2*G.V(3)).reduce().string_repr("block"))
            (-S*T^(-1)) * (V(2)^3*V(6)^2*V(3)) * (-S*T^(-1))^(-1)
            sage: (-G.V(2)^3*G.V(6)^2*G.V(3)).reduce().is_reduced()
            True

            sage: print((-G.I()).reduce().string_repr("block"))
            1
            sage: print(G.U().reduce().string_repr("block"))
            U
            sage: print((-G.S()).reduce().string_repr("block"))
            S
            sage: print((G.V(2)*G.V(3)).acton(G.U()^6).reduce().string_repr("block"))
            U
            sage: print((G.V(2)*G.V(3)).acton(G.U()^6).reduce(primitive=False).string_repr("block"))
            -U^(-1)
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        (P, R) = self.primitive_representative(method="cf")

        if primitive:
            return P
        else:
            return R.inverse().acton(self)

    def sign(self):
        r"""
        Return the sign element/matrix (+- identity) of ``self``.
        The sign is given by the sign of the trace.
        if the trace is zero it is instead given by the sign
        of the lower left entry.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: (-G.T(-1)).sign()
            [-1  0]
            [ 0 -1]
            sage: G.S().sign()
            [1 0]
            [0 1]
            sage: (-G.S()).sign()
            [-1  0]
            [ 0 -1]
            sage: (G.U()^6).sign()
            [-1  0]
            [ 0 -1]

            sage: G = HeckeTriangleGroup(n=8)
            sage: (G.U()^4).trace()
            0
            sage: (G.U()^4).sign()
            [1 0]
            [0 1]
            sage: (G.U()^(-4)).sign()
            [-1  0]
            [ 0 -1]
        """
        sgn = coerce_AA(self._matrix.trace()).sign()

        if sgn > 0:
            return self.parent().I()
        elif sgn < 0:
            return -self.parent().I()
        else:
            sgnc = coerce_AA(self.c()).sign()
            if sgnc > 0:
                return self.parent().I()
            elif sgnc < 0:
                return -self.parent().I()
            else:
                raise AssertionError("This shouldn't happen!")

    @cached_method
    def primitive_power(self, method="cf"):
        r"""
        Return the primitive power of ``self``. I.e. an integer
        ``power`` such that ``self = sign * primitive_part^power``,
        where ``sign = self.sign()`` and
        ``primitive_part = self.primitive_part(method)``.

        Warning: For the parabolic case the sign depends on
        the method: The "cf" method may return a negative power
        but the "block" method never will.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``method``  -- The method used to determine the primitive
                         power (see :meth:`primitive_representative`),
                         default: "cf". The parameter is ignored
                         for elliptic elements or +- the identity.

        OUTPUT:

        An integer. For +- the identity element ``0`` is returned,
        for parabolic and hyperbolic elements a positive integer.
        And for elliptic elements a (non-zero) integer with minimal
        absolute value such that ``primitive_part^power`` still
        has a positive sign.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.T().primitive_power()
            -1
            sage: G.V(2).acton(G.T(-3)).primitive_power()
            3
            sage: (-G.V(2)^2).primitive_power()
            2
            sage: el = (-G.V(2)*G.V(6)*G.V(3)*G.V(2)*G.V(6)*G.V(3))
            sage: el.primitive_power()
            2
            sage: (G.U()^4*G.S()*G.V(2)).acton(el).primitive_power()
            2
            sage: (G.V(2)*G.V(3)).acton(G.U()^6).primitive_power()
            -1
            sage: G.V(2).acton(G.T(-3)).primitive_power() == G.V(2).acton(G.T(-3)).primitive_power(method="block")
            True

            sage: (-G.I()).primitive_power()
            0
            sage: G.U().primitive_power()
            1
            sage: (-G.S()).primitive_power()
            1
            sage: el = (G.V(2)*G.V(3)).acton(G.U()^6)
            sage: el.primitive_power()
            -1
            sage: el.primitive_power() == (-el).primitive_power()
            True
            sage: (G.U()^(-6)).primitive_power()
            1

            sage: G = HeckeTriangleGroup(n=8)
            sage: (G.U()^4).primitive_power()
            4
            sage: (G.U()^(-4)).primitive_power()
            4
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        zero = ZZ.zero()
        one = ZZ.one()
        two = ZZ(2)

        if self.is_identity():
            return zero

        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            (data, R) = self._primitive_block_decomposition_data()
            if data[0] == 0:
                return one
            else:
                G = self.parent()
                U = G.U()
                U_power = R.inverse() * self * R

                Uj = G.I()
                for j in range(1, G.n()):
                    Uj *= U
                    if U_power == Uj:
                        #L = [one, ZZ(j)]
                        break
                    elif U_power == -Uj:
                        #L = [one, ZZ(-j)]
                        break
                else:
                    raise RuntimeError("There is a problem in the method "
                    "'primitive_power'. Please contact sage-devel@googlegroups.com")

                if abs(j) < G.n()/two:
                    return j
                elif two*j == G.n():
                    return j
                # for the cases fom here on the sign has to be adjusted to the
                # sign of self (in self._block_decomposition_data())
                elif two*j == -G.n():
                    return -j
                elif j > 0:
                    return j - G.n()
                else:
                    return j + G.n()

        primitive_part = self.primitive_part(method=method)
        if method == "cf" and self.is_parabolic():
            power_sign = coerce_AA(self.trace() * (self[1][0] - self[0][1])).sign()
        else:
            power_sign = one

        normalized_self = self.sign() * self**power_sign
        M = primitive_part

        power = 1
        while M != normalized_self:
            M *= primitive_part
            power += 1

        return power*power_sign

    def block_length(self, primitive=False):
        r"""
        Return the block length of ``self``. The block length is
        given by the number of factors used for the decomposition
        of the conjugacy representative of ``self`` described in
        :meth:`primitive_representative`. In particular the block
        length is invariant under conjugation.

        The definition is mostly used for parabolic or hyperbolic
        elements: In particular it gives a lower bound for the
        (absolute value of) the trace and the discriminant for
        primitive hyperbolic elements. Namely
        ``abs(trace) >= lambda * block_length`` and
        ``discriminant >= block_length^2 * lambda^2 - 4``.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``primitive``  -- If ``True`` then the conjugacy
                         representative of the primitive part is
                         used instead, default: ``False``.

        OUTPUT:

        An integer. For hyperbolic elements a non-negative integer.
        For parabolic elements a negative sign corresponds to taking
        the inverse. For elliptic elements a (non-trivial) integer
        with minimal absolute value is choosen. For +- the identity
        element ``0`` is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.T().block_length()
            1
            sage: G.V(2).acton(G.T(-3)).block_length()
            3
            sage: G.V(2).acton(G.T(-3)).block_length(primitive=True)
            1
            sage: (-G.V(2)).block_length()
            1

            sage: el = -G.V(2)^3*G.V(6)^2*G.V(3)
            sage: t = el.block_length()
            sage: D = el.discriminant()
            sage: trace = el.trace()
            sage: (trace, D, t)
            (-124*lam^2 - 103*lam + 68, 65417*lam^2 + 52456*lam - 36300, 6)
            sage: abs(AA(trace)) >= AA(G.lam()*t)
            True
            sage: AA(D) >= AA(t^2 * G.lam() - 4)
            True
            sage: (el^3).block_length(primitive=True) == t
            True

            sage: el = (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3))
            sage: t = el.block_length()
            sage: D = el.discriminant()
            sage: trace = el.trace()
            sage: (trace, D, t)
            (-124*lam^2 - 103*lam + 68, 65417*lam^2 + 52456*lam - 36300, 6)
            sage: abs(AA(trace)) >= AA(G.lam()*t)
            True
            sage: AA(D) >= AA(t^2 * G.lam() - 4)
            True
            sage: (el^(-2)).block_length(primitive=True) == t
            True

            sage: el = G.V(1)^5*G.V(2)*G.V(3)^3
            sage: t = el.block_length()
            sage: D = el.discriminant()
            sage: trace = el.trace()
            sage: (trace, D, t)
            (284*lam^2 + 224*lam - 156, 330768*lam^2 + 265232*lam - 183556, 9)
            sage: abs(AA(trace)) >= AA(G.lam()*t)
            True
            sage: AA(D) >= AA(t^2 * G.lam() - 4)
            True
            sage: (el^(-1)).block_length(primitive=True) == t
            True

            sage: (G.V(2)*G.V(3)).acton(G.U()^6).block_length()
            1
            sage: (G.V(2)*G.V(3)).acton(G.U()^6).block_length(primitive=True)
            1

            sage: (-G.I()).block_length()
            0
            sage: G.U().block_length()
            1
            sage: (-G.S()).block_length()
            1
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if primitive:
            L = self._primitive_block_decomposition_data()[0]
        else:
            L = self._block_decomposition_data()[0]

        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            return abs(L[1])
        else:
            return sum(abs(v[1]) for v in L)

    #@cached_method
    def _block_decomposition_data(self):
        r"""
        Return a tuple ``(L, R, sgn)`` which describes the
        decomposition of ``self`` into a specific
        conjugacy representative whose decomposition is
        further described by the tuple ``L``. The conjugation
        matrix is returned as ``R`` and since all
        factors have a positive sign, the sign ``sgn``
        of ``self`` is supplied as well as +- 1 (which corresponds
        to the sign of the sign matrix ``self.sign()``).

        The function is a generalization of
        :meth:`_primitive_block_decomposition_data`
        (see for more information) to give the decomposition data
        for possibly non-primitive elements.

        Also see :meth:`block_decomposition()` for more information
        on the block decomposition.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        OUTPUT:

        A tuple ``(L, R, sgn)``, where ``R`` is an element of
        the Hecke triangle group that conjugates the
        described representative to ``self`` up to the given sign.

        In the hyperbolic and parabolic case ``L`` is an
        ordered tuple of (tuple) data ``(j, k)``, corresponding
        to a factor ``V(j)^k``.

        If the representative is the identity then ``((1, 0))``
        is returned (consistent with the previous notation).

        In the elliptic case ``L=(a, k)``, with either ``a=0``
        corresponding to the representative ``S^k`` or with
        ``a=1`` corresponding to the representative ``U^k``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.element_repr_method("basic")

            sage: (L, R, sgn) = G.T()._block_decomposition_data()
            sage: (L, sgn)
            (((1, 1),), 1)
            sage: R
            T^(-1)
            sage: (L, R, sgn) = G.V(2).acton(G.T(-3))._block_decomposition_data()
            sage: (L, sgn)
            (((6, 3),), 1)
            sage: R
            T
            sage: (L, R, sgn) = (-G.V(2)^2)._block_decomposition_data()
            sage: (L, sgn)
            (((2, 2),), -1)
            sage: R
            T*S*T
            sage: el = (-G.V(2)*G.V(6)*G.V(3)*G.V(2)*G.V(6)*G.V(3))
            sage: (L, R, sgn) = el._block_decomposition_data()
            sage: (L, sgn)
            (((6, 1), (3, 1), (2, 1), (6, 1), (3, 1), (2, 1)), -1)
            sage: R
            T*S*T
            sage: (L, R, sgn) = (G.U()^4*G.S()*G.V(2)).acton(el)._block_decomposition_data()
            sage: (L, sgn)
            (((2, 1), (6, 1), (3, 1), (2, 1), (6, 1), (3, 1)), -1)
            sage: R
            T*S*T*S*T*S*T^2*S*T
            sage: (L, R, sgn) = (G.V(1)^5*G.V(2)*G.V(3)^3)._block_decomposition_data()
            sage: (L, sgn)
            (((3, 3), (1, 5), (2, 1)), 1)
            sage: R
            T^6*S*T

            sage: G.element_repr_method("default")
            sage: (L, R, sgn) = (-G.I())._block_decomposition_data()
            sage: (L, sgn)
            (((6, 0),), -1)
            sage: R
            [1 0]
            [0 1]
            sage: (L, R, sgn) = G.U()._block_decomposition_data()
            sage: (L, sgn)
            ((1, 1), 1)
            sage: R
            [1 0]
            [0 1]
            sage: (L, R, sgn) = (-G.S())._block_decomposition_data()
            sage: (L, sgn)
            ((0, 1), -1)
            sage: R
            [-1  0]
            [ 0 -1]
            sage: (L, R, sgn) = (G.V(2)*G.V(3)).acton(G.U()^6)._block_decomposition_data()
            sage: (L, sgn)
            ((1, -1), -1)
            sage: R
            [-2*lam^2 - 2*lam + 2 -2*lam^2 - 2*lam + 1]
            [        -2*lam^2 + 1   -2*lam^2 - lam + 2]
            sage: (L, R, sgn) = (G.U()^(-6))._block_decomposition_data()
            sage: (L, sgn)
            ((1, 1), -1)
            sage: R
            [1 0]
            [0 1]

            sage: G = HeckeTriangleGroup(n=8)
            sage: (L, R, sgn) = (G.U()^4)._block_decomposition_data()
            sage: (L, sgn)
            ((1, 4), 1)
            sage: R
            [1 0]
            [0 1]
            sage: (L, R, sgn) = (G.U()^(-4))._block_decomposition_data()
            sage: (L, sgn)
            ((1, 4), -1)
            sage: R
            [1 0]
            [0 1]
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        (L, R) = self._primitive_block_decomposition_data()
        if self.sign() == self.parent().I():
            sgn = ZZ(1)
        else:
            sgn = ZZ(-1)

        if self.is_identity():
            return (L, R, sgn)

        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            # Since L is primitive L[1] should be equal to 1
            M = (L[0], self.primitive_power())
            return (M, R, sgn)

        # If the length of L is 1, there is always at most one block
        # Note that this is includes the parabolic case
        if len(L) == 1:
            # Since L is primitive L[0][1] should be equal to 1
            # Also note that in the non-elliptic case:
            # abs(self.primitive_power()) == self.primitive_power(method="block")
            L2 = ((L[0][0], abs(self.primitive_power())),)
        else:
            L2 = L*abs(self.primitive_power())

        return (L2, R, sgn)

    def block_decomposition(self):
        r"""
        Return a tuple ``(L, R, sgn)`` such that
        ``self = sgn * R.acton(prod(L)) = sgn * R*prod(L)*R.inverse()``.

        In the parabolic and hyperbolic case the tuple entries
        in ``L`` are powers of basic block matrices:
        ``V(j) = U^(j-1)*T = self.parent().V(j)`` for ``1 <= j <= n-1``.
        In the elliptic case the tuple entries are either ``S`` or ``U``.

        This decomposition data is (also) described by
        :meth:`_block_decomposition_data`.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.element_repr_method("basic")

            sage: G.T().block_decomposition()
            ((T,), T^(-1), 1)
            sage: G.V(2).acton(G.T(-3)).block_decomposition()
            ((-S*T^(-3)*S,), T, 1)
            sage: (-G.V(2)^2).block_decomposition()
            ((T*S*T^2*S*T,), T*S*T, -1)

            sage: el = (-G.V(2)*G.V(6)*G.V(3)*G.V(2)*G.V(6)*G.V(3))
            sage: el.block_decomposition()
            ((-S*T^(-1)*S, T*S*T*S*T, T*S*T, -S*T^(-1)*S, T*S*T*S*T, T*S*T), T*S*T, -1)
            sage: (G.U()^4*G.S()*G.V(2)).acton(el).block_decomposition()
            ((T*S*T, -S*T^(-1)*S, T*S*T*S*T, T*S*T, -S*T^(-1)*S, T*S*T*S*T), T*S*T*S*T*S*T^2*S*T, -1)
            sage: (G.V(1)^5*G.V(2)*G.V(3)^3).block_decomposition()
            ((T*S*T*S*T^2*S*T*S*T^2*S*T*S*T, T^5, T*S*T), T^6*S*T, 1)

            sage: G.element_repr_method("default")
            sage: (-G.I()).block_decomposition()
            (
            ([1 0]   [1 0]  [-1  0]
            [0 1],), [0 1], [ 0 -1]
            )
            sage: G.U().block_decomposition()
            (
            ([lam  -1]   [1 0]  [1 0]
            [  1   0],), [0 1], [0 1]
            )
            sage: (-G.S()).block_decomposition()
            (
            ([ 0 -1]   [-1  0]  [-1  0]
            [ 1  0],), [ 0 -1], [ 0 -1]
            )
            sage: (G.V(2)*G.V(3)).acton(G.U()^6).block_decomposition()
            (
            ([  0   1]   [-2*lam^2 - 2*lam + 2 -2*lam^2 - 2*lam + 1]  [-1  0]
            [ -1 lam],), [        -2*lam^2 + 1   -2*lam^2 - lam + 2], [ 0 -1]
            )
            sage: (G.U()^(-6)).block_decomposition()
            (
            ([lam  -1]   [1 0]  [-1  0]
            [  1   0],), [0 1], [ 0 -1]
            )

            sage: G = HeckeTriangleGroup(n=8)
            sage: (G.U()^4).block_decomposition()
            (
            ([     lam^2 - 1 -lam^3 + 2*lam]   [1 0]  [1 0]
            [ lam^3 - 2*lam     -lam^2 + 1],), [0 1], [0 1]
            )
            sage: (G.U()^(-4)).block_decomposition()
            (
            ([     lam^2 - 1 -lam^3 + 2*lam]   [1 0]  [-1  0]
            [ lam^3 - 2*lam     -lam^2 + 1],), [0 1], [ 0 -1]
            )
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        G = self.parent()
        (L, R, sgn) = self._block_decomposition_data()
        if sgn > 0:
            sgn = G.I()
        else:
            sgn = -G.I()

        if self.is_identity():
            return ((G.I(),), R, sgn)

        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            if L[0] == 0:
                P = G.S()
            else:
                P = G.U()
            return ((P**L[1],), R, sgn)
        else:
            return (tuple(G.V(v[0])**v[1] for v in L ), R, sgn)

    def conjugacy_type(self, ignore_sign=True, primitive=False):
        r"""
        Return a unique description of the conjugacy class of ``self``
        (by default only up to a sign).

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``ignore_sign``  -- If ``True`` (default) then the conjugacy
                              classes are only considered up to a sign.

        - ``primitive``    -- If ``True`` then the conjugacy class of
                              the primitive part is considered instead
                              and the sign is ignored, default: ``False``.

        OUTPUT:

        A unique representative for the given block data (without the
        conjugation matrix) among all cyclic permutations.
        If ``ignore_sign=True`` then the sign is excluded as well.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: (-G.I()).conjugacy_type()
            ((6, 0),)
            sage: G.U().acton(G.S()).conjugacy_type()
            (0, 1)
            sage: (G.U()^4).conjugacy_type()
            (1, -3)
            sage: ((G.V(2)*G.V(3)^2*G.V(2)*G.V(3))^2).conjugacy_type()
            ((3, 2), (2, 1), (3, 1), (2, 1), (3, 2), (2, 1), (3, 1), (2, 1))

            sage: (-G.I()).conjugacy_type(ignore_sign=False)
            (((6, 0),), -1)
            sage: G.S().conjugacy_type(ignore_sign=False)
            ((0, 1), 1)
            sage: (G.U()^4).conjugacy_type(ignore_sign=False)
            ((1, -3), -1)
            sage: G.U().acton((G.V(2)*G.V(3)^2*G.V(2)*G.V(3))^2).conjugacy_type(ignore_sign=False)
            (((3, 2), (2, 1), (3, 1), (2, 1), (3, 2), (2, 1), (3, 1), (2, 1)), 1)

            sage: (-G.I()).conjugacy_type(primitive=True)
            ((6, 0),)
            sage: G.S().conjugacy_type(primitive=True)
            (0, 1)
            sage: G.V(2).acton(G.U()^4).conjugacy_type(primitive=True)
            (1, 1)
            sage: (G.V(3)^2).conjugacy_type(primitive=True)
            ((3, 1),)
            sage: G.S().acton((G.V(2)*G.V(3)^2*G.V(2)*G.V(3))^2).conjugacy_type(primitive=True)
            ((3, 2), (2, 1), (3, 1), (2, 1))
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if primitive:
            ignore_sign = True
            (L, R) = self._primitive_block_decomposition_data()
        else:
            (L, R, sgn) = self._block_decomposition_data()

        if not self.is_elliptic():
            L = tuple(cyclic_representative(L))
        return L if ignore_sign else (L, sgn)

    def reduced_elements(self):
        r"""
        Return the cycle of reduced elements in the (primitive)
        conjugacy class of ``self``.

        I.e. the set (cycle) of all reduced elements which are
        conjugate to ``self.primitive_part()``.
        E.g. ``self.primitive_representative().reduce()``.

        Also see :meth:`is_reduced`.
        In particular the result of this method only depends on the
        (primitive) conjugacy class of ``self``.

        The method assumes that ``self`` is hyperbolic or parabolic.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)
            sage: G.element_repr_method("basic")

            sage: el = G.V(1)
            sage: el.continued_fraction()
            ((0, 1), (1, 1, 2))
            sage: R = el.reduced_elements()
            sage: R
            [T*S*T*S*T^2*S, T*S*T^2*S*T*S, -T*S*T^(-1)*S*T^(-1)]
            sage: [v.continued_fraction() for v in R]
            [((), (1, 1, 2)), ((), (1, 2, 1)), ((), (2, 1, 1))]

            sage: el = G.V(3)*G.V(2)^(-1)*G.V(1)*G.V(6)
            sage: el.continued_fraction()
            ((1,), (3,))
            sage: R = el.reduced_elements()
            sage: [v.continued_fraction() for v in R]
            [((), (3,))]

            sage: G.element_repr_method("default")
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if self.is_identity() or self.is_elliptic():
            raise NotImplementedError

        def rotate(l, n):
            return tuple(l[n:] + l[:n])

        G = self.parent()
        L = []
        period = self.continued_fraction()[1]
        period_set = set()
        for k in range(len(period)):
            cur_period = rotate(period, k)
            if cur_period in period_set:
                continue
            else:
                period_set.add(cur_period)
                L.append(prod((G.T()**r*G.S() for r in cur_period), G.I()))

        return L

    def simple_elements(self):
        r"""
        Return all simple elements in the primitive conjugacy
        class of ``self``.

        I.e. the set of all simple elements which are
        conjugate to ``self.primitive_part()``.

        Also see :meth:`is_simple`.
        In particular the result of this method only depends on the
        (primitive) conjugacy class of ``self``.

        The method assumes that ``self`` is hyperbolic.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)

            sage: el = G.V(2)
            sage: el.continued_fraction()
            ((1,), (2,))
            sage: R = el.simple_elements()
            sage: R
            [
            [lam lam]
            [  1 lam]
            ]
            sage: R[0].is_simple()
            True

            sage: el = G.V(3)*G.V(2)^(-1)*G.V(1)*G.V(6)
            sage: el.continued_fraction()
            ((1,), (3,))
            sage: R = el.simple_elements()
            sage: R
            [
            [    2*lam 2*lam + 1]  [      lam 2*lam + 1]
            [        1       lam], [        1     2*lam]
            ]
            sage: [v.is_simple() for v in R]
            [True, True]

            sage: el = G.V(1)^2*G.V(2)*G.V(4)
            sage: el.discriminant()
            135*lam + 86
            sage: R = el.simple_elements()
            sage: R
            [
            [    3*lam 3*lam + 2]  [8*lam + 3 3*lam + 2]  [5*lam + 2 9*lam + 6]
            [3*lam + 4 6*lam + 3], [  lam + 2       lam], [  lam + 2 4*lam + 1],
            [2*lam + 1 7*lam + 4]
            [  lam + 2 7*lam + 2]
            ]

        This agrees with the results (p.16) from Culp-Ressler on
        binary quadratic forms for Hecke triangle groups::

            sage: [v.continued_fraction() for v in R]
            [((1,), (1, 1, 4, 2)),
            ((3,), (2, 1, 1, 4)),
            ((2,), (2, 1, 1, 4)),
            ((1,), (2, 1, 1, 4))]
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if not self.is_hyperbolic():
            return []

        G = self.parent()
        emb = self.root_extension_embedding(AA)
        R = self.reduced_elements()
        L = []

        for r in R:
            fp = r.fixed_points()[0]

            emb_res = emb(fp/G.lam())
            emb_res.simplify()
            emb_res.exactify()
            for j in range(1, emb_res.floor()+1):
                L.append(G.T(-j).acton(r))

        return L

    def simple_fixed_point_set(self, extended=True):
        r"""
        Return a set of all attracting fixed points in the
        conjugacy class of the primitive part of ``self``.

        If ``extended=True`` (default) then also
        ``S.acton(alpha)`` are added for ``alpha`` in the set.

        This is a so called `irreducible system of poles`
        for rational period functions for the parent group.
        I.e. the fixed points occur as a irreducible part
        of the non-zero pole set of some rational period
        function and all pole sets are given as a union
        of such irreducible systems of poles.

        The method assumes that ``self`` is hyperbolic.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)

            sage: el = G.V(2)
            sage: el.simple_fixed_point_set()
            {1/2*e, (-1/2*lam + 1/2)*e}
            sage: el.simple_fixed_point_set(extended=False)
            {1/2*e}

            sage: el = G.V(3)*G.V(2)^(-1)*G.V(1)*G.V(6)
            sage: el.simple_fixed_point_set()
            {(-lam + 3/2)*e + 1/2*lam - 1, (-lam + 3/2)*e - 1/2*lam + 1, 1/2*e - 1/2*lam, 1/2*e + 1/2*lam}

            sage: el.simple_fixed_point_set(extended=False)
            {1/2*e - 1/2*lam, 1/2*e + 1/2*lam}
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if self.is_identity() or self.is_elliptic():
            raise NotImplementedError

        from sage.sets.set import Set

        R = self.simple_elements()
        FPS = Set(v.fixed_points()[0] for v in R)

        if not extended:
            return FPS

        S = self.parent().S()
        FPS2 = Set(S.acton(v) for v in FPS)

        return FPS.union(FPS2)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: V = HeckeTriangleGroup(17).V(13)
            sage: latex(V)
            \begin{pmatrix} \mathit{\lambda}^{3} - 2 \mathit{\lambda} & \mathit{\lambda}^{2} - 1 \\ \mathit{\lambda}^{4} - 3 \mathit{\lambda}^{2} + 1 & \mathit{\lambda}^{3} - 2 \mathit{\lambda} \end{pmatrix}
        """
        latex_out = r"\begin{pmatrix} %s & %s \\ %s & %s \end{pmatrix}" % (latex(self.a()), latex(self.b()), latex(self.c()), latex(self.d()))
        return latex_out.replace("lam", r"\lambda")

    def __neg__(self):
        r"""
        Return the group element corresponding to the negative of the underlying matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: -U
            [-lam    1]
            [  -1    0]
        """
        return self.parent().element_class(self.parent(), -self._matrix, check=False)

    def __getitem__(self, key):
        r"""
        Return the corresponding rows/entries of the underlying matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U[0]
            (lam, -1)
            sage: U[0].parent()
            Ambient free module of rank 2 over the principal ideal domain Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1 with lam = 1.801937735804839?
            sage: U[1][0]
            1
            sage: U[1][0].parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1 with lam = 1.801937735804839?
        """
        return self._matrix[key]

    def a(self):
        r"""
        Return the upper left entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.a()
            lam
            sage: U.a().parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1 with lam = 1.801937735804839?
        """
        return self._matrix[0][0]

    def b(self):
        r"""
        Return the upper right entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.b()
            -1
            sage: U.b().parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1 with lam = 1.801937735804839?
        """
        return self._matrix[0][1]

    def c(self):
        r"""
        Return the lower left entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.c()
            1
        """
        return self._matrix[1][0]

    def d(self):
        r"""
        Return the lower right of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.d()
            0
        """
        return self._matrix[1][1]

    def trace(self):
        r"""
        Return the trace of ``self``, which is the sum of the diagonal entries.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.U().trace()
            lam
            sage: G.S().trace()
            0
        """
        return self._matrix.trace()

    def discriminant(self):
        r"""
        Return the discriminant of ``self`` which corresponds to
        the discriminant of the corresponding quadratic form of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.V(3).discriminant()
            4*lam^2 + 4*lam - 4
            sage: AA(G.V(3).discriminant())
            16.19566935808922?
        """

        return self.trace()**2 - 4

    def is_translation(self, exclude_one=False):
        r"""
        Return whether ``self`` is a translation. If ``exclude_one = True``,
        then the identity map is not considered as a translation.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: (-HeckeTriangleGroup(n=7).T(-4)).is_translation()
            True
            sage: (-HeckeTriangleGroup(n=7).I()).is_translation()
            True
            sage: (-HeckeTriangleGroup(n=7).I()).is_translation(exclude_one=True)
            False
        """
        a,b,c,d = self._matrix.list()

        if not (c.is_zero() and a == d and (a.is_one() or (-a).is_one())):
            return False
        elif exclude_one and b.is_zero():
            return False
        else:
            return True

    def is_reflection(self):
        r"""
        Return whether ``self`` is the usual reflection on the unit circle.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: (-HeckeTriangleGroup(n=7).S()).is_reflection()
            True
            sage: HeckeTriangleGroup(n=7).U().is_reflection()
            False
        """
        return self == self.parent().S() or self == -self.parent().S()

    def is_hyperbolic(self):
        r"""
        Return whether ``self`` is a hyperbolic matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_hyperbolic() for k in range(1,8) ]
            [False, True, True, True, True, False, False]
            sage: G.U().is_hyperbolic()
            False
        """
        return coerce_AA(self.discriminant()) > 0

    def is_parabolic(self, exclude_one=False):
        r"""
        Return whether ``self`` is a parabolic matrix.

        If ``exclude_one`` is set, then +- the identity
        element is not considered parabolic.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_parabolic() for k in range(1,8) ]
            [True, False, False, False, False, True, False]
            sage: G.U().is_parabolic()
            False
            sage: G.V(6).is_parabolic(exclude_one=True)
            True
            sage: G.V(7).is_parabolic(exclude_one=True)
            False
        """
        if exclude_one and self.is_identity():
            return False

        return self.discriminant() == 0

    def is_identity(self):
        r"""
        Return whether ``self`` is the identity or minus the identity.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_identity() for k in range(1,8) ]
            [False, False, False, False, False, False, False]
            sage: G.U().is_identity()
            False
        """
        if self == self.parent().I() or self == -self.parent().I():
            return True
        else:
            return False

    def is_elliptic(self):
        r"""
        Return whether ``self`` is an elliptic matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_elliptic() for k in range(1,8) ]
            [False, False, False, False, False, False, True]
            sage: G.U().is_elliptic()
            True
        """
        return coerce_AA(self.discriminant()) < 0

    def is_primitive(self):
        r"""
        Returns whether ``self`` is primitive. We call an element
        primitive if (up to a sign and taking inverses) it generates
        the full stabilizer subgroup of the corresponding fixed point.
        In the non-elliptic case this means that primitive elements
        cannot be written as a `non-trivial` power of another element.

        The notion is mostly used for hyperbolic and parabolic elements.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.V(2).acton(G.T(-1)).is_primitive()
            True
            sage: G.T(3).is_primitive()
            False
            sage: (-G.V(2)^2).is_primitive()
            False
            sage: (G.V(1)^5*G.V(2)*G.V(3)^3).is_primitive()
            True

            sage: (-G.I()).is_primitive()
            True
            sage: (-G.U()).is_primitive()
            True
            sage: (-G.S()).is_primitive()
            True
            sage: (G.U()^6).is_primitive()
            True

            sage: G = HeckeTriangleGroup(n=8)
            sage: (G.U()^2).is_primitive()
            False
            sage: (G.U()^(-4)).is_primitive()
            False
            sage: (G.U()^(-3)).is_primitive()
            True
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        pow = self.primitive_power()

        if self.is_elliptic():
            if self.parent().n() == infinity:
                raise NotImplementedError

            # if this is not up-to-sign then a factor 2 should
            # be added before (the second) self.parent().n()
            return (pow % (2*self.parent().n())).gcd(self.parent().n()) == 1
        else:
            return abs(pow) <= 1

    def is_reduced(self, require_primitive=True, require_hyperbolic=True):
        r"""
        Returns whether ``self`` is reduced. We call an element
        reduced if the associated lambda-CF is purely periodic.

        I.e. (in the hyperbolic case) if the associated hyperbolic
        fixed point (resp. the associated hyperbolic binary quadratic form)
        is reduced.

        Note that if ``self`` is reduced then the element corresponding
        to the cyclic permutation of the lambda-CF (which is conjugate
        to the original element) is again reduced. In particular the
        reduced elements in the conjugacy class of ``self`` form a
        finite cycle.

        Elliptic elements and +- identity are not considered reduced.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        INPUT:

        - ``require_primitive``  -- If ``True`` (default) then non-primitive elements
                                    are not considered reduced.

        - ``require_hyperbolic`` -- If ``True`` (default) then non-hyperbolic elements
                                    are not considered reduced.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.I().is_reduced(require_hyperbolic=False)
            False
            sage: G.U().reduce().is_reduced(require_hyperbolic=False)
            False
            sage: G.T().reduce().is_reduced()
            False
            sage: G.T().reduce().is_reduced(require_hyperbolic=False)
            True
            sage: (G.V(5)^2).reduce(primitive=False).is_reduced()
            False
            sage: (G.V(5)^2).reduce(primitive=False).is_reduced(require_primitive=False)
            True
            sage: G.V(5).reduce().is_reduced()
            True
            sage: (-G.V(2)).reduce().is_reduced()
            True
            sage: (-G.V(2)^3*G.V(6)^2*G.V(3)).reduce().is_reduced()
            True
            sage: (G.U()^4*G.S()*G.V(2)).acton(-G.V(2)^3*G.V(6)^2*G.V(3)).reduce().is_reduced()
            True
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if self.is_identity() or self.is_elliptic():
            return False
        elif require_hyperbolic and not self.is_hyperbolic():
            return False
        elif require_primitive and not self.is_primitive():
            return False
        else:
            return self.continued_fraction()[0] == ()

    def is_simple(self):
        r"""
        Return whether ``self`` is simple. We call an element
        simple if it is hyperbolic, primitive, has positive sign
        and if the associated hyperbolic fixed points satisfy:
        ``alpha' < 0 < alpha`` where ``alpha`` is the attracting
        fixed point for the element.

        I.e. if the associated hyperbolic fixed point (resp. the
        associated hyperbolic binary quadratic form) is simple.

        There are only finitely many simple elements for a given
        discriminant. They can be used to provide explicit
        descriptions of rational period functions.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)

            sage: el = G.V(2)
            sage: el.is_simple()
            True
            sage: R = el.simple_elements()
            sage: [v.is_simple() for v in R]
            [True]
            sage: (fp1, fp2) = R[0].fixed_points(embedded=True)
            sage: (fp1, fp2)
            (1.272019649514069?, -1.272019649514069?)
            sage: fp2 < 0 < fp1
            True

            sage: el = G.V(3)*G.V(2)^(-1)*G.V(1)*G.V(6)
            sage: el.is_simple()
            False
            sage: R = el.simple_elements()
            sage: [v.is_simple() for v in R]
            [True, True]
            sage: (fp1, fp2) = R[1].fixed_points(embedded=True)
            sage: fp2 < 0 < fp1
            True

            sage: el = G.V(1)^2*G.V(2)*G.V(4)
            sage: el.is_simple()
            True
            sage: R = el.simple_elements()
            sage: el in R
            True
            sage: [v.is_simple() for v in R]
            [True, True, True, True]
            sage: (fp1, fp2) = R[2].fixed_points(embedded=True)
            sage: fp2 < 0 < fp1
            True
        """
        if self != self.primitive_part():
            return False

        # The last condition is/should be equivalent to:
        a,b,c,d = self._matrix.list()
        return (coerce_AA(a) > 0 and coerce_AA(b) > 0 and coerce_AA(c) > 0 and coerce_AA(d) > 0)

    def is_hecke_symmetric(self):
        r"""
        Return whether the conjugacy class of the primitive part of
        ``self``, denoted by ``[gamma]`` is `Hecke-symmetric`:
        I.e. if ``[gamma] == [gamma^(-1)]``.

        This is equivalent to ``self.simple_fixed_point_set()`` being
        equal with it's `Hecke-conjugated` set (where each fixed point
        is replaced by the other (`Hecke-conjugated`) fixed point.

        It is also equivalent to ``[Q] == [-Q]`` for the corresponding
        hyperbolic binary quadratic form ``Q``.

        The method assumes that ``self`` is hyperbolic.

        .. WARNING::

            The case ``n=infinity`` is not verified at all and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)

            sage: el = G.V(2)
            sage: el.is_hecke_symmetric()
            False
            sage: (el.simple_fixed_point_set(), el.inverse().simple_fixed_point_set())
            ({1/2*e, (-1/2*lam + 1/2)*e}, {-1/2*e, (1/2*lam - 1/2)*e})

            sage: el = G.V(3)*G.V(2)^(-1)*G.V(1)*G.V(6)
            sage: el.is_hecke_symmetric()
            False
            sage: el.simple_fixed_point_set() == el.inverse().simple_fixed_point_set()
            False

            sage: el = G.V(2)*G.V(3)
            sage: el.is_hecke_symmetric()
            True
            sage: sorted(el.simple_fixed_point_set(), key=str)
            [(-lam + 3/2)*e + 1/2*lam - 1,
             (-lam + 3/2)*e - 1/2*lam + 1,
             (lam - 3/2)*e + 1/2*lam - 1,
             (lam - 3/2)*e - 1/2*lam + 1]
            sage: el.simple_fixed_point_set() == el.inverse().simple_fixed_point_set()
            True
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if self.is_identity() or self.is_elliptic():
            raise NotImplementedError

        return self.conjugacy_type() == self.inverse().conjugacy_type()

    def rational_period_function(self, k):
        r"""
        The method assumes that ``self`` is hyperbolic.

        Return the rational period function of weight ``k`` for
        the primitive conjugacy class of ``self``.

        A `rational period function` of weight ``k`` is a
        rational function ``q`` which satisfies:
        ``q + q|S == 0`` and ``q + q|U + q|U^2 + ... + q|U^(n-1) == 0``,
        where ``S = self.parent().S()``, ``U = self.parent().U()`` and
        ``|`` is the usual `slash-operator` of weight `k`.
        Note that if ``k < 0`` then ``q`` is a polynomial.

        This method returns a very basic rational period function
        associated with the primitive conjugacy class of ``self``.
        The (strong) expectation is that all rational period functions
        are formed by linear combinations of such functions.

        There is also a close relation with modular integrals of
        weight ``2-k`` and sometimes ``2-k`` is used for the weight
        instead of ``k``.

        Warning: The case ``n=infinity`` is not verified at all
        and probably wrong!

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=5)
            sage: S = G.S()
            sage: U = G.U()

            sage: def is_rpf(f, k=None):
            ....:     if not f + S.slash(f, k=k) == 0:
            ....:         return False
            ....:     if not sum([(U^m).slash(f, k=k) for m in range(G.n())]) == 0:
            ....:         return False
            ....:     return True

            sage: z = PolynomialRing(G.base_ring(), 'z').gen()
            sage: [is_rpf(1 - z^(-k), k=k) for k in range(-6, 6, 2)]  # long time
            [True, True, True, True, True, True]
            sage: [is_rpf(1/z, k=k) for k in range(-6, 6, 2)]
            [False, False, False, False, True, False]

            sage: el = G.V(2)
            sage: el.is_hecke_symmetric()
            False
            sage: rpf = el.rational_period_function(-4)
            sage: is_rpf(rpf) == is_rpf(rpf, k=-4)
            True
            sage: is_rpf(rpf)
            True
            sage: is_rpf(rpf, k=-6)
            False
            sage: is_rpf(rpf, k=2)
            False
            sage: rpf
            -lam*z^4 + lam
            sage: rpf = el.rational_period_function(-2)
            sage: is_rpf(rpf)
            True
            sage: rpf
            (lam + 1)*z^2 - lam - 1
            sage: el.rational_period_function(0) == 0
            True
            sage: rpf = el.rational_period_function(2)
            sage: is_rpf(rpf)
            True
            sage: rpf
            ((lam + 1)*z^2 - lam - 1)/(lam*z^4 + (-lam - 2)*z^2 + lam)

            sage: el = G.V(3)*G.V(2)^(-1)*G.V(1)*G.V(6)
            sage: el.is_hecke_symmetric()
            False
            sage: rpf = el.rational_period_function(-6)
            sage: is_rpf(rpf)
            True
            sage: rpf
            (68*lam + 44)*z^6 + (-24*lam - 12)*z^4 + (24*lam + 12)*z^2 - 68*lam - 44
            sage: rpf = el.rational_period_function(-2)
            sage: is_rpf(rpf)
            True
            sage: rpf
            (4*lam + 4)*z^2 - 4*lam - 4
            sage: el.rational_period_function(0) == 0
            True
            sage: rpf = el.rational_period_function(2)
            sage: is_rpf(rpf) == is_rpf(rpf, k=2)
            True
            sage: is_rpf(rpf)
            True
            sage: rpf.denominator()
            (8*lam + 5)*z^8 + (-94*lam - 58)*z^6 + (199*lam + 124)*z^4 + (-94*lam - 58)*z^2 + 8*lam + 5

            sage: el = G.V(2)*G.V(3)
            sage: el.is_hecke_symmetric()
            True
            sage: el.rational_period_function(-4) == 0
            True
            sage: rpf = el.rational_period_function(-2)
            sage: is_rpf(rpf)
            True
            sage: rpf
            (8*lam + 4)*z^2 - 8*lam - 4
            sage: el.rational_period_function(0) == 0
            True
            sage: rpf = el.rational_period_function(2)
            sage: is_rpf(rpf)
            True
            sage: rpf.denominator()
            (144*lam + 89)*z^8 + (-618*lam - 382)*z^6 + (951*lam + 588)*z^4 + (-618*lam - 382)*z^2 + 144*lam + 89
            sage: el.rational_period_function(4) == 0
            True
        """
        if self.parent().n() == infinity:
            from warnings import warn
            warn("The case n=infinity here is not verified at all and probably wrong!")

        if self.is_identity() or self.is_elliptic():
            raise NotImplementedError("This method is not implemented for the identity or elliptic element")

        try:
            k = ZZ(k)
            if k % 2:
                raise TypeError
        except TypeError:
            raise ValueError("k={} must be an even integer!".format(k))

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        P = PolynomialRing(self.parent().base_ring(), 'z')
        z = P.gen()

        s = P.zero()

        #L1 = []
        for v in self.simple_elements():
            a,b,c,d = v._matrix.list()
            Q = c*z**2 + (d - a)*z - b
            s += Q**(-k/ZZ(2))

        for v in self.inverse().simple_elements():
            a,b,c,d = v._matrix.list()
            Q = c*z**2 + (d - a)*z - b
            s -= ZZ(-1)**(k/ZZ(2)) * Q**(-k/ZZ(2))

        return s

    def linking_number(self):
        r"""
        Let ``g`` denote a holomorphic primitive of ``E2`` in the sense:
        ``lambda/(2*pi*i) d/dz g = E2``. Let ``gamma=self`` and let
        ``M_gamma(z)`` be ``Log((c*z+d) * sgn(a+d))`` if ``c, a+d > 0``,
        resp. ``Log((c*z+d) / i*sgn(c))`` if ``a+d = 0, c!=0``,
        resp. ``0`` if ``c=0``. Let ``k=4 * n / (n-2)``, then:

        ``g(gamma.acton(z) - g(z) - k*M_gamma(z)`` is equal to
        ``2*pi*i / (n-2) * self.linking_number()``.

        In particular it is independent of ``z`` and a conjugacy invariant.

        If ``self`` is hyperbolic then in the classical case ``n=3``
        this is the linking number of the closed geodesic
        (corresponding to ``self``) with the trefoil knot.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms

            sage: def E2_primitive(z, n=3, prec=10, num_prec=53):
            ....:     G = HeckeTriangleGroup(n=n)
            ....:     MF = QuasiModularForms(group=G, k=2, ep=-1)
            ....:     q = MF.get_q(prec=prec)
            ....:     int_series = integrate((MF.E2().q_expansion(prec=prec) - 1) / q)
            ....:     t_const = (2*pi*i/G.lam()).n(num_prec)
            ....:     d = MF.get_d(fix_d=True, d_num_prec=num_prec)
            ....:     q = exp(t_const * z)
            ....:     return t_const*z + sum([(int_series.coefficients()[m]).subs(d=d) * q**int_series.exponents()[m] for m in range(len(int_series.coefficients()))])

            sage: def M(gamma, z, num_prec=53):
            ....:     a = ComplexField(num_prec)(gamma.a())
            ....:     b = ComplexField(num_prec)(gamma.b())
            ....:     c = ComplexField(num_prec)(gamma.c())
            ....:     d = ComplexField(num_prec)(gamma.d())
            ....:     if c == 0:
            ....:         return 0
            ....:     elif a + d == 0:
            ....:         return log(-i.n(num_prec)*(c*z + d)*sign(c))
            ....:     else:
            ....:         return log((c*z+d)*sign(a+d))

            sage: def num_linking_number(A, z, n=3, prec=10, num_prec=53):
            ....:     z = z.n(num_prec)
            ....:     k = 4 * n / (n - 2)
            ....:     return (n-2) / (2*pi*i).n(num_prec) * (E2_primitive(A.acton(z), n=n, prec=prec, num_prec=num_prec) - E2_primitive(z, n=n, prec=prec, num_prec=num_prec) - k*M(A, z, num_prec=num_prec))

            sage: G = HeckeTriangleGroup(8)
            sage: z = i
            sage: for A in [G.S(), G.T(), G.U(), G.U()^(G.n()//2), G.U()^(-3)]:
            ....:     print("A={}: ".format(A.string_repr("conj")))
            ....:     num_linking_number(A, z, G.n())
            ....:     A.linking_number()
            A=[S]:
            0.000000000000...
            0
            A=[V(1)]:
            6.000000000000...
            6
            A=[U]:
            -2.00000000000...
            -2
            A=[U^4]:
            0.596987639289... + 0.926018962976...*I
            0
            A=[U^(-3)]:
            5.40301236071... + 0.926018962976...*I
            6

            sage: z = ComplexField(1000)(- 2.3 + 3.1*i)
            sage: B = G.I()
            sage: for A in [G.S(), G.T(), G.U(), G.U()^(G.n()//2), G.U()^(-3)]:
            ....:     print("A={}: ".format(A.string_repr("conj")))
            ....:     num_linking_number(B.acton(A), z, G.n(), prec=100, num_prec=1000).n(53)
            ....:     B.acton(A).linking_number()
            A=[S]:
            6.63923483989...e-31 + 2.45195568651...e-30*I
            0
            A=[V(1)]:
            6.000000000000...
            6
            A=[U]:
            -2.00000000000... + 2.45195568651...e-30*I
            -2
            A=[U^4]:
            0.00772492873864... + 0.00668936643212...*I
            0
            A=[U^(-3)]:
            5.99730551444... + 0.000847636355069...*I
            6

            sage: z = ComplexField(5000)(- 2.3 + 3.1*i)
            sage: B = G.U()
            sage: for A in [G.S(), G.T(), G.U(), G.U()^(G.n()//2), G.U()^(-3)]:    # long time
            ....:     print("A={}: ".format(A.string_repr("conj")))
            ....:     num_linking_number(B.acton(A), z, G.n(), prec=200, num_prec=5000).n(53)
            ....:     B.acton(A).linking_number()
            A=[S]:
            -7.90944791339...e-34 - 9.38956758807...e-34*I
            0
            A=[V(1)]:
            5.99999997397... - 5.96520311160...e-8*I
            6
            A=[U]:
            -2.00000000000... - 1.33113963568...e-61*I
            -2
            A=[U^4]:
            -2.32704571946...e-6 + 5.91899385948...e-7*I
            0
            A=[U^(-3)]:
            6.00000032148... - 1.82676936467...e-7*I
            6

            sage: A = G.V(2)*G.V(3)
            sage: B = G.I()
            sage: num_linking_number(B.acton(A), z, G.n(), prec=200, num_prec=5000).n(53)    # long time
            6.00498424588... - 0.00702329345176...*I
            sage: A.linking_number()
            6

            The numerical properties for anything larger are basically
            too bad to make nice further tests...
        """
        if self.is_identity():
            return ZZ.zero()

        (L, R, sgn) = self._block_decomposition_data()
        n = self.parent().n()

        if self.is_elliptic():
            if L[0] == 0:
                return ZZ(0)
            elif 2*L[1] == n:
                return ZZ(0)
            else:
                return ZZ(-2*L[1])
        else:
            t = sum(v[1] for v in L)
            u = sum((v[0]-1) for v in L)

            return ZZ((n-2)*t - 2*u)

    def root_extension_field(self):
        r"""
        Return a field extension which contains the fixed points of ``self``.
        Namely the root extension field of the parent for the discriminant of ``self``.
        Also see the parent method ``root_extension_field(D)`` and
        :meth:`root_extension_embedding` (which provides the correct embedding).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)
            sage: G.V(3).discriminant()
            32
            sage: G.V(3).root_extension_field() == G.root_extension_field(32)
            True
            sage: G.T().root_extension_field() == G.root_extension_field(G.T().discriminant()) == G.base_field()
            True
            sage: (G.S()).root_extension_field() == G.root_extension_field(G.S().discriminant())
            True

            sage: G = HeckeTriangleGroup(n=7)
            sage: D = G.V(3).discriminant()
            sage: D
            4*lam^2 + 4*lam - 4
            sage: G.V(3).root_extension_field() == G.root_extension_field(D)
            True
            sage: G.U().root_extension_field() == G.root_extension_field(G.U().discriminant())
            True
            sage: G.V(1).root_extension_field() == G.base_field()
            True
        """
        return self.parent().root_extension_field(self.discriminant())

    def root_extension_embedding(self, K=None):
        r"""
        Return the correct embedding from the root extension field to ``K``.

        INPUT:

        - ``K`` -- A field to which we want the (correct) embedding.
                   If ``K=None`` (default) then ``AlgebraicField()`` is
                   used for elliptic elements and ``AlgebraicRealField()``
                   otherwise.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)

            sage: fp = (-G.S()).fixed_points()[0]
            sage: alg_fp = (-G.S()).root_extension_embedding()(fp)
            sage: alg_fp
            1*I
            sage: alg_fp == (-G.S()).fixed_points(embedded=True)[0]
            True

            sage: fp = (-G.V(2)).fixed_points()[1]
            sage: alg_fp = (-G.V(2)).root_extension_embedding()(fp)
            sage: alg_fp
            -1.732050807568...?
            sage: alg_fp == (-G.V(2)).fixed_points(embedded=True)[1]
            True

            sage: fp = (-G.V(2)).fixed_points()[0]
            sage: alg_fp = (-G.V(2)).root_extension_embedding()(fp)
            sage: alg_fp
            1.732050807568...?
            sage: alg_fp == (-G.V(2)).fixed_points(embedded=True)[0]
            True

            sage: G = HeckeTriangleGroup(n=7)

            sage: fp = (-G.S()).fixed_points()[1]
            sage: alg_fp = (-G.S()).root_extension_embedding()(fp)
            sage: alg_fp
            0.?... - 1.000000000000...?*I
            sage: alg_fp == (-G.S()).fixed_points(embedded=True)[1]
            True

            sage: fp = (-G.U()^4).fixed_points()[0]
            sage: alg_fp = (-G.U()^4).root_extension_embedding()(fp)
            sage: alg_fp
            0.9009688679024...? + 0.4338837391175...?*I
            sage: alg_fp == (-G.U()^4).fixed_points(embedded=True)[0]
            True

            sage: (-G.U()^4).root_extension_embedding(CC)(fp)
            0.900968867902... + 0.433883739117...*I
            sage: (-G.U()^4).root_extension_embedding(CC)(fp).parent()
            Complex Field with 53 bits of precision

            sage: fp = (-G.V(5)).fixed_points()[1]
            sage: alg_fp = (-G.V(5)).root_extension_embedding()(fp)
            sage: alg_fp
            -0.6671145837954...?
            sage: alg_fp == (-G.V(5)).fixed_points(embedded=True)[1]
            True
        """
        return self.parent().root_extension_embedding(self.discriminant(), K)

    def fixed_points(self, embedded=False, order="default"):
        r"""
        Return a pair of (mutually conjugate) fixed points of ``self``
        in a possible quadratic extension of the base field.

        INPUT:

        - ``embedded`` -- If ``True`` the fixed points are embedded into
                          ``AlgebraicRealField`` resp. ``AlgebraicField``.
                          Default: ``False``.

        - ``order``    -- If ``order="none"`` the fixed points are choosen
                          and ordered according to a fixed formula.

                          If ``order="sign"`` the fixed points are always ordered
                          according to the sign in front of the square root.

                          If ``order="default"`` (default) then in case the fixed
                          points are hyperbolic they are ordered according to the
                          sign of the trace of ``self`` instead, such that the
                          attracting fixed point comes first.

                          If ``order="trace"`` the fixed points are always ordered
                          according to the sign of the trace of ``self``.
                          If the trace is zero they are ordered by the sign in
                          front of the square root. In particular the fixed_points
                          in this case remain the same for ``-self``.

        OUTPUT:

        If ``embedded=True`` an element of either ``AlgebraicRealField`` or
        ``AlgebraicField`` is returned. Otherwise an element of a relative field
        extension over the base field of (the parent of) ``self`` is returned.

        Warning: Relative field extensions don't support default embeddings.
        So the correct embedding (which is the positive resp. imaginary positive
        one) has to be choosen.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)
            sage: (-G.T(-4)).fixed_points()
            (+Infinity, +Infinity)
            sage: (-G.S()).fixed_points()
            (1/2*e, -1/2*e)
            sage: p = (-G.S()).fixed_points(embedded=True)[0]
            sage: p
            I
            sage: (-G.S()).acton(p) == p
            True
            sage: (-G.V(2)).fixed_points()
            (1/2*e, -1/2*e)
            sage: (-G.V(2)).fixed_points() == G.V(2).fixed_points()
            True
            sage: p = (-G.V(2)).fixed_points(embedded=True)[1]
            sage: p
            -1.732050807568878?
            sage: (-G.V(2)).acton(p) == p
            True

            sage: G = HeckeTriangleGroup(n=7)
            sage: (-G.S()).fixed_points()
            (1/2*e, -1/2*e)
            sage: p = (-G.S()).fixed_points(embedded=True)[1]
            sage: p
            -I
            sage: (-G.S()).acton(p) == p
            True
            sage: (G.U()^4).fixed_points()
            ((1/2*lam^2 - 1/2*lam - 1/2)*e + 1/2*lam, (-1/2*lam^2 + 1/2*lam + 1/2)*e + 1/2*lam)
            sage: pts = (G.U()^4).fixed_points(order="trace")
            sage: (G.U()^4).fixed_points() == [pts[1], pts[0]]
            False
            sage: (G.U()^4).fixed_points(order="trace") == (-G.U()^4).fixed_points(order="trace")
            True
            sage: (G.U()^4).fixed_points() == (G.U()^4).fixed_points(order="none")
            True
            sage: (-G.U()^4).fixed_points() == (G.U()^4).fixed_points()
            True
            sage: (-G.U()^4).fixed_points(order="none") == pts
            True
            sage: p = (G.U()^4).fixed_points(embedded=True)[1]
            sage: p
            0.9009688679024191? - 0.4338837391175581?*I
            sage: (G.U()^4).acton(p) == p
            True
            sage: (-G.V(5)).fixed_points()
            ((1/2*lam^2 - 1/2*lam - 1/2)*e, (-1/2*lam^2 + 1/2*lam + 1/2)*e)
            sage: (-G.V(5)).fixed_points() == G.V(5).fixed_points()
            True
            sage: p = (-G.V(5)).fixed_points(embedded=True)[0]
            sage: p
            0.6671145837954892?
            sage: (-G.V(5)).acton(p) == p
            True
        """
        if self.c() == 0:
            return (infinity, infinity)
        else:
            D = self.discriminant()
            if D.is_square():
                e = D.sqrt()
            else:
                e = self.root_extension_field().gen()

            a,b,c,d = self._matrix.list()

            if order == "none":
                sgn = ZZ(1)
            elif order == "sign":
                sgn = coerce_AA(c).sign()
            elif order == "default":
                if self.is_elliptic() or self.trace() == 0:
                    sgn = coerce_AA(c).sign()
                else:
                    sgn = coerce_AA(self.trace()).sign()
            elif order == "trace":
                if self.trace() == 0:
                    sgn = coerce_AA(c).sign()
                else:
                    sgn = coerce_AA(self.trace()).sign()
            else:
                raise NotImplementedError

            if embedded:
                e = coerce_AA(D).sqrt()
                e.simplify()
                a = coerce_AA(a)
                d = coerce_AA(d)
                c = coerce_AA(c)

            root1 = (a-d)/(2*c) + sgn*e/(2*c)
            root2 = (a-d)/(2*c) - sgn*e/(2*c)

            if embedded:
                root1.simplify()
                root1.exactify()
                root2.simplify()
                root2.exactify()

            return (root1, root2)

    def acton(self, tau):
        r"""
        Return the image of ``tau`` under the action of ``self``
        by linear fractional transformations or by conjugation
        in case ``tau`` is an element of the parent of ``self``.

        It is possible to act on points of ``HyperbolicPlane()``.

        .. NOTE::

            There is a 1-1 correspondence between hyperbolic
            fixed points and the corresponding primitive element
            in the stabilizer. The action in the two cases above
            is compatible with this correspondence.

        INPUT:

        - ``tau``   -- Either an element of ``self`` or any
                       element to which a linear fractional
                       transformation can be applied in
                       the usual way.

                       In particular ``infinity`` is a possible
                       argument and a possible return value.

                       As mentioned it is also possible to use
                       points of ``HyperbolicPlane()``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(5)
            sage: G.S().acton(SR(1 + i/2))
            2/5*I - 4/5
            sage: G.S().acton(SR(1 + i/2)).parent()
            Symbolic Ring
            sage: G.S().acton(QQbar(1 + i/2))
            2/5*I - 4/5
            sage: G.S().acton(QQbar(1 + i/2)).parent()
            Algebraic Field

            sage: G.S().acton(i + exp(-2))
            -1/(e^(-2) + I)
            sage: G.S().acton(i + exp(-2)).parent()
            Symbolic Ring

            sage: G.T().acton(infinity) == infinity
            True
            sage: G.U().acton(infinity)
            lam
            sage: G.V(2).acton(-G.lam()) == infinity
            True

            sage: G.V(2).acton(G.U()) == G.V(2)*G.U()*G.V(2).inverse()
            True
            sage: G.V(2).inverse().acton(G.U())
            [  0  -1]
            [  1 lam]

            sage: p = HyperbolicPlane().PD().get_point(-I/2+1/8)
            sage: G.V(2).acton(p)
            Point in PD -((-(47*I + 161)*sqrt(5) - 47*I - 161)/(145*sqrt(5) + 94*I + 177) + I)/(I*(-(47*I + 161)*sqrt(5) - 47*I - 161)/(145*sqrt(5) + 94*I + 177) + 1)
            sage: bool(G.V(2).acton(p).to_model('UHP').coordinates() == G.V(2).acton(p.to_model('UHP').coordinates()))
            True

            sage: p = HyperbolicPlane().PD().get_point(I)
            sage: G.U().acton(p)
            Boundary point in PD 1/2*(sqrt(5) - 2*I + 1)/(-1/2*I*sqrt(5) - 1/2*I + 1)
            sage: G.U().acton(p).to_model('UHP') == HyperbolicPlane().UHP().get_point(G.lam())
            True
            sage: G.U().acton(p) == HyperbolicPlane().UHP().get_point(G.lam()).to_model('PD')
            True
        """

        if tau.parent() == self.parent():
            return self*tau*self.inverse()

        # if tau is a point of HyperbolicPlane then we use it's coordinates in the UHP model
        model = None
        if (tau in HyperbolicPlane()):
            model = tau.model()
            tau = tau.to_model('UHP').coordinates()

        a,b,c,d = self._matrix.list()

        if tau == infinity:
            if c.is_zero():
                result = infinity
            else:
                result = a/c
        elif c*tau + d == 0:
            result = infinity
        else:
            result = (a*tau + b) / (c*tau + d)

        if model is None:
            return result
        else:
            return HyperbolicPlane().UHP().get_point(result).to_model(model)

    def _act_on_(self, other, self_on_left):
        r"""
        Define the action by linear fractional transformation of Hecke triangle group
        elements on complex points (using :meth:`acton`).

        For the action on matrices by conjugation :meth:`acton` has to be used explicitly
        (to avoid confusion/ambiguity in expressions of the form gamma1*gamma2*z).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(5)
            sage: p = HyperbolicPlane().PD().get_point(I)
            sage: G.U()*p
            Boundary point in PD 1/2*(sqrt(5) - 2*I + 1)/(-1/2*I*sqrt(5) - 1/2*I + 1)
            sage: G.S()*G.U()*p == G.S()*(G.U()*p)
            True
            sage: G.S()*G.U()*p == (G.S()*G.U())*p
            True
            sage: (G.S()*G.U())*p == G.S()*(G.U()*p)
            True

            sage: p = G.lam()
            sage: G.U()*G.T()*p
            1/2*lam + 1/2
            sage: p = QQbar(i*sqrt(2))
            sage: G.U()*p
            1.618033988749895? + 0.7071067811865475?*I
            sage: p = CC(-i + sqrt(2))
            sage: G.U()*p
            1.14662946795886 - 0.333333333333333*I
            sage: p = infinity
            sage: G.U()*p
            lam
        """

        if (self_on_left):
            if (other == infinity or other in CC or other in HyperbolicPlane()):
                return self.acton(other)
        return None

    def slash(self, f, tau=None, k=None):
        r"""
        Return the `slash-operator` of weight ``k`` to applied to ``f``,
        evaluated at ``tau``. I.e. ``(f|_k[self])(tau)``.

        INPUT:

        - ``f``   -- A function in ``tau`` (or an object for which
          evaluation at ``self.acton(tau)`` makes sense.

        - ``tau`` -- Where to evaluate the result.
          This should be a valid argument for :meth:`acton`.

          If ``tau`` is a point of ``HyperbolicPlane()`` then
          its coordinates in the upper half plane model are used.

          Default: ``None`` in which case ``f`` has to be
          a rational function / polynomial in one variable and
          the generator of the polynomial ring is used for ``tau``.
          That way ``slash`` acts on rational functions / polynomials.

        - ``k``   -- An even integer.

          Default: ``None`` in which case ``f`` either
          has to be a rational function / polynomial in one
          variable (then -degree is used).
          Or ``f`` needs to have a ``weight`` attribute which
          is then used.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: G = HeckeTriangleGroup(n=5)
            sage: E4 = ModularForms(group=G, k=4, ep=1).E4()
            sage: z = CC(-1/(-1/(2*i+30)-1))
            sage: (G.S()).slash(E4, z)
            32288.0558881... - 118329.856601...*I
            sage: (G.V(2)*G.V(3)).slash(E4, z)
            32288.0558892... - 118329.856603...*I
            sage: E4(z)
            32288.0558881... - 118329.856601...*I

            sage: z = HyperbolicPlane().PD().get_point(CC(-I/2 + 1/8))
            sage: (G.V(2)*G.V(3)).slash(E4, z)
            -(21624.437... - 12725.035...*I)/((0.610... + 0.324...*I)*sqrt(5) + 2.720... + 0.648...*I)^4

            sage: z = PolynomialRing(G.base_ring(), 'z').gen()
            sage: rat = z^2 + 1/(z-G.lam())
            sage: dr = rat.numerator().degree() - rat.denominator().degree()
            sage: G.S().slash(rat) == G.S().slash(rat, tau=None, k=-dr)
            True
            sage: G.S().slash(rat)
            (z^6 - lam*z^4 - z^3)/(-lam*z^4 - z^3)
            sage: G.S().slash(rat, k=0)
            (z^4 - lam*z^2 - z)/(-lam*z^4 - z^3)
            sage: G.S().slash(rat, k=-4)
            (z^8 - lam*z^6 - z^5)/(-lam*z^4 - z^3)
        """

        if k is None:
            if hasattr(f, 'weight'):
                k = f.weight()
            else:
                try:
                    par = f.numerator().parent()
                    degf = par(f.numerator()).degree() - par(f.denominator()).degree()
                except (ValueError, TypeError, AttributeError):
                    raise ValueError("The weight k could not be determined automatically and needs to be specified manually!")
                k = -degf

        try:
            k = ZZ(k)
            if k % 2:
                raise TypeError
        except TypeError:
            raise ValueError("k={} must be an even integer!".format(k))

        if tau is None:
            try:
                tau = f.numerator().parent().gen()
            except (ValueError, TypeError, AttributeError):
                raise ValueError("f={} is not a rational function or a polynomial in one variable, so tau has to be specified explicitly!".format(f))

        if (tau in HyperbolicPlane()):
            tau = tau.to_model('UHP').coordinates()

        return (self.c()*tau + self.d())**(-k) * f(self.acton(tau))

    def as_hyperbolic_plane_isometry(self, model="UHP"):
        r"""
        Return ``self`` as an isometry of ``HyperbolicPlane()`` (in the upper half plane model).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: el = HeckeTriangleGroup(7).V(4)
            sage: el.as_hyperbolic_plane_isometry()
            Isometry in UHP
            [lam^2 - 1       lam]
            [lam^2 - 1 lam^2 - 1]
            sage: el.as_hyperbolic_plane_isometry().parent()
            Set of Morphisms from Hyperbolic plane in the Upper Half Plane Model to Hyperbolic plane in the Upper Half Plane Model in Category of hyperbolic models of Hyperbolic plane
            sage: el.as_hyperbolic_plane_isometry("KM").parent()
            Set of Morphisms from Hyperbolic plane in the Klein Disk Model to Hyperbolic plane in the Klein Disk Model in Category of hyperbolic models of Hyperbolic plane
        """
        from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane
        return HyperbolicPlane().UHP().get_isometry(self._matrix).to_model(model)
