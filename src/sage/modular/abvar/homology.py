# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint sage.libs.pari
r"""
Homology of modular abelian varieties

Sage can compute with homology groups associated to modular abelian
varieties with coefficients in any commutative ring. Supported
operations include computing matrices and characteristic
polynomials of Hecke operators, rank, and rational decomposition as
a direct sum of factors (obtained by cutting out kernels of Hecke
operators).

AUTHORS:

- William Stein (2007-03)

EXAMPLES::

    sage: J = J0(43)
    sage: H = J.integral_homology()
    sage: H
    Integral Homology of Abelian variety J0(43) of dimension 3
    sage: H.hecke_matrix(19)
    [ 0  0 -2  0  2  0]
    [ 2 -4 -2  0  2  0]
    [ 0  0 -2 -2  0  0]
    [ 2  0 -2 -4  2 -2]
    [ 0  2  0 -2 -2  0]
    [ 0  2  0 -2  0  0]
    sage: H.base_ring()
    Integer Ring
    sage: d = H.decomposition(); d
    [
    Submodule of rank 2 of Integral Homology of Abelian variety J0(43) of dimension 3,
    Submodule of rank 4 of Integral Homology of Abelian variety J0(43) of dimension 3
    ]
    sage: a = d[0]
    sage: a.hecke_matrix(5)
    [-4  0]
    [ 0 -4]
    sage: a.T(7)
    Hecke operator T_7 on
     Submodule of rank 2 of Integral Homology of Abelian variety J0(43) of dimension 3
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.commutative_rings import CommutativeRings
from sage.modular.hecke.module import HeckeModule_free_module
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.richcmp import richcmp_method, richcmp, richcmp_not_equal

# TODO: we will probably also need homology that is *not* a Hecke module.


@richcmp_method
class Homology(HeckeModule_free_module):
    """
    A homology group of an abelian variety, equipped with a Hecke
    action.
    """
    def hecke_polynomial(self, n, var='x'):
        r"""
        Return the `n`-th Hecke polynomial in the given variable.

        INPUT:

        - ``n`` -- positive integer

        - ``var`` -- string (default: ``'x'``); the variable name

        OUTPUT: a polynomial over `\ZZ` in the given variable

        EXAMPLES::

            sage: H = J0(43).integral_homology(); H
            Integral Homology of Abelian variety J0(43) of dimension 3
            sage: f = H.hecke_polynomial(3); f
            x^6 + 4*x^5 - 16*x^3 - 12*x^2 + 16*x + 16
            sage: parent(f)
            Univariate Polynomial Ring in x over Integer Ring
            sage: H.hecke_polynomial(3,'w')
            w^6 + 4*w^5 - 16*w^3 - 12*w^2 + 16*w + 16
        """
        return self.hecke_matrix(n).charpoly(var)


class Homology_abvar(Homology):
    """
    The homology of a modular abelian variety.
    """
    def __init__(self, abvar, base):
        """
        This is an abstract base class, so it is called implicitly in the
        following examples.

        EXAMPLES::

            sage: H = J0(43).integral_homology()
            sage: type(H)
            <class 'sage.modular.abvar.homology.IntegralHomology_with_category'>

        TESTS::

            sage: H = J0(43).integral_homology()
            sage: loads(dumps(H)) == H
            True
        """
        if base not in CommutativeRings():
            raise TypeError("base ring must be a commutative ring")
        HeckeModule_free_module.__init__(
            self, base, abvar.level(), weight=2)
        self.__abvar = abvar

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: J0(37).integral_homology() == J0(41).integral_homology()
            False
            sage: J0(37).integral_homology() == J0(37).rational_homology()
            False
            sage: J0(37).integral_homology() == loads(dumps(J0(37).integral_homology()))
            True
        """
        if not isinstance(other, Homology_abvar):
            return NotImplemented
        else:
            return richcmp((self.abelian_variety(), self.base_ring()),
                           (other.abelian_variety(), other.base_ring()), op)

    def _repr_(self):
        """
        Return string representation of ``self``. This must be defined in the
        derived class.

        EXAMPLES::

            sage: H = J0(43).integral_homology()
            sage: from sage.modular.abvar.homology import Homology_abvar
            sage: Homology_abvar._repr_(H)
            Traceback (most recent call last):
            ...
            NotImplementedError: please override this in the derived class
        """
        raise NotImplementedError("please override this in the derived class")

    def gens(self):
        """
        Return generators of ``self``.

        This is not yet implemented!

        EXAMPLES::

            sage: H = J0(37).homology()
            sage: H.gens()    # this will change
            Traceback (most recent call last):
            ...
            NotImplementedError: homology classes not yet implemented
        """
        raise NotImplementedError("homology classes not yet implemented")

    def gen(self, n):
        """
        Return `n`-th generator of ``self``.

        This is not yet implemented!

        EXAMPLES::

            sage: H = J0(37).homology()
            sage: H.gen(0)    # this will change
            Traceback (most recent call last):
            ...
            NotImplementedError: homology classes not yet implemented
        """
        raise NotImplementedError("homology classes not yet implemented")

    def abelian_variety(self):
        """
        Return the abelian variety that this is the homology of.

        EXAMPLES::

            sage: H = J0(48).homology()
            sage: H.abelian_variety()
            Abelian variety J0(48) of dimension 3
        """
        return self.__abvar

    def ambient_hecke_module(self):
        """
        Return the ambient Hecke module that this homology is contained
        in.

        EXAMPLES::

            sage: H = J0(48).homology(); H
            Integral Homology of Abelian variety J0(48) of dimension 3
            sage: H.ambient_hecke_module()
            Integral Homology of Abelian variety J0(48) of dimension 3
        """
        return self

    def free_module(self):
        """
        Return the underlying free module of this homology group.

        EXAMPLES::

            sage: H = J0(48).homology()
            sage: H.free_module()
            Ambient free module of rank 6 over the principal ideal domain Integer Ring
        """
        try:
            return self.__free_module
        except AttributeError:
            M = self.base_ring()**self.rank()
            self.__free_module = M
            return M

    def hecke_bound(self):
        r"""
        Return bound on the number of Hecke operators needed to generate
        the Hecke algebra as a `\ZZ`-module acting on this
        space.

        EXAMPLES::

            sage: J0(48).homology().hecke_bound()
            16
            sage: J1(15).homology().hecke_bound()
            32
        """
        return self.__abvar.modular_symbols(sign=1).hecke_bound()

    def hecke_matrix(self, n):
        """
        Return the matrix of the `n`-th Hecke operator acting on this
        homology group.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: a matrix over the coefficient ring of this homology group

        EXAMPLES::

            sage: H = J0(23).integral_homology()
            sage: H.hecke_matrix(3)
            [-1 -2  2  0]
            [ 0 -3  2 -2]
            [ 2 -4  3 -2]
            [ 2 -2  0  1]

        The matrix is over the coefficient ring::

            sage: J = J0(23)
            sage: J.homology(QQ[I]).hecke_matrix(3).parent()
            Full MatrixSpace of 4 by 4 dense matrices over
             Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        raise NotImplementedError

    def rank(self):
        """
        Return the rank as a module or vector space of this homology
        group.

        EXAMPLES::

            sage: H = J0(5077).homology(); H
            Integral Homology of Abelian variety J0(5077) of dimension 422
            sage: H.rank()
            844
        """
        return self.__abvar.dimension() * 2

    def submodule(self, U, check=True):
        r"""
        Return the submodule of this homology group given by `U`,
        which should be a submodule of the free module associated to this
        homology group.

        INPUT:

        - ``U`` -- submodule of ambient free module (or
          something that defines one)

        - ``check`` -- currently ignored


        .. NOTE::

           We do *not* check that U is invariant under all Hecke
           operators.

        EXAMPLES::

            sage: H = J0(23).homology(); H
            Integral Homology of Abelian variety J0(23) of dimension 2
            sage: F = H.free_module()
            sage: U = F.span([[1,2,3,4]])
            sage: M = H.submodule(U); M
            Submodule of rank 1 of Integral Homology of Abelian variety J0(23) of dimension 2

        Note that the submodule command doesn't actually check that the
        object defined is a homology group or is invariant under the Hecke
        operators. For example, the fairly random `M` that we just
        defined is not invariant under the Hecke operators, so it is not a
        Hecke submodule - it is only a `\ZZ`-submodule.

        ::

            sage: M.hecke_matrix(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: subspace is not invariant under matrix
        """
        return Homology_submodule(self, U)


class IntegralHomology(Homology_abvar):
    r"""
    The integral homology `H_1(A,\ZZ)` of a modular
    abelian variety.
    """
    def __init__(self, abvar):
        """
        Create the integral homology of a modular abelian variety.

        INPUT:

        - ``abvar`` -- a modular abelian variety

        EXAMPLES::

            sage: H = J0(23).integral_homology(); H
            Integral Homology of Abelian variety J0(23) of dimension 2
            sage: type(H)
            <class 'sage.modular.abvar.homology.IntegralHomology_with_category'>

        TESTS::

            sage: loads(dumps(H)) == H
            True
        """
        Homology_abvar.__init__(self, abvar, ZZ)

    def _repr_(self):
        """
        String representation of the integral homology.

        EXAMPLES::

            sage: J0(23).integral_homology()._repr_()
            'Integral Homology of Abelian variety J0(23) of dimension 2'
        """
        return "Integral Homology of %s" % self.abelian_variety()

    def hecke_matrix(self, n):
        """
        Return the matrix of the `n`-th Hecke operator acting on this
        homology group.

        EXAMPLES::

            sage: J0(48).integral_homology().hecke_bound()
            16
            sage: t = J1(13).integral_homology().hecke_matrix(3); t
            [-2  2  2 -2]
            [-2  0  2  0]
            [ 0  0  0 -2]
            [ 0  0  2 -2]
            sage: t.base_ring()
            Integer Ring
        """
        n = Integer(n)
        return self.abelian_variety()._integral_hecke_matrix(n)

    def hecke_polynomial(self, n, var='x'):
        """
        Return the `n`-th Hecke polynomial on this integral homology group.

        EXAMPLES::

            sage: f = J0(43).integral_homology().hecke_polynomial(2)
            sage: f.base_ring()
            Integer Ring
            sage: factor(f)
            (x + 2)^2 * (x^2 - 2)^2
        """
        n = Integer(n)
        M = self.abelian_variety().modular_symbols(sign=1)
        f = (M.hecke_polynomial(n, var)**2).change_ring(ZZ)
        return f


class RationalHomology(Homology_abvar):
    r"""
    The rational homology `H_1(A,\QQ)` of a modular
    abelian variety.
    """
    def __init__(self, abvar):
        """
        Create the rational homology of a modular abelian variety.

        INPUT:

        - ``abvar`` -- a modular abelian variety

        EXAMPLES::

            sage: H = J0(23).rational_homology(); H
            Rational Homology of Abelian variety J0(23) of dimension 2

        TESTS::

            sage: loads(dumps(H)) == H
            True
        """
        Homology_abvar.__init__(self, abvar, QQ)

    def _repr_(self):
        """
        Return string representation of the rational homology.

        EXAMPLES::

            sage: J0(23).rational_homology()._repr_()
            'Rational Homology of Abelian variety J0(23) of dimension 2'
        """
        return "Rational Homology of %s" % self.abelian_variety()

    def hecke_matrix(self, n):
        """
        Return the matrix of the `n`-th Hecke operator acting on this
        homology group.

        EXAMPLES::

            sage: t = J1(13).homology(QQ).hecke_matrix(3); t
            [-2  2  2 -2]
            [-2  0  2  0]
            [ 0  0  0 -2]
            [ 0  0  2 -2]
            sage: t.base_ring()
            Rational Field
            sage: t = J1(13).homology(GF(3)).hecke_matrix(3); t
            [1 2 2 1]
            [1 0 2 0]
            [0 0 0 1]
            [0 0 2 1]
            sage: t.base_ring()
            Finite Field of size 3
        """
        n = Integer(n)
        return self.abelian_variety()._rational_hecke_matrix(n)

    def hecke_polynomial(self, n, var='x'):
        """
        Return the `n`-th Hecke polynomial on this rational homology group.

        EXAMPLES::

            sage: f = J0(43).rational_homology().hecke_polynomial(2)
            sage: f.base_ring()
            Rational Field
            sage: factor(f)
            (x + 2) * (x^2 - 2)
        """
        f = self.hecke_operator(n).matrix().characteristic_polynomial(var)
        _, poly = f.is_square(True)
        return poly


class Homology_over_base(Homology_abvar):
    r"""
    The homology over a modular abelian variety over an arbitrary base
    commutative ring (not `\ZZ` or `\QQ`).
    """
    def __init__(self, abvar, base_ring):
        r"""
        Called when creating homology with coefficients not
        `\ZZ` or `\QQ`.

        INPUT:

        - ``abvar`` -- a modular abelian variety

        - ``base_ring`` -- a commutative ring

        EXAMPLES::

            sage: H = J0(23).homology(GF(5)); H
            Homology with coefficients in Finite Field of size 5 of Abelian variety J0(23) of dimension 2
            sage: type(H)
            <class 'sage.modular.abvar.homology.Homology_over_base_with_category'>

        TESTS::

            sage: loads(dumps(H)) == H
            True
        """
        Homology_abvar.__init__(self, abvar, base_ring)

    def _repr_(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: H = J0(23).homology(GF(5))
            sage: H._repr_()
            'Homology with coefficients in Finite Field of size 5 of Abelian variety J0(23) of dimension 2'
        """
        return "Homology with coefficients in %s of %s" % (self.base_ring(), self.abelian_variety())

    def hecke_matrix(self, n):
        """
        Return the matrix of the `n`-th Hecke operator acting on this
        homology group.

        EXAMPLES::

            sage: t = J1(13).homology(GF(3)).hecke_matrix(3); t
            [1 2 2 1]
            [1 0 2 0]
            [0 0 0 1]
            [0 0 2 1]
            sage: t.base_ring()
            Finite Field of size 3
        """
        n = Integer(n)
        return self.abelian_variety()._integral_hecke_matrix(n).change_ring(self.base_ring())


class Homology_submodule(Homology):
    """
    A submodule of the homology of a modular abelian variety.
    """
    def __init__(self, ambient, submodule):
        """
        Create a submodule of the homology of a modular abelian variety.

        INPUT:

        - ``ambient`` -- the homology of some modular abelian
          variety with ring coefficients

        - ``submodule`` -- a submodule of the free module
          underlying ambient

        EXAMPLES::

            sage: H = J0(37).homology()
            sage: H.submodule([[1,0,0,0]])
            Submodule of rank 1 of Integral Homology of Abelian variety J0(37) of dimension 2

        TESTS::

            sage: loads(dumps(H)) == H
            True
        """
        if not isinstance(ambient, Homology_abvar):
            raise TypeError("ambient must be the homology of a modular abelian variety")
        self.__ambient = ambient
        submodule = ambient.free_module().submodule(submodule)
        self.__submodule = submodule
        HeckeModule_free_module.__init__(
            self, ambient.base_ring(), ambient.level(), weight=2)

    def _repr_(self):
        """
        String representation of this submodule of homology.

        EXAMPLES::

            sage: H = J0(37).homology()
            sage: G = H.submodule([[1, 2, 3, 4]])
            sage: G._repr_()
            'Submodule of rank 1 of Integral Homology of Abelian variety J0(37) of dimension 2'
        """
        return "Submodule of rank %s of %s" % (self.rank(), self.__ambient)

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: J0(37).homology().decomposition() # indirect doctest
            [
            Submodule of rank 2 of Integral Homology of Abelian variety J0(37) of dimension 2,
            Submodule of rank 2 of Integral Homology of Abelian variety J0(37) of dimension 2
            ]
        """
        if not isinstance(other, Homology_submodule):
            return NotImplemented
        lx = self.__ambient
        rx = other.__ambient
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)
        return richcmp(self.__submodule, other.__submodule, op)

    def ambient_hecke_module(self):
        """
        Return the ambient Hecke module that this homology is contained
        in.

        EXAMPLES::

            sage: H = J0(48).homology(); H
            Integral Homology of Abelian variety J0(48) of dimension 3
            sage: d = H.decomposition(); d
            [
            Submodule of rank 2 of Integral Homology of Abelian variety J0(48) of dimension 3,
            Submodule of rank 4 of Integral Homology of Abelian variety J0(48) of dimension 3
            ]
            sage: d[0].ambient_hecke_module()
            Integral Homology of Abelian variety J0(48) of dimension 3
        """
        return self.__ambient

    def free_module(self):
        """
        Return the underlying free module of the homology group.

        EXAMPLES::

            sage: H = J0(48).homology()
            sage: K = H.decomposition()[1]; K
            Submodule of rank 4 of Integral Homology of Abelian variety J0(48) of dimension 3
            sage: K.free_module()
            Free module of degree 6 and rank 4 over Integer Ring
            Echelon basis matrix:
            [ 1  0  0  0  0  0]
            [ 0  1  0  0  1 -1]
            [ 0  0  1  0 -1  1]
            [ 0  0  0  1  0 -1]
        """
        return self.__submodule

    def hecke_bound(self):
        """
        Return a bound on the number of Hecke operators needed to generate
        the Hecke algebra acting on this homology group.

        EXAMPLES::

            sage: d = J0(43).homology().decomposition(2); d
            [
            Submodule of rank 2 of Integral Homology of Abelian variety J0(43) of dimension 3,
            Submodule of rank 4 of Integral Homology of Abelian variety J0(43) of dimension 3
            ]

        Because the first factor has dimension 2 it corresponds to an
        elliptic curve, so we have a Hecke bound of 1.

        ::

            sage: d[0].hecke_bound()
            1
            sage: d[1].hecke_bound()
            8
        """
        if self.rank() <= 2:
            return ZZ(1)
        else:
            return self.__ambient.hecke_bound()

    def hecke_matrix(self, n):
        """
        Return the matrix of the `n`-th Hecke operator acting on this
        homology group.

        EXAMPLES::

            sage: d = J0(125).homology(GF(17)).decomposition(2); d
            [
            Submodule of rank 4 of Homology with coefficients in Finite Field of size 17
             of Abelian variety J0(125) of dimension 8,
            Submodule of rank 4 of Homology with coefficients in Finite Field of size 17
             of Abelian variety J0(125) of dimension 8,
            Submodule of rank 8 of Homology with coefficients in Finite Field of size 17
             of Abelian variety J0(125) of dimension 8
            ]
            sage: t = d[0].hecke_matrix(17); t
            [16 15 15  0]
            [ 0  5  0  2]
            [ 2  0  5 15]
            [ 0 15  0 16]
            sage: t.base_ring()
            Finite Field of size 17
            sage: t.fcp()
            (x^2 + 13*x + 16)^2
        """
        n = Integer(n)
        try:
            return self.__hecke_matrix[n]
        except AttributeError:
            self.__hecke_matrix = {}
        except KeyError:
            pass
        t = self.__ambient.hecke_matrix(n)
        s = t.restrict(self.__submodule)
        self.__hecke_matrix[n] = s
        return s

    def rank(self):
        """
        Return the rank of this homology group.

        EXAMPLES::

            sage: d = J0(43).homology().decomposition(2)
            sage: [H.rank() for H in d]
            [2, 4]
        """
        return self.__submodule.rank()
