# sage.doctest: needs sage.libs.singular
r"""
Free resolutions

Let `R` be a commutative ring. A finite free resolution of an `R`-module `M`
is a chain complex of free `R`-modules

.. MATH::

    0 \rightarrow R^{n_k} \xrightarrow{d_k}
    \cdots \xrightarrow{d_2} R^{n_1} \xrightarrow{d_1} R^{n_0}

terminating with a zero module at the end that is exact (all homology groups
are zero) such that the image of `d_1` is `M`.

EXAMPLES::

    sage: from sage.homology.free_resolution import FreeResolution
    sage: S.<x,y,z,w> = PolynomialRing(QQ)
    sage: m = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z]).transpose()
    sage: r = FreeResolution(m, name='S'); r
    S^1 <-- S^3 <-- S^2 <-- 0

    sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
    sage: r = I.free_resolution(); r
    S^1 <-- S^3 <-- S^2 <-- 0

::

    sage: S.<x,y,z,w> = PolynomialRing(QQ)
    sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
    sage: r = I.graded_free_resolution(); r
    S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0

An example of a minimal free resolution from [CLO2005]_::

    sage: R.<x,y,z,w> = QQ[]
    sage: I = R.ideal([y*z - x*w, y^3 - x^2*z, x*z^2 - y^2*w, z^3 - y*w^2])
    sage: r = I.free_resolution();  r
    S^1 <-- S^4 <-- S^4 <-- S^1 <-- 0
    sage: len(r)
    3
    sage: r.matrix(2)
    [-z^2 -x*z  y*w -y^2]
    [   y    0   -x    0]
    [  -w    y    z    x]
    [   0    w    0    z]

AUTHORS:

- Kwankyu Lee (2022-05-13): initial version
- Travis Scrimshaw (2022-08-23): refactored for free module inputs
"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#                 (C) 2022 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.abstract_method import abstract_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.sage_object import SageObject
from sage.structure.element import Matrix
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.integral_domains import IntegralDomains
from sage.modules.free_module_element import vector
from sage.modules.free_module import FreeModule
from sage.modules.free_module import Module_free_ambient, FreeModule_generic
from sage.rings.ideal import Ideal_generic

from copy import copy


class FreeResolution(SageObject, metaclass=ClasscallMetaclass):
    r"""
    A free resolution.

    Let `R` be a commutative ring. A *free resolution* of an `R`-module `M`
    is a (possibly infinite) chain complex of free `R`-modules

    .. MATH::

        \cdots \rightarrow R^{n_k} \xrightarrow{d_k}
        \cdots \xrightarrow{d_2} R^{n_1} \xrightarrow{d_1} R^{n_0}

    that is exact (all homology groups are zero) such that the image
    of `d_1` is `M`.
    """
    @staticmethod
    def __classcall_private__(cls, module, *args, graded=False, degrees=None, shifts=None, **kwds):
        """
        Dispatch to the correct constructor.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: m = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z]).transpose()
            sage: r = FreeResolution(m, name='S')
            sage: type(r)
            <class 'sage.homology.free_resolution.FiniteFreeResolution_singular'>

            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(I)
            sage: type(r)
            <class 'sage.homology.free_resolution.FiniteFreeResolution_singular'>

            sage: R.<x> = QQ[]
            sage: M = R^3
            sage: v = M([x^2, 2*x^2, 3*x^2])
            sage: w = M([0, x, 2*x])
            sage: S = M.submodule([v, w])
            sage: r = FreeResolution(S)
            sage: type(r)
            <class 'sage.homology.free_resolution.FiniteFreeResolution_free_module'>

            sage: I = R.ideal([x^4 + 3*x^2 + 2])
            sage: r = FreeResolution(I)
            sage: type(r)
            <class 'sage.homology.free_resolution.FiniteFreeResolution_free_module'>

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([x^2, y^3])
            sage: Q = R.quo(I)
            sage: Q.is_integral_domain()
            False
            sage: xb, yb = Q.gens()                                                     # needs sage.rings.function_field
            sage: FreeResolution(Q.ideal([xb]))  # has torsion                          # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            NotImplementedError: the ring must be a polynomial ring using Singular
        """
        if degrees is not None or shifts is not None:
            graded = True

        if isinstance(module, Matrix):
            is_free_module = False
            S = module.base_ring()
            if S in PrincipalIdealDomains():
                module = module.echelon_form()
                if module.nrows() > module.rank():
                    module = module.submatrix(nrows=module.rank())
                    module.set_immutable()
                is_free_module = True
            if not module.is_immutable():
                # We need to make an immutable copy of the matrix
                module = copy(module)
                module.set_immutable()
            if is_free_module:
                if graded:
                    from sage.homology.graded_resolution import GradedFiniteFreeResolution_free_module
                    return GradedFiniteFreeResolution_free_module(module,
                                                                  *args,
                                                                  degrees=degrees,
                                                                  shifts=shifts,
                                                                  **kwds)
                return FiniteFreeResolution_free_module(module, *args, **kwds)

            from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            if not isinstance(S, MPolynomialRing_libsingular):
                raise NotImplementedError("the matrix must be over a PID or a "
                                          " polynomial ring that is using Singular")

            if graded:
                # We are computing a graded resolution
                from sage.homology.graded_resolution import GradedFiniteFreeResolution_singular
                return GradedFiniteFreeResolution_singular(module, *args, degrees=degrees,
                                                           shifts=shifts, **kwds)

            return FiniteFreeResolution_singular(module, **kwds)

        if graded:
            return module.graded_free_resolution(*args, **kwds)
        return module.free_resolution(*args, **kwds)

    def __init__(self, module, name='S', **kwds):
        """
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a ring
        - ``name`` -- (default: ``'S'``) the name of the ring for printing

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: m1 = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z])
            sage: r = FreeResolution(m1, name='S')
            sage: TestSuite(r).run(skip=['_test_pickling'])
        """
        if isinstance(module, Ideal_generic):
            S = module.ring()
        else:  # module or matrix
            S = module.base_ring()

        self._base_ring = S
        self._name = name
        self._module = module

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: m1 = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z])
            sage: r = FreeResolution(m1, name='S')
            sage: print(FreeResolution._repr_(r))
            Free resolution of the row space of the matrix:
            [z^2 - y*w y*z - x*w y^2 - x*z]
        """
        if isinstance(self._module, Matrix):
            return f"Free resolution of the row space of the matrix:\n{self._module}"
        return f"Free resolution of {self._module}"

    def _repr_module(self, i):
        r"""
        Return the string form of the `i`-th free module.

        INPUT:

        - ``i`` -- positive integer

        EXAMPLES::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: m = matrix(S, 1, [y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(m.transpose(), name='S')
            sage: r._repr_module(2)
            'S^2'
            sage: r  # indirect doctest
            S^1 <-- S^3 <-- S^2 <-- 0

        TESTS::

            sage: S.<x,y,z> = PolynomialRing(QQ)
            sage: I = S.ideal(0)
            sage: C = I.free_resolution()
            sage: C
            S^1 <-- 0
        """
        if i == 0:
            if self._length > 0:
                r = self._maps[0].nrows()
            else:
                r = self._initial_differential.domain().dimension()
            s = f'{self._name}^{r}'
            return s
        elif i > self._length:
            s = '0'
        else:
            r = self._maps[i - 1].ncols()
            if r > 0:
                s = f'{self._name}^{r}'
            else:
                s = '0'
        return s

    @abstract_method
    def differential(self, i):
        r"""
        Return the `i`-th differential map.

        INPUT:

        - ``i`` -- positive integer

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: m1 = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z])
            sage: r = FreeResolution(m1, name='S')
            sage: FreeResolution.differiental(r, 1)
            Traceback (most recent call last):
            ...
            AttributeError: type object 'FreeResolution' has no attribute 'differiental'...
        """

    def target(self):
        r"""
        Return the codomain of the `0`-th differential map.

        The codomain of the `0`-th differential map is the cokernel of
        the first differential map.

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r
            S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
            sage: r.target()
            Quotient module by
             Submodule of Ambient free module of rank 1 over the integral domain
              Multivariate Polynomial Ring in x, y, z, w over Rational Field
              Generated by the rows of the matrix:
              [-z^2 + y*w]
              [ y*z - x*w]
              [-y^2 + x*z]
        """
        return self.differential(0).codomain()


class FiniteFreeResolution(FreeResolution):
    r"""
    Finite free resolutions.

    The matrix at index `i` in the list defines the differential map from
    `(i + 1)`-th free module to the `i`-th free module over the base ring by
    multiplication on the left. The number of matrices in the list is the
    length of the resolution. The number of rows and columns of the matrices
    define the ranks of the free modules in the resolution.

    Note that the first matrix in the list defines the differential map at
    homological index `1`.

    A subclass must provide a ``_maps`` attribute that contains a list of the
    maps defining the resolution.

    A subclass can define ``_initial_differential`` attribute that
    contains the `0`-th differential map whose codomain is the target
    of the free resolution.

    EXAMPLES::

        sage: from sage.homology.free_resolution import FreeResolution
        sage: S.<x,y,z,w> = PolynomialRing(QQ)
        sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
        sage: r = FreeResolution(I)
        sage: r.differential(0)
        Coercion map:
          From: Ambient free module of rank 1 over the integral domain
                Multivariate Polynomial Ring in x, y, z, w over Rational Field
          To:   Quotient module by
                Submodule of Ambient free module of rank 1 over the integral domain
                Multivariate Polynomial Ring in x, y, z, w over Rational Field
                Generated by the rows of the matrix:
                [-z^2 + y*w]
                [ y*z - x*w]
                [-y^2 + x*z]
    """
    @lazy_attribute
    def _length(self):
        """
        The length of ``self``.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(I)
            sage: r._length
            2
        """
        return len(self._maps)

    def _repr_(self):
        """
        Return the string form of this resolution.

        INPUT:

        - ``i`` -- positive integer

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution(); r
            S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
        """
        s = self._repr_module(0)
        for i in range(1, self._length + 1):
            s += ' <-- ' + self._repr_module(i)
        s += ' <-- 0'
        return s

    def __len__(self):
        r"""
        Return the length of this resolution.

        The length of a free resolution is the index of the last nonzero free module.

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution(); r
            S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
            sage: len(r)
            2
        """
        return len(self._maps)

    def __getitem__(self, i):
        r"""
        Return the `i`-th free module of this resolution.

        INPUT:

        - ``i`` -- positive integer

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution(); r
            S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
            sage: r.target()
            Quotient module by Submodule of Ambient free module of rank 1 over the integral domain
            Multivariate Polynomial Ring in x, y, z, w over Rational Field
            Generated by the rows of the matrix:
            [-z^2 + y*w]
            [ y*z - x*w]
            [-y^2 + x*z]
        """
        if i < 0:
            raise IndexError('invalid index')
        elif i > self._length:
            F = FreeModule(self._base_ring, 0)
        elif i == 0:
            F = self.differential(0).domain()
        elif i == self._length:
            F = FreeModule(self._base_ring, self._maps[i - 1].ncols())
        else:
            F = FreeModule(self._base_ring, self._maps[i].nrows())
        return F

    def differential(self, i):
        r"""
        Return the `i`-th differential map.

        INPUT:

        - ``i`` -- positive integer

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r
            S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
            sage: r.differential(3)
            Free module morphism defined as left-multiplication by the matrix
             []
             Domain:   Ambient free module of rank 0 over the integral domain
                       Multivariate Polynomial Ring in x, y, z, w over Rational Field
             Codomain: Ambient free module of rank 2 over the integral domain
                       Multivariate Polynomial Ring in x, y, z, w over Rational Field
            sage: r.differential(2)
            Free module morphism defined as left-multiplication by the matrix
              [-y  x]
              [ z -y]
              [-w  z]
              Domain:   Ambient free module of rank 2 over the integral domain
                        Multivariate Polynomial Ring in x, y, z, w over Rational Field
              Codomain: Ambient free module of rank 3 over the integral domain
                        Multivariate Polynomial Ring in x, y, z, w over Rational Field
            sage: r.differential(1)
            Free module morphism defined as left-multiplication by the matrix
              [z^2 - y*w y*z - x*w y^2 - x*z]
              Domain:   Ambient free module of rank 3 over the integral domain
                        Multivariate Polynomial Ring in x, y, z, w over Rational Field
              Codomain: Ambient free module of rank 1 over the integral domain
                        Multivariate Polynomial Ring in x, y, z, w over Rational Field
            sage: r.differential(0)
            Coercion map:
              From: Ambient free module of rank 1 over the integral domain
                    Multivariate Polynomial Ring in x, y, z, w over Rational Field
              To:   Quotient module by
                    Submodule of Ambient free module of rank 1 over the integral domain
                    Multivariate Polynomial Ring in x, y, z, w over Rational Field
                    Generated by the rows of the matrix:
                    [-z^2 + y*w]
                    [ y*z - x*w]
                    [-y^2 + x*z]

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: S = P2.coordinate_ring()
            sage: I = S.ideal(0)
            sage: C = I.graded_free_resolution(); C
            S(0) <-- 0
            sage: C[1]
            Ambient free module of rank 0 over the integral domain
             Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: C[0]
            Ambient free module of rank 1 over the integral domain
             Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: C.differential(1)
            Free module morphism defined as left-multiplication by the matrix
            []
            Domain: Ambient free module of rank 0 over the integral domain
             Multivariate Polynomial Ring in x, y, z over Rational Field
            Codomain: Ambient free module of rank 1 over the integral domain
             Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: C.differential(1).matrix()
            []
            sage: C.differential(1).matrix().dimensions()
            (1, 0)
        """
        if i < 0:
            raise IndexError('invalid index')
        elif i == 0:
            try:
                return self._initial_differential
            except AttributeError:
                raise ValueError('0th differential map undefined')
        elif i > self._length + 1:
            s = FreeModule(self._base_ring, 0)
            t = FreeModule(self._base_ring, 0)
            m = s.hom(0, t, side='right')
        elif i == self._length + 1:
            s = FreeModule(self._base_ring, 0)
            if self._length > 0:
                t = FreeModule(self._base_ring, self._maps[i - 2].ncols())
            else:
                t = self._initial_differential.domain()
            m = s.hom(0, t, side='right')
        else:
            s = FreeModule(self._base_ring, self._maps[i - 1].ncols())
            t = FreeModule(self._base_ring, self._maps[i - 1].nrows())
            m = s.hom(self._maps[i - 1], t, side='right')
        return m

    def matrix(self, i):
        r"""
        Return the matrix representing the `i`-th differential map.

        INPUT:

        - ``i`` -- positive integer

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution(); r
            S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
            sage: r.matrix(3)
            []
            sage: r.matrix(2)
            [-y  x]
            [ z -y]
            [-w  z]
            sage: r.matrix(1)
            [z^2 - y*w y*z - x*w y^2 - x*z]
        """
        if i <= 0:
            raise IndexError('invalid index')
        elif i <= self._length:
            return self._maps[i - 1]
        else:
            return self.differential(i).matrix()

    def chain_complex(self):
        r"""
        Return this resolution as a chain complex.

        A chain complex in Sage has its own useful methods.

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: unicode_art(r.chain_complex())
                                                           ⎛-y  x⎞
                                                           ⎜ z -y⎟
                       (z^2 - y*w y*z - x*w y^2 - x*z)     ⎝-w  z⎠
             0 <── C_0 <────────────────────────────── C_1 <────── C_2 <── 0
        """
        from sage.homology.chain_complex import ChainComplex
        mats = {}
        for i in range(self._length, 0, -1):
            mats[i] = self.matrix(i)
        return ChainComplex(mats, degree_of_differential=-1)

    @lazy_attribute
    def _initial_differential(self):
        r"""
        Define the `0`-th differential map of this resolution.

        EXAMPLES::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(I)
            sage: r._initial_differential
            Coercion map:
              From: Ambient free module of rank 1 over the integral domain
                    Multivariate Polynomial Ring in x, y, z, w over Rational Field
              To:   Quotient module by
                    Submodule of Ambient free module of rank 1 over the integral domain
                    Multivariate Polynomial Ring in x, y, z, w over Rational Field
            Generated by the rows of the matrix:
            [-z^2 + y*w]
            [ y*z - x*w]
            [-y^2 + x*z]
        """
        module = self._module
        if isinstance(module, Ideal_generic):
            S = module.ring()
            M = FreeModule(S, 1)
            N = M.submodule([vector([g]) for g in module.gens()])
        elif isinstance(module, Module_free_ambient):
            S = module.base_ring()
            M = module.ambient_module()
            N = module
        elif isinstance(module, Matrix):
            S = module.base_ring()
            N = module.row_space()
            M = N.ambient_module()
        Q = M.quotient(N)
        return Q.coerce_map_from(M)

    def _m(self):
        r"""
        Return the matrix whose column space is ``self._module``.

        If ``self._module`` is an ideal, then just the ideal is returned.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(I)
            sage: r._m()
            Ideal (-z^2 + y*w, y*z - x*w, -y^2 + x*z) of
             Multivariate Polynomial Ring in x, y, z, w over Rational Field

            sage: m = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z]).transpose()
            sage: r = FreeResolution(m, name='S')
            sage: r._m()
            [z^2 - y*w y*z - x*w y^2 - x*z]

            sage: M = m.image()
            sage: r = FreeResolution(M, name='S')
            sage: r._m()
            [z^2 - y*w y*z - x*w y^2 - x*z]
        """
        if isinstance(self._module, Ideal_generic):
            return self._module
        if isinstance(self._module, Module_free_ambient):
            return self._module.matrix().transpose()
        if isinstance(self._module, Matrix):
            return self._module.transpose()
        raise ValueError("unable to create a matrix/ideal to build the resolution")


class FiniteFreeResolution_free_module(FiniteFreeResolution):
    r"""
    Free resolutions of a free module.

    INPUT:

    - ``module`` -- a free module or ideal over a PID
    - ``name`` -- the name of the base ring

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: M = R^3
        sage: v = M([x^2, 2*x^2, 3*x^2])
        sage: w = M([0, x, 2*x])
        sage: S = M.submodule([v, w]); S
        Free module of degree 3 and rank 2 over
         Univariate Polynomial Ring in x over Rational Field
         Echelon basis matrix:
         [  x^2 2*x^2 3*x^2]
         [    0     x   2*x]
        sage: res = S.free_resolution(); res
        S^3 <-- S^2 <-- 0
        sage: ascii_art(res.chain_complex())
                    [  x^2     0]
                    [2*x^2     x]
                    [3*x^2   2*x]
         0 <-- C_0 <-------------- C_1 <-- 0

        sage: R.<x> = PolynomialRing(QQ)
        sage: I = R.ideal([x^4 + 3*x^2 + 2])
        sage: res = I.free_resolution(); res
        S^1 <-- S^1 <-- 0
    """
    @lazy_attribute
    def _maps(self):
        r"""
        Return the maps that define ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: M = R^3
            sage: v = M([x^2, 2*x^2, 3*x^2])
            sage: w = M([0, x, 2*x])
            sage: S = M.submodule([v, w])
            sage: res = S.free_resolution(); res
            S^3 <-- S^2 <-- 0
            sage: ascii_art(res.chain_complex())
                        [  x^2     0]
                        [2*x^2     x]
                        [3*x^2   2*x]
             0 <-- C_0 <-------------- C_1 <-- 0
            sage: res._maps
            [
            [  x^2     0]
            [2*x^2     x]
            [3*x^2   2*x]
            ]

            sage: R.<x> = PolynomialRing(QQ)
            sage: I = R.ideal([x^4 + 3*x^2 + 2])
            sage: res = I.free_resolution()
            sage: res._maps
            [[x^4 + 3*x^2 + 2]]

            sage: from sage.homology.free_resolution import FreeResolution
            sage: M = matrix([[x^2, 2],
            ....:             [3*x^2, 5],
            ....:             [5*x^2, 4]])
            sage: res = FreeResolution(M.transpose()); res
            S^3 <-- S^2 <-- 0
            sage: res._m()
            [     1      0]
            [   5/2   -x^2]
            [     2 -6*x^2]
            sage: res._maps
            [
            [     1      0]
            [   5/2   -x^2]
            [     2 -6*x^2]
            ]

        An overdetermined system over a PID::

            sage: res = FreeResolution(M); res
            S^2 <-- S^2 <-- 0
            sage: res._m()
            [x^2   0]
            [  2  -1]
            sage: res._maps
            [
            [x^2   0]
            [  2  -1]
            ]
        """
        if isinstance(self._module, Ideal_generic):
            from sage.matrix.constructor import matrix
            return [matrix([[self._module.gen()]])]
        return [self._m()]


class FiniteFreeResolution_singular(FiniteFreeResolution):
    r"""
    Minimal free resolutions of ideals or submodules of free modules
    of multivariate polynomial rings implemented in Singular.

    INPUT:

    - ``module`` -- a submodule of a free module `M` of rank `n` over `S` or
      an ideal of a multi-variate polynomial ring

    - ``name`` -- string (optional); name of the base ring

    - ``algorithm`` -- (default: ``'heuristic'``) Singular algorithm
      to compute a resolution of ``ideal``

    OUTPUT: a minimal free resolution of the ideal

    If ``module`` is an ideal of `S`, it is considered as a submodule of a
    free module of rank `1` over `S`.

    The available algorithms and the corresponding Singular commands
    are shown below:

    ============= ============================
    algorithm     Singular commands
    ============= ============================
    ``minimal``   ``mres(ideal)``
    ``shreyer``   ``minres(sres(std(ideal)))``
    ``standard``  ``minres(nres(std(ideal)))``
    ``heuristic`` ``minres(res(std(ideal)))``
    ============= ============================

    EXAMPLES::

        sage: from sage.homology.free_resolution import FreeResolution
        sage: S.<x,y,z,w> = PolynomialRing(QQ)
        sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
        sage: r = FreeResolution(I); r
        S^1 <-- S^3 <-- S^2 <-- 0
        sage: len(r)
        2

    ::

        sage: FreeResolution(I, algorithm='minimal')
        S^1 <-- S^3 <-- S^2 <-- 0
        sage: FreeResolution(I, algorithm='shreyer')
        S^1 <-- S^3 <-- S^2 <-- 0
        sage: FreeResolution(I, algorithm='standard')
        S^1 <-- S^3 <-- S^2 <-- 0
        sage: FreeResolution(I, algorithm='heuristic')
        S^1 <-- S^3 <-- S^2 <-- 0

    We can also construct a resolution by passing in a matrix defining
    the initial differential::

        sage: m = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z]).transpose()
        sage: r = FreeResolution(m, name='S'); r
        S^1 <-- S^3 <-- S^2 <-- 0
        sage: r.matrix(1)
        [z^2 - y*w y*z - x*w y^2 - x*z]

    An additional construction is using a submodule of a free module::

        sage: M = m.image()
        sage: r = FreeResolution(M, name='S'); r
        S^1 <-- S^3 <-- S^2 <-- 0

    A nonhomogeneous ideal::

        sage: I = S.ideal([z^2 - y*w, y*z - x*w, y^2 - x])
        sage: R = FreeResolution(I); R
        S^1 <-- S^3 <-- S^3 <-- S^1 <-- 0
        sage: R.matrix(2)
        [ y*z - x*w    y^2 - x          0]
        [-z^2 + y*w          0    y^2 - x]
        [         0 -z^2 + y*w -y*z + x*w]
        sage: R.matrix(3)
        [   y^2 - x]
        [-y*z + x*w]
        [ z^2 - y*w]
    """
    def __init__(self, module, name='S', algorithm='heuristic', **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(I)
            sage: TestSuite(r).run(skip=['_test_pickling'])

            sage: m = matrix(S, 1, [z^2 - y*w, y*z - x*w, y^2 - x*z]).transpose()
            sage: r = FreeResolution(m, name='S')
            sage: TestSuite(r).run(skip=['_test_pickling'])

            sage: M = m.image()
            sage: r = FreeResolution(M, name='S')
            sage: TestSuite(r).run(skip=['_test_pickling'])
        """
        self._algorithm = algorithm
        super().__init__(module, name=name, **kwds)

    @lazy_attribute
    def _maps(self):
        r"""
        Return the maps that define ``self``.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = FreeResolution(I)
            sage: r._maps
            [
                                             [-y  x]
                                             [ z -y]
            [z^2 - y*w y*z - x*w y^2 - x*z], [-w  z]
            ]
        """
        from sage.libs.singular.singular import si2sa_resolution
        from sage.libs.singular.function import singular_function

        # This ensures the first component of the Singular resolution to be a
        # module, like the later components. This is important when the
        # components are converted to Sage modules.
        module = singular_function("module")
        mod = module(self._m())

        if self._algorithm == 'minimal':
            mres = singular_function('mres')  # syzygy method
            r = mres(mod, 0)
        elif self._algorithm == 'shreyer':
            std = singular_function('std')
            sres = singular_function('sres')  # Shreyer method
            minres = singular_function('minres')
            r = minres(sres(std(mod), 0))
        elif self._algorithm == 'standard':
            nres = singular_function('nres')  # standard basis method
            minres = singular_function('minres')
            r = minres(nres(mod, 0))
        elif self._algorithm == 'heuristic':
            std = singular_function('std')
            res = singular_function('res')    # heuristic method
            minres = singular_function('minres')
            r = minres(res(std(mod), 0))

        return si2sa_resolution(r)
