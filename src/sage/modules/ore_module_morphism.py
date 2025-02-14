r"""
Morphisms between Ore modules

Let `R` be a commutative ring, `\theta : K \to K` by a ring
endomorphism and `\partial : K \to K` be a `\theta`-derivation.
Let also `\mathcal S = K[X; \theta, \partial]` denote the
associated Ore polynomial ring.

By definition, a Ore module is a module over `\mathcal S`. In
SageMath, there are rather represented as modules over `R`
equipped with the map giving the action of the Ore variable `X`.
We refer to :mod:`sage.modules.ore_module` for more details.

A morphism of Ore modules is a `R`-linear morphism commuting
with the Ore action, or equivalenty a `\mathcal S`-linear.

.. RUBRIC:: Construction of morphisms

There are several ways for creating Ore modules morphisms in SageMath.
First of all, one can use the method :meth:`hom`, passing to it the
matrix (in the canonical bases) of the morphism we want to build::

    sage: K.<z> = GF(5^3)
    sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
    sage: M.<e0,e1> = S.quotient_module(X^2 + X + z)

    sage: mat = matrix(2, 2, [z,     3*z^2 + z + 2,
    ....:                     z + 1,       4*z + 4])
    sage: f = M.hom(mat)
    sage: f
    Ore module endomorphism of Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5

Clearly, this method is not optimal: typing all the entries of the
defining matrix is long and is a potential source of errors.

Instead, one can use a dictionary encoding the values taken by the
morphism on a set of generators; the morphism is then automatically
prolonged by `\mathcal S`-linearity.
Actually here, `f` was just the multiplication by `X^3` on `M`.
We can then redefine it simply as follows::

    sage: g = M.hom({e0: X^3*e0})
    sage: g
    Ore module endomorphism of Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5

One can then recover the matrix by using the method :meth:`matrix`::

    sage: g.matrix()
    [            z 3*z^2 + z + 2]
    [        z + 1       4*z + 4]

Alternatively, one can use the method :meth:`multiplication_map`::

    sage: h = M.multiplication_map(X^3)
    sage: h
    Ore module endomorphism of Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: g == h
    True

Of course, the method :meth:`multiplication_map` also accepts
values in the base ring::

    sage: h = M.multiplication_map(2)
    sage: h
    Ore module endomorphism of Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: h.matrix()
    [2 0]
    [0 2]

Be careful that scalar multiplications do not always properly
define a morphism of Ore modules::

    sage: M.multiplication_map(z)
    Traceback (most recent call last):
    ...
    ValueError: does not define a morphism of Ore modules

.. RUBRIC: Kernels, images and related things

SageMath provides methods to compute kernels, cokernels,
images and coimages. In order to illustrate this, we will
build the sequence

.. MATH::

    0 \to \mathcal S/ \mathcal S P
    \to \mathcal S/ \mathcal S PQ
    \to \mathcal S/ \mathcal S Q \to 0

and check that it is exact.
We first build the Ore modules::

    sage: P = X^2 + z*X + 1
    sage: Q = X^3 + z^2*X^2 + X + z
    sage: U = S.quotient_module(P, names='u')
    sage: U.inject_variables()
    Defining u0, u1
    sage: V = S.quotient_module(P*Q, names='v')
    sage: V.inject_variables()
    Defining v0, v1, v2, v3, v4
    sage: W = S.quotient_module(Q, names='w')
    sage: W.inject_variables()
    Defining w0, w1, w2

Next, we build the morphisms::

    sage: f = U.hom({u0: Q*v0})
    sage: g = V.hom({v0: w0})

We can now check that `f` is injective by computing its
kernel::

    sage: f.kernel()
    Ore module of rank 0 over Finite Field in z of size 5^3 twisted by z |--> z^5

We see on the output that it has dimension `0`; so it
vanishes. Instead of reading the output, one can check
programmatically the vanishing of a Ore module using the
method :meth:`is_zero`::

    sage: f.kernel().is_zero()
    True

Actually, in our use case, one can, more simply, use the
method :meth:`is_injective`::

    sage: f.is_injective()
    True

Similarly, one checks that `g` is surjective::

    sage: g.is_surjective()
    True

or equivalently::

    sage: g.cokernel().is_zero()
    True

Now, we need to check that the kernel of `g` equals the
image of `f`. For this, we compute both and compare the
results::

    sage: ker = g.kernel()
    sage: im = f.image()
    sage: ker == im
    True

As a sanity check, one can also verity that the composite
`g \circ f` vanishes::

    sage: h = g * f
    sage: h
    Ore module morphism:
      From: Ore module <u0, u1> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module <w0, w1, w2> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: h.is_zero()
    True

Let us now consider another morphism `f` and build the
canonical isomorphism `\text{coim }f \to \text{im }f`
that it induces.
We start by defining `f`::

    sage: A = X + z
    sage: B = X + z + 1
    sage: P = X^2 + X + z
    sage: U = S.quotient_module(B*P, names='u')
    sage: U.inject_variables()
    Defining u0, u1, u2
    sage: V = S.quotient_module(P*A, names='v')
    sage: V.inject_variables()
    Defining v0, v1, v2

    sage: f = U.hom({u0: A*v0})
    sage: f
    Ore module morphism:
      From: Ore module <u0, u1, u2> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module <v0, v1, v2> over Finite Field in z of size 5^3 twisted by z |--> z^5

Now we compute the image and the coimage::

    sage: I = f.image(names='im')
    sage: I
    Ore module <im0, im1> over Finite Field in z of size 5^3 twisted by z |--> z^5

    sage: C = f.coimage(names='co')
    sage: C
    Ore module <co0, co1> over Finite Field in z of size 5^3 twisted by z |--> z^5

We can already check that the image and the coimage have the
same rank. We now want to construct the isomorphism between
them. For this, we first need to corestrict `f` to its image.
This is achieved via the method :meth:`morphism_corestriction`
of the Ore module::

    sage: g = I.morphism_corestriction(f)
    sage: g
    Ore module morphism:
      From: Ore module <u0, u1, u2> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module <im0, im1> over Finite Field in z of size 5^3 twisted by z |--> z^5

Next, we want to factor `g` by the coimage. We proceed as follows::

    sage: h = C.morphism_quotient(g)
    sage: h
    Ore module morphism:
      From: Ore module <co0, co1> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module <im0, im1> over Finite Field in z of size 5^3 twisted by z |--> z^5

We have found the morphism we were looking for: it is `h`.
We can now check that it is an isomorphism::

    sage: h.is_isomorphism()
    True

As a shortcut, we can use explicit conversions as follows::

    sage: H = Hom(C, I)  # the hom space
    sage: h2 = H(f)
    sage: h2
    Ore module morphism:
      From: Ore module <co0, co1> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module <im0, im1> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: h == h2
    True

.. RUBRIC: Determinants and characteristic polynomials

For endomorphisms, one can compute classical invariants as
determinants and characteristic polynomials.
To illustrate this, we check on an example that the characteristic
polynomial of the multiplication by `X^3` on the quotient
`\mathcal S / \mathcal S P` is the reduced norm of `P`::

    sage: P = X^5 + z*X^4 + z^2*X^2 + z + 1
    sage: M = S.quotient_module(P)
    sage: f = M.multiplication_map(X^3)
    sage: f.charpoly()
    x^5 + x^4 + x^3 + x^2 + 1

    sage: P.reduced_norm('x')
    x^5 + x^4 + x^3 + x^2 + 1

AUTHOR:

- Xavier Caruso (2024-10)
"""

# ***************************************************************************
#    Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.misc.latex import latex
from sage.structure.element import Element
from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.categories.map import Map
from sage.categories.morphism import Morphism
from sage.modules.ore_module import OreModule, OreSubmodule, OreQuotientModule

class OreModuleMorphism(Morphism):
    r"""
    Generic class for morphism between Ore modules.
    """
    def __init__(self, parent, im_gens, check=True):
        r"""
        Initialize this Ore module.

        INPUT:

        - ``parent`` -- the hom space

        - ``im_gens`` -- the image of the generators (formatted as
          a list, a tuple, a dictionary or a matrix) or a Ore modules
          morphism

        - ``check`` (default: ``True``) -- a boolean, whether we
          should check if the given data correctly defined a morphism
          of Ore modules

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: type(f)
            <class 'sage.modules.ore_module_homspace.OreModule_homspace_with_category.element_class'>

            sage: TestSuite(f).run()
        """
        Morphism.__init__(self, parent)
        domain = parent.domain()
        codomain = parent.codomain()
        base = domain.base_ring()
        MS = parent.matrix_space()
        if (isinstance(im_gens, Element)
            and base.has_coerce_map_from(im_gens.parent())):
            self._matrix = MS(im_gens)
        elif isinstance(im_gens, Matrix):
            self._matrix = MS(im_gens)
        elif isinstance(im_gens, OreModuleMorphism):
            # Not optimal: too many intermediate morphisms constructed
            f = im_gens
            if f.domain() is not domain:
                f = domain._hom_change_domain(f)
            if f.codomain() is not codomain:
                f = codomain._hom_change_codomain(f)
            self._matrix = f._matrix
        elif isinstance(im_gens, (list, tuple)):
            if len(im_gens) != domain.rank():
                raise ValueError("wrong number of generators")
            self._matrix = matrix([codomain(v).list() for v in im_gens])
        elif isinstance(im_gens, dict):
            zero = parent.base_ring().zero()
            dimd = domain.rank()
            dimc = codomain.rank()
            d = dimc + dimd
            vs = [domain(x).list() + codomain(y).list() for x, y in im_gens.items()]
            if len(vs) < 2*d:
                vs += (2*d - len(vs)) * [d * [zero]]
            M = matrix(vs)
            M.echelonize()
            oldr = 0
            r = M.rank()
            iter = 1
            fd = domain._pseudohom
            fc = codomain._pseudohom
            while r > oldr:
                for i in range(r):
                    row = M.row(i).list()
                    x = row[:dimd]
                    y = row[dimd:]
                    for _ in range(iter):
                        x = fd(x)
                        y = fc(y)
                    v = x.list() + y.list()
                    for j in range(d):
                        M[i+r,j] = v[j]
                M.echelonize()
                oldr = r
                r = M.rank()
                iter *= 2
            if list(M.pivots()) != list(range(dimd)):
                raise ValueError("does not define a morphism of Ore modules")
            self._matrix = M.submatrix(0, dimd, dimd, dimc)
        else:
            raise ValueError("cannot construct a morphism from the given data")
        if check:
            for x in parent.domain().basis():
                if self._call_(x.image()) != self._call_(x).image():
                    raise ValueError("does not define a morphism of Ore modules")

    def _repr_type(self):
        r"""
        Return a string with the type of this morphism.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: f._repr_type()
            'Ore module'
        """
        return "Ore module"

    def _latex_(self):
        r"""
        Return a LaTeX representation of this morphism.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: latex(f)
            \begin{array}{l}
            \text{\texttt{Ore module morphism:}} \\
            \text{\texttt{{ }{ }From:}}\hspace{1ex} \texttt{Ore module of rank } 2\texttt{ over } \Bold{F}_{5^{3}} \texttt{ twisted by } z \mapsto z^{5} \\
            \text{\texttt{{ }{ }To:}}\hspace{3ex} \texttt{Ore module of rank } 2\texttt{ over } \Bold{F}_{5^{3}} \texttt{ twisted by } z \mapsto z^{5}
            \end{array}

        ::

            sage: Me = M.rename_basis('e')
            sage: fe = Me.multiplication_map(X^3)
            sage: latex(fe)
            \begin{array}{l}
            \text{\texttt{Ore module morphism:}} \\
            \text{\texttt{{ }{ }From:}}\hspace{1ex} \left<e_{0}, e_{1}\right>_{\Bold{F}_{5^{3}} , z \mapsto z^{5} } \\
            \text{\texttt{{ }{ }To:}}\hspace{3ex} \left<e_{0}, e_{1}\right>_{\Bold{F}_{5^{3}} , z \mapsto z^{5} }
            \end{array}
        """
        s = "\\begin{array}{l}\n"
        s += "\\text{\\texttt{%s morphism:}} \\\\\n" % self._repr_type()
        s += "\\text{\\texttt{{ }{ }From:}}\\hspace{1ex} %s \\\\\n" % latex(self.domain())
        s += "\\text{\\texttt{{ }{ }To:}}\\hspace{3ex} %s\n" % latex(self.codomain())
        s += "\\end{array}"
        return s

    def matrix(self):
        r"""
        Return the matrix defining this morphism.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(2)
            sage: f.matrix()
            [2 0]
            [0 2]

            sage: g = M.multiplication_map(X^3)
            sage: g.matrix()
            [            0 3*z^2 + z + 1]
            [      2*z + 1             0]
        """
        return self._matrix.__copy__()

    def _call_(self, x):
        r"""
        Return the image of `x` by this morphism.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module(X^2 + z*X + 1)
            sage: f = M.multiplication_map(X^3)
            sage: f(v)
            (2*z^2 + 4*z + 4)*v + (4*z^2 + 3*z + 3)*w

        We check that it is the correct answer::

            sage: X^3 * v
            (2*z^2 + 4*z + 4)*v + (4*z^2 + 3*z + 3)*w
        """
        return self.codomain()(x * self._matrix)

    def is_zero(self):
        r"""
        Return ``True`` if this morphism is zero.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2 + 1
            sage: M = S.quotient_module(P^2, names='e')
            sage: M.inject_variables()
            Defining e0, e1, e2, e3, e4, e5

            sage: f = M.hom({e0: P*e0})
            sage: f.is_zero()
            False
            sage: (f*f).is_zero()
            True
        """
        return self._matrix.is_zero()

    def is_identity(self):
        r"""
        Return ``True`` if this morphism is the identity.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module(X^2 + z)
            sage: f = M.hom({v: v})
            sage: f.is_identity()
            True

            sage: f = M.hom({v: 2*v})
            sage: f.is_identity()
            False
        """
        return self.domain() is self.codomain() and self._matrix.is_one()

    def _add_(self, other):
        r"""
        Return the sum of this morphism and ``other``.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: g = M.multiplication_map(X^6)
            sage: h = f + g
            sage: h == M.multiplication_map(X^3 + X^6)
            True
        """
        if not isinstance(other, OreModuleMorphism):
            raise ValueError("the morphism is not a morphism of Ore modules")
        H = self.parent()
        return H(self._matrix + other._matrix, check=False)

    def _neg_(self):
        r"""
        Return the oppositive of this morphism.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: g = -f
            sage: g == M.multiplication_map(-X^3)
            True
        """
        H = self.parent()
        return H(-self._matrix, check=False)

    def _sub_(self, other):
        r"""
        Return the different between this morphism and ``other``.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: g = M.multiplication_map(X^6)
            sage: h = f - g
            sage: h == M.multiplication_map(X^3 - X^6)
            True
        """
        if not isinstance(other, OreModuleMorphism):
            raise ValueError("the morphism is not a morphism of Ore modules")
        H = self.parent()
        return H(self._matrix - other._matrix, check=False)

    def _rmul_(self, a):
        r"""
        Return the product of the scalar `a` by this morphism.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: g = 2*f
            sage: g == M.multiplication_map(2*X^3)
            True
        """
        H = self.parent()
        return H(a*self._matrix, check=False)

    def __eq__(self, other):
        r"""
        Return ``True`` if this morphism is equal to ``other``.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = 2*M.multiplication_map(X^3)
            sage: g = M.multiplication_map(2*X^3)
            sage: f == g
            True
        """
        if not isinstance(other, OreModuleMorphism):
            try:
                other = self.parent()(other)
            except ValueError:
                return False
        return self._matrix == other._matrix

    def is_injective(self):
        r"""
        Return ``True`` if this morphism is injective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2
            sage: M = S.quotient_module(P^2, names='m')
            sage: M.inject_variables()
            Defining m0, m1, m2, m3, m4, m5
            sage: N = S.quotient_module(P, names='n')
            sage: N.inject_variables()
            Defining n0, n1, n2

            sage: f = N.hom({n0: P*m0})
            sage: f.is_injective()
            True

            sage: g = M.hom({m0: n0})
            sage: g.is_injective()
            False
        """
        return self._matrix.rank() == self.domain().rank()

    def is_surjective(self):
        r"""
        Return ``True`` if this morphism is surjective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2
            sage: M = S.quotient_module(P^2, names='m')
            sage: M.inject_variables()
            Defining m0, m1, m2, m3, m4, m5
            sage: N = S.quotient_module(P, names='n')
            sage: N.inject_variables()
            Defining n0, n1, n2

            sage: f = N.hom({n0: P*m0})
            sage: f.is_surjective()
            False

            sage: g = M.hom({m0: n0})
            sage: g.is_surjective()
            True
        """
        return self._matrix.rank() == self.codomain().rank()

    def is_bijective(self):
        r"""
        Return ``True`` if this morphism is bijective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: f.is_bijective()
            True

            sage: N = S.quotient_module(X^2)
            sage: g = N.multiplication_map(X^3)
            sage: g.is_bijective()
            False
        """
        return self.is_injective() and self.is_surjective()

    def is_isomorphism(self):
        r"""
        Return ``True`` if this morphism is an isomorphism.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: f.is_isomorphism()
            True

            sage: N = S.quotient_module(X^2)
            sage: g = N.multiplication_map(X^3)
            sage: g.is_isomorphism()
            False
        """
        return self.is_bijective()

    def _composition_(self, other, homset):
        r"""
        Return the composite ``other`` `\circ` ``self``.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: g = f * f
            sage: g == M.multiplication_map(X^6)
            True
        """
        if not isinstance(other, OreModuleMorphism):
            raise ValueError("the morphism is not a morphism of Ore modules")
        return homset(other._matrix * self._matrix, check=False)

    def inverse(self):
        r"""
        Return the inverse of this morphism.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^2 + z
            sage: M.<e0,e1,e2,e3> = S.quotient_module(P^2)

            sage: f = M.multiplication_map(X^3)
            sage: g = f.inverse()
            sage: (f*g).is_identity()
            True
            sage: (g*f).is_identity()
            True

        If the morphism is not invertible, an error is raised::

            sage: h = M.hom({e0: P*e0})
            sage: h.inverse()
            Traceback (most recent call last):
            ...
            ValueError: this morphism is not invertible
        """
        if not self.is_isomorphism():
            raise ValueError("this morphism is not invertible")
        H = self.parent()
        return H(self._matrix.inverse(), check=False)

    __invert__ = inverse

    def kernel(self, names=None):
        r"""
        Return ``True`` if this morphism is injective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2
            sage: M = S.quotient_module(P^2, names='m')
            sage: M.inject_variables()
            Defining m0, m1, m2, m3, m4, m5
            sage: N = S.quotient_module(P, names='n')
            sage: N.inject_variables()
            Defining n0, n1, n2

            sage: f = M.hom({m0: n0})
            sage: ker = f.kernel()
            sage: ker
            Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: ker.basis()
            [m0 + (2*z^2 + 3*z + 1)*m3 + (4*z^2 + 3*z + 3)*m4 + (2*z^2 + 3*z)*m5,
             m1 + (z + 3)*m3 + (z^2 + z + 4)*m4,
             m2 + (2*z^2 + 4*z + 2)*m4 + (2*z^2 + z + 1)*m5]
        """
        ker = self._matrix.left_kernel_matrix()
        return OreSubmodule(self.domain(), ker, names)

    def image(self, names=None):
        r"""
        Return ``True`` if this morphism is injective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2
            sage: M = S.quotient_module(P^2, names='m')
            sage: M.inject_variables()
            Defining m0, m1, m2, m3, m4, m5
            sage: N = S.quotient_module(P, names='n')
            sage: N.inject_variables()
            Defining n0, n1, n2

            sage: f = N.hom({n0: P*m0})
            sage: im = f.image()
            sage: im
            Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: im.basis()
            [m0 + (2*z^2 + 3*z + 1)*m3 + (4*z^2 + 3*z + 3)*m4 + (2*z^2 + 3*z)*m5,
             m1 + (z + 3)*m3 + (z^2 + z + 4)*m4,
             m2 + (2*z^2 + 4*z + 2)*m4 + (2*z^2 + z + 1)*m5]
        """
        return OreSubmodule(self.codomain(), self._matrix, names)

    def cokernel(self, names=None):
        r"""
        Return ``True`` if this morphism is injective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2
            sage: M = S.quotient_module(P^2, names='m')
            sage: M.inject_variables()
            Defining m0, m1, m2, m3, m4, m5
            sage: N = S.quotient_module(P, names='n')
            sage: N.inject_variables()
            Defining n0, n1, n2

            sage: f = N.hom({n0: P*m0})
            sage: coker = f.cokernel()
            sage: coker
            Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: coker.basis()
            [m3, m4, m5]
        """
        return OreQuotientModule(self.codomain(), self._matrix, names)

    def coimage(self, names=None):
        r"""
        Return ``True`` if this morphism is injective.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + z^2
            sage: M = S.quotient_module(P^2, names='m')
            sage: M.inject_variables()
            Defining m0, m1, m2, m3, m4, m5
            sage: N = S.quotient_module(P, names='n')
            sage: N.inject_variables()
            Defining n0, n1, n2

            sage: f = M.hom({m0: n0})
            sage: coim = f.coimage()
            sage: coim
            Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: coim.basis()
            [m3, m4, m5]
        """
        ker = self._matrix.left_kernel_matrix()
        return OreQuotientModule(self.domain(), ker, names)

    def determinant(self):
        r"""
        Return the determinant of this endomorphism.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<m0,m1> = S.quotient_module(X^2 + z)
            sage: f = M.multiplication_map(X^3)
            sage: f.determinant()
            2

        If the domain differs from the codomain (even if they have
        the same rank), an error is raised::

            sage: N.<n0,n1> = S.quotient_module(X^2 + z^25)
            sage: g = M.hom({z*m0: n0})
            sage: g.determinant()
            Traceback (most recent call last):
            ...
            ValueError: determinants are only defined for endomorphisms
        """
        if self.domain() is not self.codomain():
            raise ValueError("determinants are only defined for endomorphisms")
        return self._matrix.determinant()

    det = determinant

    def characteristic_polynomial(self, var='x'):
        r"""
        Return the determinant of this endomorphism.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 + (z^2 + 3)*X + z^5
            sage: M = S.quotient_module(P)
            sage: f = M.multiplication_map(X^3)
            sage: f.characteristic_polynomial()
            x^3 + x^2 + 2*x + 2

        We check that the latter is equal to the reduced norm
        of `P`::

            sage: P.reduced_norm('x')
            x^3 + x^2 + 2*x + 2

        TESTS::

            sage: M.<m0,m1> = S.quotient_module(X^2 + z)
            sage: N.<n0,n1> = S.quotient_module(X^2 + z^25)
            sage: g = M.hom({z*m0: n0})
            sage: g.characteristic_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: characteristic polynomials are only defined for endomorphisms
        """
        if self.domain() is not self.codomain():
            raise ValueError("characteristic polynomials are only defined for endomorphisms")
        return self._matrix.charpoly(var)

    charpoly = characteristic_polynomial

class OreModuleRetraction(Map):
    r"""
    Conversion (partially defined) map from an ambient module
    to one of its submodule.
    """
    def _call_(self, y):
        r"""
        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^2 + z*X + 1
            sage: M = S.quotient_module(P^2, names='e')
            sage: M.inject_variables()
            Defining e0, e1, e2, e3
            sage: N = M.span(P*e0, names='u')
            sage: N(P*e0)  # indirect doctest
            u0 + z*u1

            sage: N(e0)
            Traceback (most recent call last):
            ...
            ValueError: not in the submodule
        """
        X = self.codomain()
        try:
            xs = X._basis.solve_left(y)
        except ValueError:
            raise ValueError("not in the submodule")
        return X(xs)

class OreModuleSection(Map):
    r"""
    Section map of the projection onto a quotient.
    It is not necessarily compatible with the Ore action.
    """
    def _call_(self, y):
        r"""
        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^2 + z*X + 1
            sage: M = S.quotient_module(P^2, names='e')
            sage: M.inject_variables()
            Defining e0, e1, e2, e3
            sage: N = M.quo(P*e0, names='u')
            sage: N.inject_variables()
            Defining u0, u1
            sage: M(u0)  # indirect doctest
            e2
        """
        X = self.codomain()
        Y = self.domain()
        indices = Y._indices
        zero = X.base_ring().zero()
        xs = X.rank() * [zero]
        for i in range(Y.rank()):
            xs[indices[i]] = y[i]
        return X(xs)
