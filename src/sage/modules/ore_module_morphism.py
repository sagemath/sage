r"""
Morphisms between Ore modules

By definition, a morphism of Ore module is a `R`-linear morphism
commuting with the Ore action, or equivalenty a `\mathcal S`-linear
map.

.. RUBRIC:: Construction of morphisms

There are several ways for creating Ore modules morphisms in SageMath.
First of all, one can use the method :meth:`hom`, passing in to it the
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
morphism on a set of generators; the morphism is automatically
prolonged by `\mathcal S`-linearity.
Actually, `f` was just the multiplication by `X^3` on `M`. We can
then redefine it simply as follows:

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
build the sequence (where `\S` is the corresponding Ore
polynomial ring)::

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
method :meth:`is_zero`.

Of course, here, one can, more simply, use the method
:meth:`is_injective`::

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
`g \circ f` vanishes:

    sage: h = g*f
    sage: h
    Ore module morphism:
      From: Ore module <u0, u1> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module <w0, w1, w2> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: h.is_zero()
    True

Let us now consider another morphism `f` and build the
canonical isomorphism ``\text{coim }f \to \text{im }f`
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

Now we want to factor `g` by the coimage. We proceed as follows::

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

    sage: H = C.Hom(I)
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
`\mathcal S / \mathcal S P` is the reduced norm of P.

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
from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.categories.map import Map
from sage.categories.morphism import Morphism
from sage.modules.ore_module import OreModule, OreSubmodule, OreQuotientModule

class OreModuleMorphism(Morphism):
    def __init__(self, parent, im_gens, check=True):
        Morphism.__init__(self, parent)
        domain = parent.domain()
        codomain = parent.codomain()
        MS = parent.matrix_space()
        if isinstance(im_gens, Matrix):
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
        if check:
            for x in parent.domain().basis():
                if self._call_(x.image()) != self._call_(x).image():
                    raise ValueError("does not define a morphism of Ore modules")

    def _repr_type(self):
        return "Ore module"

    def _latex_(self):
        s = "\\begin{array}{l}\n"
        s += "\\text{\\texttt{%s morphism:}} \\\\\n" % self._repr_type()
        s += "\\text{\\texttt{{ }{ }From:}}\\hspace{1ex} %s \\\\\n" % latex(self.domain())
        s += "\\text{\\texttt{{ }{ }To:}}\\hspace{3ex} %s \n" % latex(self.codomain())
        s += "\\end{array}"
        return s

    def matrix(self):
        return self._matrix.__copy__()

    def is_zero(self):
        return self._matrix.is_zero()

    def is_identity(self):
        return self.domain() is self.codmain() and self._matrix.is_one()

    def determinant(self):
        if self.domain() is not self.codomain():
            raise ValueError("determinants are only defined for endomorphisms")
        return self._matrix.determinant()

    det = determinant

    def characteristic_polynomial(self, var='x'):
        if self.domain() is not self.codomain():
            raise ValueError("characteristic polynomials are only defined for endomorphisms")
        return self._matrix.charpoly(var)

    charpoly = characteristic_polynomial

    def is_injective(self):
        return self._matrix.rank() == self.domain().rank()

    def is_surjective(self):
        return self._matrix.rank() == self.codomain().rank()

    def is_bijective(self):
        return self.is_injective() and self.is_surjective()

    def is_isomorphism(self):
        return self.is_bijective()

    def _call_(self, x):
        return self.codomain()(x * self._matrix)

    def _composition_(self, other, homset):
        if not isinstance(other, OreModuleMorphism):
            raise ValueError(str(other))
        return homset(other._matrix * self._matrix)

    def kernel(self, names=None):
        ker = self._matrix.left_kernel_matrix()
        return OreSubmodule(self.domain(), ker, names)

    def image(self, names=None):
        return OreSubmodule(self.codomain(), self._matrix, names)

    def cokernel(self, names=None):
        return OreQuotientModule(self.codomain(), self._matrix, names)

    def coimage(self, names=None):
        ker = self._matrix.left_kernel_matrix()
        return OreQuotientModule(self.domain(), ker, names)

class OreModuleRetraction(Map):
    def _call_(self, y):
        X = self.codomain()
        xs = X._basis.solve_left(y)
        return X(xs)

class OreModuleSection(Map):
    def _call_(self, y):
        X = self.codomain()
        Y = self.domain()
        indices = Y._indices
        zero = X.base_ring().zero()
        xs = X.rank() * [zero]
        for i in range(Y.rank()):
            xs[indices[i]] = y[i]
        return X(xs)
