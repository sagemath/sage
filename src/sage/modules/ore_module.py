r"""
Ore modules.

Let `R` be a commutative ring, `\theta : K \to K` by a ring
endomorphism and `\partial : K \to K` be a `\theta`-derivation,
that is an additive map satisfying the following axiom

.. MATH::

    \partial(x y) = \theta(x) \partial(y) + \partial(x) y

A Ore module over `(R, \theta, \partial)` is a `R`-module `M`
equipped with a additive `f : M \to M` such that

.. MATH::

    f(a x) = \theta(a) f(x) + \partial(a) x

Such a map `f` is called a pseudomorphism.

Equivalently, a Ore module is a module over the (noncommutative)
Ore polynomial ring `\mathcal S = R[X; \theta, \partial]`.

.. RUBRIC:: Defining Ore modules

SageMath provides support for creating and manipulating Ore
modules that are finite free over the base ring `R`.

To start with, the method :meth:`quotient_module` creates the
quotient `\mathcal S/ \mathcal S P`, endowed with its structure
of `\mathcal S`-module, that is its structure of Ore module::

    sage: K.<z> = GF(5^3)
    sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
    sage: M = S.quotient_module(X^2 + z)
    sage: M
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

Classical methods are available and we can work with elements in
`M` as we do usually for vectors in finite free modules::

    sage: M.basis()
    [(1, 0), (0, 1)]

    sage: v = M((z, z^2)); v
    (z, z^2)
    sage: z*v
    (z^2, 2*z + 2)

The Ore action (or equivalently the structure of `mathcal A`-module)
is also easily accessible)::

    sage: X*v
    (3*z^2 + 2*z, 2*z^2 + 4*z + 4)

The method :meth:`pseudohom` returns the map `f` defining the action
of `X`::

    sage: M.pseudohom()
    Free module pseudomorphism (twisted by z |--> z^5) defined by the matrix
    [  0   1]
    [4*z   0]
    Domain: Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    Codomain: Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

A useful feature is the possibility to give chosen names to the vectors
of the canonical basis. This is easily done as follows:

    sage: N.<u,v,w> = S.quotient_module(X^3 + z*X + 1)
    sage: N
    Ore module <u, v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: N.basis()
    [u, v, w]

Alternatively, one can pass in the argument ``names``; this could be
useful in particular when we want to name them `e_0, e_1, \ldots`::

    sage: A = S.quotient_module(X^11 + z, names='e')
    sage: A
    Ore module <e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: A.basis()
    [e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10]

Do not forget to use the method :meth:`inject_variables` to get the
`e_i` in your namespace::

    sage: e0
    Traceback (most recent call last):
    ...
    NameError: name 'e0' is not defined
    sage: A.inject_variables()
    Defining e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10
    sage: e0
    e0

.. RUBRIC:: Submodules and quotients

SageMath provides facilities for creating submodules and quotient
modules of Ore modules.
First of all, we define the Ore module `\mathcal S/\mathcal S P^2`
(for some Ore polynomials `P`), which is obviously not simple::

    sage: P = X^2 + z*X + 1
    sage: U = S.quotient_module(P^2, names='u')
    sage: U.inject_variables()
    Defining u0, u1, u2, u3

We now build the submodule `\mathcal S Q / \mathcal S PQ` using
the method :meth:`span`::

    sage: V = U.span(P*u0)
    sage: V
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: V.basis()
    [u0 + (z^2 + 2*z + 2)*u2 + 4*z*u3,
     u1 + (2*z^2 + 4*z + 4)*u2 + u3]

We underline that the span is really the `\mathcal S`-span and
not the `R`-span (as otherwise, it will not be a Ore module).

As before, one can use the attributes ``names`` to give explicit
names to the basis vectors::

    sage: V = U.span(P*u0, names='v')
    sage: V
    Ore module <v0, v1> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: V.inject_variables()
    Defining v0, v1
    sage: v0
    v0
    sage: U(v0)
    u0 + (z^2 + 2*z + 2)*u2 + 4*z*u3

A coercion map from `V` to `U` is automatically created and set up.
Hence, we can safely combine vectors in `V` and vectors in `U` in a
single expression::

    sage: v0 - u0
    (z^2 + 2*z + 2)*u2 + 4*z*u3

We can create the quotient `U/V` using a similar syntax::

    sage: W = U.quo(P*u0)
    sage: W
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: W.basis()
    [u2, u3]

We see that SageMath reuses by default the names of the representatives
to denote the vectors in the quotient `U/V`. This behaviour can be overrided
by providing explicit names using the attributes ``names``.

Shortcuts for creating quotients are also available::

    sage: U / (P*u0)
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: U/V
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

.. RUBRIC:: Morphisms of Ore modules

By definition, a morphism of Ore module is a `R`-linear morphism
commuting with the Ore action, or equivalenty a `\mathcal S`-linear
map.

There are several ways for creating Ore modules morphisms in SageMath.
First of all, one can use the method :meth:`hom`, passing in to it the
matrix (in the canonical bases) of the morphism we want to build::

    sage: mat = matrix(2, 2, [3*z^2 + z + 2, 2*z,
    ....:                     3*z^2 + z + 1, 4])
    sage: f = V.hom(mat, codomain=W)
    sage: f
    Ore module morphism:
      From: Ore module <v0, v1> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

Clearly, this method is not optimal: typing all the entries of the
defining matrix is long and is a potential source of errors.
Instead, one can use a dictionary encoding the values taken by the
morphism on a set of generators; the morphism is automatically
prolonged by `\mathcal S`-linearity::

    sage: g = V.hom({P*u0: W(u0)})
    sage: g
    Ore module morphism:
      From: Ore module <v0, v1> over Finite Field in z of size 5^3 twisted by z |--> z^5
      To:   Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

One can then recover the matrix of `g` by using the method :meth:`matrix`::

    sage: g.matrix()
    [3*z^2 + z + 2           2*z]
    [3*z^2 + z + 1             4]

The method :meth:`multiplication_map` can also be used to define
scalar endomorphisms::

    sage: h = U.multiplication_map(2)
    sage: h
    Ore module endomorphism of Ore module <u0, u1, u2, u3> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: h.matrix()
    [2 0 0 0]
    [0 2 0 0]
    [0 0 2 0]
    [0 0 0 2]

Be careful that scalar multiplications do not always properly
define a morphism of Ore modules::

    sage: U.multiplication_map(z)
    Traceback (most recent call last):
    ...
    ValueError: does not define a morphism of Ore modules

The method :meth:`multiplication_map` also accepts a Ore
polynomial::

    sage: h = U.multiplication_map(X^3)
    sage: h
    Ore module endomorphism of Ore module <u0, u1, u2, u3> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: h.matrix()
    [              0               0               0               1]
    [              4             3*z   z^2 + 2*z + 4 2*z^2 + 4*z + 4]
    [      2*z^2 + 4           z + 2         2*z + 3               3]
    [              2             3*z               3               3]

For endomorphisms, one can compute classical invariants as
determinants and characteristic polynomials::

    sage: h.det()
    1
    sage: h.charpoly()
    x^4 + 4*x^3 + x^2 + 4*x + 1

One can check that the latter is the same than the reduced
norm of `P^2`::

    sage: (P^2).reduced_norm()
    z^4 + 4*z^3 + z^2 + 4*z + 1

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

import operator
from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.categories.action import Action
from sage.categories.fields import Fields
from sage.categories.ore_modules import OreModules

from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.modules.free_module import FreeModule_ambient
from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.modules.ore_module_element import OreModuleElement

# Action of Ore polynomials on Ore modules
##########################################

class OreAction(Action):
    def _act_(self, P, x):
        M = x.parent()
        ans = P[0]*x
        y = x
        for i in range(1, P.degree() + 1):
            y = y.image()
            ans += y._rmul_(P[i])
        return ans

class ScalarAction(Action):
    def _act_(self, a, x):
        return x._rmul_(a)

# Generic class for Ore modules
###############################

class OreModule(FreeModule_ambient):
    # TODO: ensure uniqueness of parents
    Element = OreModuleElement

    def __init__(self, mat, twist, names=None, category=None):
        r"""
        Initialize this Ore module.

        INPUT:

        - ``mat`` -- the matrix defining the action of the Ore variable

        - ``twist`` -- the twisting morphism/derivation

        - ``names`` (default: ``None``) -- a string of a list of stings,
          the names of the vector of the canonical basis; if ``None``,
          elements are represented as vectors in `K^d`

        - ``category`` (default: ``None``) -- the category of this
          Ore module

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)  # indirect doctest
            sage: type(M)
            <class 'sage.modules.ore_module.OreModule_with_category'>

        """
        base = mat.base_ring()
        if category is None:
            category = OreModules(base, twist)
        if base not in Fields():
            raise NotImplementedError("Ore modules are only implemented over fields")
        rank = mat.nrows()
        if mat.ncols() != rank:
            raise ValueError("matrix must be square")
        FreeModule_ambient.__init__(self, base, rank, category=category)
        self.register_action(ScalarAction(base, self, True, operator.mul))
        self._ore = category._ore
        self._pseudohom = FreeModule_ambient.pseudohom(self, mat, self._ore, codomain=self)
        if names is None:
            pass
        elif isinstance(names, (list, tuple)):
            if rank != len(names):
                raise ValueError
            names = [str(name) for name in names]
        elif isinstance(names, str):
            names = [ names + str(i) for i in range(rank) ]
        else:
            raise ValueError
        self._names = names
        if names is not None:
            self._latex_names = [latex_variable_name(name) for name in names]
        self._submodule_class = OreSubmodule
        self._quotientModule_class = OreQuotientModule

    def _repr_(self):
        s = "Ore module "
        if self._names is None:
            s += "of rank %s " % self.rank()
        else:
            s += "<" + ", ".join(self._names) + "> "
        s += "over %s %s" % (self.base_ring(), self._ore._repr_twist())
        return s

    def _latex_(self):
        if self._names is None:
            s = "\\texttt{Ore module of rank } %s" % self.rank()
            s += "\\texttt{ over } %s" % latex(self.base_ring())
            twist = self._ore._latex_twist()
            if twist == "":
                s += "\\texttt{ untwisted}"
            else:
                s += "\\texttt{ twisted by }" + twist
        else:
            s = "\\left<" + ", ".join(self._latex_names) + "\\right>"
            s += "_{%s" % latex(self.base_ring())
            twist = self._ore._latex_twist()
            if twist != "":
                s += "," + twist
            s += "}"
        return s

    def _repr_element(self, x):
        return FreeModuleElement_generic_dense._repr_(x)

    def _latex_element(self, x):
        return FreeModuleElement_generic_dense._latex_(x)

    def pseudohom(self):
        return self._pseudohom

    def ore_ring(self, names='x', action=True):
        S = self.category().ore_ring(names)
        if action:
            self._unset_coercions_used()
            self.register_action(OreAction(S, self, True, operator.mul))
        return S

    def matrix(self):
        return self._pseudohom.matrix()

    def basis(self):
        rank = self.rank()
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        coeffs = [zero] * rank
        B = [ ]
        for i in range(rank):
            coeffs[i] = one
            B.append(self(coeffs))
            coeffs[i] = zero
        return B

    def gens(self):
        return self.basis()

    def gen(self, i):
        rank = self.rank()
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        coeffs = [zero] * rank
        coeffs[i] = one
        return self(coeffs)

    def random_element(self, *args, **kwds):
        K = self.base_ring()
        r = self.rank()
        vs = [K.random_element(*args, **kwds) for _ in range(r)]
        return self(vs)

    def module(self):
        return self.base_ring() ** self.rank()

    def _Hom_(self, codomain, category):
        from sage.modules.ore_module_homspace import OreModule_homspace
        return OreModule_homspace(self, codomain)

    def hom(self, im_gens, codomain=None):
        from sage.modules.ore_module_morphism import OreModuleMorphism
        if codomain is None:
            if isinstance(im_gens, Matrix):
                codomain = self
            elif isinstance(im_gens, OreModuleMorphism):
                codomain = im_gens.codomain()
            elif isinstance(im_gens, (list, tuple)):
                codomain = Sequence(im_gens).universe()
            elif isinstance(im_gens, dict):
                codomain = Sequence(im_gens.values()).universe()
            else:
                raise ValueError("im_gens must be a list, a tuple, a dictionary, a matrix or a Ore module morphism")
        H = self.Hom(codomain)
        if isinstance(im_gens, Matrix):
            return H(im_gens)
        elif isinstance(im_gens, OreModuleMorphism):
            f = im_gens
            if f.domain() is not self:
                f = self._hom_change_domain(f)
            if f.codomain() is not codomain:
                f = codomain._hom_change_codomain(f)
            return f
        elif isinstance(im_gens, (list, tuple)):
            if len(im_gens) != self.rank():
                raise ValueError("wrong number of generators")
            M = matrix([codomain(v).list() for v in im_gens])
            return H(M)
        elif isinstance(im_gens, dict):
            zero = self.base_ring().zero()
            dimd = self.rank()
            dimc = codomain.rank()
            d = dimc + dimd
            vs = [self(x).list() + codomain(y).list() for x, y in im_gens.items()]
            if len(vs) < 2*d:
                vs += (2*d - len(vs)) * [d * [zero]]
            M = matrix(vs)
            M.echelonize()
            oldr = 0
            r = M.rank()
            iter = 1
            fd = self._pseudohom
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
            M = M.submatrix(0, dimd, dimd, dimc)
            return H(M)
        else:
            raise ValueError("im_gens must be a list, a tuple, a dictionary, a matrix or a Ore module morphism")

    def multiplication_map(self, P):
        if isinstance(P, OrePolynomial):
            S = P.parent()
            ore = self.category()._ore
            if S._morphism != ore._morphism or S._derivation != ore._derivation:
                raise ValueError("twist does not match")
            action = OreAction(S, self, True, operator.mul)
            M = matrix([action._act_(P, x).list() for x in self.basis()])
        else:
            P = self.base_ring()(P)
            r = self.rank()
            M = matrix(r, r, P)
        H = self.Hom(self)
        return H(M)

    def identity_morphism(self):
        H = self.Hom(self)
        M = identity_matrix(self.base_ring(), self.rank())
        return H(M)

    def _span(self, gens):
        base = self.base_ring()
        rank = self.rank()
        f = self._pseudohom
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        rows = []
        for gen in gens:
            if isinstance(gen, OreModule):
                incl = self.coerce_map_from(gen)
                if incl is None:
                    raise ValueError("not canonically a submodule")
                rows += incl._matrix.rows()
            elif isinstance(gen, OreModuleElement):
                rows.append(self(gen).list())
        if len(rows) < 2*rank:
            zero = rank * [base.zero()]
            rows += (2*rank - len(rows)) * [rank*[0]]
        M = matrix(base, rows)
        M.echelonize()
        oldr = 0
        r = M.rank()
        iter = 1
        while r > oldr:
            for i in range(r):
                v = M.row(i)
                for _ in range(iter):
                    v = f(v)
                v = v.list()
                for j in range(rank):
                    M[i+r,j] = v[j]
            M.echelonize()
            oldr = r
            r = M.rank()
            iter *= 2
        return M.matrix_from_rows(range(r))

    def span(self, gens, names=None):
        gens = self._span(gens)
        return self._submodule_class(self, gens, names=names)

    def quotient(self, sub, names=None, check=True):
        gens = self._span(sub)
        return self._quotientModule_class(self, gens, names=names)

    quo = quotient

    def __eq__(self, other):
        return self is other

    def __hash__(self):
        return id(self)

# Submodules
############

class OreSubmodule(OreModule):
    def __init__(self, M, gens, names):
        from sage.modules.ore_module_morphism import OreModuleRetraction
        self._M = M
        base = M.base_ring()
        if isinstance(gens, Matrix):
            basis = gens
        else:
            basis = matrix(base, gens)
        basis = basis.echelon_form()
        rank = basis.rank()
        if basis.nrows() != rank:
            basis = basis.matrix_from_rows(range(rank))
        rows = [basis.solve_left(M(x).image()) for x in basis.rows()]
        OreModule.__init__(self, matrix(base, rows), M.ore_ring(action=False), names)
        self._inject = coerce = self.hom(basis, codomain=M)
        self._basis = basis
        M.register_coercion(coerce)
        self.register_conversion(OreModuleRetraction(M, self))

    def _repr_element(self, x):
        return self._M(x)._repr_()

    def _latex_element(self, x):
        return self._M(x)._latex_()

    def ambient(self):
        return self._M

    def injection_morphism(self):
        return self._inject

    def _hom_change_domain(self, f):
        if f.domain() is not self._M:
            raise ValueError
        return f * self._inject

    def _hom_change_codomain(self, f):
        if f.codomain() is not self._M:
            raise ValueError
        rows = []
        basis = self._basis
        rows = [basis.solve_left(y) for y in f._matrix.rows()]
        return f.domain().hom(rows, codomain=self)

# Quotients
###########

class OreQuotientModule(OreModule):
    def __init__(self, M, gens, names):
        self._M = M
        d = M.rank()
        base = M.base_ring()
        if isinstance(gens, Matrix):
            basis = gens
        else:
            basis = matrix(base, gens)
        self._ker = basis = basis.echelon_form()
        pivots = basis.pivots()
        r = basis.rank()
        coerce = matrix(base, d, d-r)
        indices = []
        i = 0
        for j in range(d):
            if i < r and pivots[i] == j:
                i += 1
            else:
                indices.append(j)
                coerce[j,j-i] = base.one()
        for i in range(r):
            for j in range(d-r):
                coerce[pivots[i],j] = -basis[i,indices[j]]
        rows = [M.gen(i).image() * coerce for i in indices]
        OreModule.__init__(self, matrix(base, rows), M.ore_ring(action=False), names)
        self._indices = indices
        self._project = coerce = M.hom(coerce, codomain=self)
        self.register_coercion(coerce)

    def _repr_element(self, x):
        M = self._M
        indices = self._indices
        base = self.base_ring()
        coords = M.rank() * [base.zero()]
        for i in range(self.rank()):
            coords[indices[i]] = x[i]
        return M(coords)._repr_()

    def _latex_element(self, x):
        M = self._M
        indices = self._indices
        base = self.base_ring()
        coords = M.rank() * [base.zero()]
        for i in range(self.rank()):
            coords[indices[i]] = x[i]
        return M(coords)._latex_()

    def dividend(self):
        return self._M

    @cached_method
    def divisor(self, names=None):
        return self._submodule_class(self._M, self._ker, names=names)

    def projection_morphism(self):
        return self._project

    def _hom_change_domain(self, f):
        if f.domain() is not self._M:
            raise ValueError
        Z = self._ker * f._matrix
        if not Z.is_zero():
            raise ValueError
        mat = f._matrix.matrix_from_rows(self._indices)
        return self.hom(mat, codomain=f.codomain())

    def _hom_change_codomain(self, f):
        if f.codomain() is not self._M:
            raise ValueError
        return self._project * f
