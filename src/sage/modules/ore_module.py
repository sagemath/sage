r"""
Ore modules

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

    def __init__(self, f, twist=None, names=None, category=None):
        base = f.base_ring()
        if category is None:
            category = OreModules(base, twist)
        if base not in Fields():
            raise NotImplementedError("Ore modules are only implemented over fields")
        rank = f.nrows()
        if f.ncols() != rank:
            raise ValueError("matrix must be square")
        FreeModule_ambient.__init__(self, base, rank, category=category)
        self.register_action(ScalarAction(base, self, True, operator.mul))
        self._ore = category._ore
        self._pseudohom = FreeModule_ambient.pseudohom(self, f, self._ore, codomain=self)
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
        from sage.modules.ore_module_morphism import OreModule_morphism
        if codomain is None:
            if isinstance(im_gens, Matrix):
                codomain = self
            elif isinstance(im_gens, OreModule_morphism):
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
        elif isinstance(im_gens, OreModule_morphism):
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
