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

from sage.categories.fields import Fields
from sage.categories.ore_modules import OreModules

from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix

from sage.modules.free_module import FreeModule_ambient
from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.modules.ore_module_element import OreModule_element


class OreModule(FreeModule_ambient):
    Element = OreModule_element

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

    def _repr_element(self, x):
        return FreeModuleElement_generic_dense._repr_(x)

    def _Hom_(self, codomain, category):
        from sage.modules.ore_module_homspace import OreModule_homspace
        return OreModule_homspace(self, codomain)

    def hom(self, f, codomain=None):
        from sage.modules.ore_module_morphism import OreModule_morphism
        if codomain is None:
            codomain = self
        H = self.Hom(codomain)
        if isinstance(f, OreModule_morphism):
            if f.domain() is not self:
                f = self._hom_change_domain(f)
            if f.codomain() is not codomain:
                f = codomain._hom_change_codomain(f)
            return f
        return H(f)

    def pseudohom(self):
        return self._pseudohom

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
            B.append(self.element_class(self, coeffs))
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
        return self.element_class(self, coeffs)

    def module(self):
        return self.base_ring() ** self.rank()

    def _span(self, gens):
        if not isinstance(gens, (list, tuple)):
            raise ValueError("not a list of generators")
        f = self._pseudohom
        rank = self.rank()
        d = len(gens)
        M = matrix(self.base_ring(), max(d, 2*rank), rank)
        for i in range(d):
            v = gens[i].list()
            for j in range(rank):
                M[i,j] = v[j]
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
        if isinstance(sub, OreSubmodule):
            if sub._M is not self:
                raise ValueError("not a submodule")
            gens = sub._basis
        else:
            gens = self._span(sub)
        return self._quotientModule_class(self, gens, names=names)

    quo = quotient

    def __eq__(self, other):
        return self is other

    def __hash__(self):
        return id(self)


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
        OreModule.__init__(self, matrix(base, rows), M.ore_ring(), names)
        self._inject = coerce = self.hom(basis, codomain=M)
        self._basis = basis
        M.register_coercion(coerce)

    def _repr_element(self, x):
        return self._M(x)._repr_()

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
        OreModule.__init__(self, matrix(base, rows), M.ore_ring(), names)
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

    def dividend(self):
        return self._M

    #@cached_method
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
