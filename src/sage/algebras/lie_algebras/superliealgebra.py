from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras import Algebras
from sage.misc.cachefunc import cached_method
from sage.rings.ring import Algebra
from sage.misc.functional import is_odd, is_even
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.rings.integer_ring import ZZ
from sage.structure.indexed_generators import (IndexedGenerators, standardize_names_index_set)
from sage.sets.family import Family
from sage.modules.free_module import FreeModule
from sage.rings.quotient_ring_element import QuotientRingElement
import unittest
from sage.modules.free_module_morphism import FreeModuleMorphism
from sage.modules.with_basis.morphism import ModuleMorphism
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_attribute import lazy_attribute

class SuperLieAlgebra(CombinatorialFreeModule, UniqueRepresentation):

    @staticmethod
    def __classcall_private__(cls, base, s_coeff, names=None, index_set=None, degrees=None, category=None, **kwargs):

        if names is None:
            if degrees is None:
                raise ValueError("You must specify names or degrees")
            else:
                n = len(degrees)
            names = tuple('x{}'.format(i) for i in range(n))

        elif isinstance(names, str):
            names = tuple(names.split(','))
            n = len(names)

        else:
            n = len(names)
            names = tuple(names)

        if degrees is None:
            degrees = tuple([1] * n)

        else:
            degrees = tuple(degrees)

        names, index_set = standardize_names_index_set(names, index_set)

        if names is not None and names != tuple(index_set):
            d = {x: index_set[i] for i,x in enumerate(names)}
            get_pairs = lambda X: X.items() if isinstance(X, dict) else X
            try:
                s_coeff = {(d[k[0]], d[k[1]]): [(d[x], y) for x,y in get_pairs(s_coeff[k])]
                           for k in s_coeff}
            except (KeyError, ValueError):
                pass

        s_coeff = SuperLieAlgebra._standardize_s_coeff(names, s_coeff, index_set, degrees)

        return super().__classcall__(cls, base=base, names=names, s_coeff=s_coeff, index_set=index_set,
                                     degrees=degrees, category=category, **kwargs)

    #Standardizing the bracket dictionary:

    def _standardize_s_coeff(names, s_coeff, index_set, degrees):

        index_to_pos = {k: i for i,k in enumerate(index_set)}

        sc = {}

        for k in s_coeff.keys():
            v = s_coeff[k]
            if isinstance(v, dict):
                v = v.items()

            if index_to_pos[k[0]] > index_to_pos[k[1]]:
                li=names.index(k[0])
                lj=names.index(k[1])
                i = degrees[li]
                j=degrees[lj]
                key = (k[1], k[0])
                vals = tuple((g, (-1)**(i*j+1)*val) for g, val in v if val != 0)
            else:
                key = tuple(k)
                vals = tuple((g, val) for g, val in v if val != 0)

            if key in sc.keys() and sorted(sc[key]) != sorted(vals):
                raise ValueError("two distinct values given for one and the same bracket")

            if vals:
                sc[key] = vals
        return Family(sc)

    def __init__(self, base, s_coeff, names, index_set, degrees, category=None):


        self._index_to_pos = {k: i for i,k in enumerate(index_set)}

        self._names = names
        self._ngens = len(self._names)
        self._degrees = degrees
        self.name_degree_map = dict(zip(names, degrees))
        self._M = FreeModule(base, len(index_set))
        self.s_coeff = s_coeff
        self.vector_presentation = dict(zip(names, self._M.basis()))

        def to_vector(tuples):
            vec = [base.zero()]*len(index_set)
            for k,c in tuples:
                vec[self._index_to_pos[k]] = c
            vec = self._M(vec)
            vec.set_immutable()
            return vec

        self._s_coeff = {(self._index_to_pos[k[0]], self._index_to_pos[k[1]]): to_vector(s_coeff[k]) for k in s_coeff.keys()}

        base_cat = LieAlgebras(base).WithBasis().FiniteDimensional().or_subcategory(category)
        category = base_cat.or_subcategory(category, join=True)

        CombinatorialFreeModule.__init__(self, base, index_set)

    def some_elements(self):

        return list(self.basis())

    def base_module(self):

        return VectorSpace(QQ, len(list(self.basis())))

    @cached_method
    def zero(self):

        return self.element_class(self, {})

    def bracket(self, lt, rt):

        return lt.bracket(rt)


    #Morphisms of graded Lie algebras:

    def gradedmorphism(self, on_generators, domain, codomain, argument):

        from sage.categories.lie_algebras import LieAlgebras
        from itertools import combinations

        m = domain.base_module()
        cm = codomain.base_module()

        spanning_set = [X.to_vector() for X in list(on_generators)]
        im_gens = [Y.to_vector() for Y in list(on_generators.values())]

        def solve_linear_system(A, b):
            R = QQ
            A_inv = A.solve_left(matrix.identity(A.ncols()))
            M = A * A_inv
            for Mi, bk in zip(M.rows(), b):
                test_bk = sum((R(Mij) * bj for Mij,bj in zip(Mi,b)), cm.zero())
                if test_bk != bk:
                    raise ValueError("contradictory linear system")

            return [sum((R(Aij) * bk for Aij,bk in zip(Ai,b)), cm.zero())
                    for Ai in A_inv.rows()]

        bracketlength = 1
        n = 0
        while True:
            sm = m.submodule(spanning_set)
            A = matrix(sm.base_ring(), [sm.coordinate_vector(X) for X in spanning_set])
            try:
                im_gens = solve_linear_system(A, im_gens)
            except ValueError:
                raise ValueError("this does not define a graded Lie algebra morphism")

            spanning_set = list(sm.basis())
            if n == len(spanning_set):

                break

            bracketlength += 1
            n = len(spanning_set)
            for i,j in combinations(range(n), 2):

                Z = list(on_generators.keys())[i].bracket(list(on_generators.keys())[j])
                imZ = list(on_generators.values())[i].bracket(list(on_generators.values())[j])
                spanning_set.append(Z.to_vector())
                im_gens.append(imZ.to_vector())

        A = matrix(m.base_ring(), spanning_set)
        im_gens = solve_linear_system(A, im_gens)

        bh = lambda t: t
        return cm.sum(bh(c) * im_gens[i] for i, c in (argument.to_vector()).items())


    class Element(SuperLieAlgebra.Element):
         
        #Returns the degree of the highest-degree monomial in any expression:
        def degree(self):

            l = [m for m in self]
            monomials = [x[0] for x in l]
            dictionary = self.parent().name_degree_map
            degrees=self.parent()._degrees

            values = [dictionary[key] for key in monomials if key in dictionary]
            if values:
                return max(values)
            else:
                return 0
            
        #Returns the degree of a homogeneous expression:
        def homdegree(self):

            l = [m for m in self]
            monomials = [x[0] for x in l]
            dictionary = self.parent().name_degree_map
            degrees=self.parent()._degrees

            values = [dictionary[key] for key in monomials if key in dictionary]
            if values:
                if len(values)!=1:
                    return "Not homogeneous"
                else:
                    return values[0]
            else:
                return 0

        def bracket(self, rt):
            names=self.parent()._names
            degrees=self.parent()._degrees
            dictionary = self.parent().vector_presentation
            basis=self.parent().basis()

            structurecoefficients = self.parent()._s_coeff
            d = self.parent().dimension()
            ret = [0]*d

            if self == 0:
                return 0
            if rt == 0:
                return 0

            l1 = [m for m in self]
            monomials1 = [x[0] for x in l1]
            values1 = [dictionary[key] for key in monomials1 if key in dictionary]

            vector1 = 0
            for i in range(len(l1)):
                vector1 += self.coefficients()[i]*values1[i]

            l2 = [m for m in rt]
            monomials2 = [x[0] for x in l2]
            values2 = [dictionary[key] for key in monomials2 if key in dictionary]

            vector2 = 0
            for i in range(len(l2)):
                vector2 += rt.coefficients()[i]*values2[i]


            for i1 in range(d):
                c1 = vector1[i1]
                if not c1:
                    continue
                for i2 in range(d):
                    c2 = vector2[i2]
                    if not c2:
                        continue
                    prod_c1_c2 = c1 * c2
                    if (i1, i2) in structurecoefficients:
                        v = structurecoefficients[i1, i2]
                        for i3 in range(d):
                            ret[i3] += prod_c1_c2 * v[i3]
                    elif (i2, i1) in structurecoefficients:
                        v = structurecoefficients[i2, i1]
                        for i3 in range(d):
                            i=degrees[i2]
                            j=degrees[i1]
                            ret[i3] += (-1)**(i*j+1)*prod_c1_c2 * v[i3]

            output = self.parent()._M(ret)
            final = 0
            for i in range(d):
                final += basis.values()[i]*output[i]
            return final

    def _bracket_(self, lhs, rhs):
         return self(lhs).bracket(self(rhs))

    def _homdegree_(self, lhs):
         return self(lhs).homdegree()

    #The following tests ensure that the bracket satisfies the relevant properties.

    def _test_antisymmetry(self, **options):
        tester = self._tester(**options)
        elts = tester.some_elements()
        zero = self.zero()
        for x in elts:
            tester.assertEqual(self._bracket_(x, y), self._bracket_(y,x)*(-1)**(x.degree()*y.degree()+1))

    def _test_distributivity(self, **options):
            tester = self._tester(**options)
            S = tester.some_elements()
            from sage.misc.misc import some_tuples
            for x,y,z in some_tuples(S, 3, tester._max_runs):

                tester.assertEqual(self._bracket_(x, (y + z)),
                                   self._bracket_(x, y) + self._bracket_(x, z))

                tester.assertEqual(self._bracket_((x + y), z),
                                   self._bracket_(x, z) + self._bracket_(y, z))


    def _test_degree(self, **options):
        tester = self._tester(**options)
        elts = tester.some_elements()
        for x in elts:
                for y in elts:
                    if self._bracket_(x,y).homdegree() != 0:
                        tester.assertEqual(self._homdegree_(x)+self._homdegree_(y), (self._bracket_(x,y)).homdegree())



    def _test_jacobi_identity(self, **options):
            tester = self._tester(**options)
            elts = tester.some_elements()
            zero = self.zero()
            for x in elts:
                for y in elts:
                    if x == y:
                        continue
                    for z in elts:
                        a = (-1)**(x.degree()*z.degree())*self._bracket_(x, self._bracket_(y, z))
                        b = (-1)**(y.degree()*x.degree())*self._bracket_(y, self._bracket_(z, x))
                        c =  (-1)**(z.degree()*y.degree())*self._bracket_(z, self._bracket_(x, y))
                        tester.assertEqual(a+b+c, zero)


#Example:
#d = {('x','y'):{'z':1}}
#L.<x,y,z>=SuperLieAlgebra(QQ, s_coeff=d, degrees=(1,1,2))
#L.gradedmorphism({x:x, y:y, z:z}, L, L, x.bracket(y))
