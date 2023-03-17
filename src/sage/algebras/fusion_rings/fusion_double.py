"""
The Fusion Ring of the Drinfeld Double of a Finite Group
"""
# ****************************************************************************
#  Copyright (C) 2023 Wenqi Li
#                     Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.all import Algebras, AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.misc.misc import inject_variable
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ideal import Ideal

class FusionDouble(CombinatorialFreeModule):
    r"""
    This constructs the Fusion Ring of the modular
    tensor category of modules over the Drinfeld
    Double of a finite group. Usage is similar
    to :class:`FusionRing`::.

    INPUT:

    - ``G`` -- a finite group
    - ``prefix`` (optional: defaults to 's') a prefix for the names of simple objects.
    - ``inject_varables`` (optional): set TRUE to create variables for the simple objects.

    REFERENCES:

    - [BaKi2001]_ Chapter 3
    - [Mas1995]_

    EXAMPLES ::

        sage: G = DihedralGroup(5)
        sage: H = FusionDouble(G, inject_variables=True)
        sage: H.basis()
        Finite family {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8, 9: s9, 10: s10, 11: s11, 12: s12, 13: s13, 14: s14, 15: s15}
        sage: for x in H.basis():
        ....:     print ("%s : %s"%(x,x^2))
        ....:
        s0 : s0
        s1 : s0
        s2 : s0 + s1 + s3
        s3 : s0 + s1 + s2
        s4 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s5 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s6 : s0 + s1 + s11
        s7 : s0 + s1 + s13
        s8 : s0 + s1 + s15
        s9 : s0 + s1 + s12
        s10 : s0 + s1 + s14
        s11 : s0 + s1 + s6
        s12 : s0 + s1 + s9
        s13 : s0 + s1 + s7
        s14 : s0 + s1 + s10
        s15 : s0 + s1 + s8
        sage: s4*s5
        s1 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        sage: s4.ribbon()
        1
        sage: s5.ribbon()
        -1
        sage: s8.ribbon()
        zeta5^3
    """
    def __init__(self, G, prefix="s",inject_variables=False):
        self._G = G
        self._prefix = prefix
        self._names = {}
        self._elt = {}
        self._chi = {}
        count = 0
        for g in G.conjugacy_classes_representatives():
            for chi in G.centralizer(g).irreducible_characters():
                self._names[count] = "%s%s"%(prefix, count)
                self._elt[count] = g
                self._chi[count] = chi
                count += 1
        self._rank = count
        self._cyclotomic_order = G.exponent()
        cat = AlgebrasWithBasis(ZZ).Subobjects()
        CombinatorialFreeModule.__init__(self, ZZ, [k for k in self._names], category=cat)
        if inject_variables:
            self.inject_variables()

    def _repr_(self):
        return "The Fusion Ring of the Drinfeld Double of %s"%self._G

    def __call__(self, *args):
        if len(args) > 1:
            args = (args,)
        return super(GAlg, self).__call__(*args)

    def _element_constructor(self, k):
        return self.monomial(k)

    def inject_variables(self):
        for i in range(self._rank):
            inject_variable(self._names[i],self.monomial(i))

    @cached_method
    def s_ij(self, i, j):
        sum = 0
        G = self._G
        [i, j] = [x.support_of_term() for x in [i,j]]
        [a, chi_1] = [self._elt[i], self._chi[i]]
        [b, chi_2] = [self._elt[j], self._chi[j]]
        for g in G:
            if a*g*b*g.inverse() == g*b*g.inverse()*a:
                sum += chi_1(g*b*g.inverse()) * chi_2(g.inverse()*a*g)
        coef = 1 / (G.centralizer(a).order() * G.centralizer(b).order())
        return coef * sum

    def s_ijconj(self, i, j):
        return self.s_ij(i, j).conjugate()

    @cached_method
    def N_ijk(self, i, j, k):
        """
        The symmetric invariant of three simple objects,
        this returns the dimension of

        .. MATH::
           Hom(i \\otimes j\\otimes k, s_0)

        where `s_0` is the unit element (assuming prefix='s').
        Method of computation is through the Verlinde formula
        """
        sz = self.one()
        return ZZ(sum(self.s_ij(i, r) * self.s_ij(j, r) * self.s_ij(k, r)/self.s_ij(sz, r) for r in self.basis()))

    def Nk_ij(self, i, j, k):
        r"""
        Returns the fusion coefficient `N^k_{ij}`, computed using the Verlinde formula:

        .. MATH::
            N^k_{ij} = \sum_l \frac{s(i, \ell)\, s(j, \ell)\, \overline{s(k, \ell)}}{s(I, \ell)},

        """
        return self.N_ijk(i, j, self.dual(k))

    def char_Nk_ij(self, i, j, k):
        r"""
        Use character theoretic method to compute the fusion coefficient `N_{ij}^k`.
        Each simple object, for example `i` corresponds to a conjugacy class `\mathcal{C}_i`
        of the underlying group `G`, and an irreducible character `\chi_i` of the
        centralizer `C(g_i)` of a representative `g_i` of `\mathcal{C}_i`. In addition
        to a fixed representative `g_i` and `g_j` of the classes `\mathcal{C}_i`
        and `\mathcal{C}_j`, the formula will make use of variable elements `h_i` and
        `h_j` that are subject to the condition `h_ih_j=g_k`.

        .. MATH::

            \frac{|\mathcal{C}_i|}{|G|}\sum_{\substack{h_i\in\mathcal{C}_i\\ h_j\in\mathcal{C}_j\\ h_ih_j=g_k}}|C(h_i)\cap C(h_j)|\,
            \langle\chi_i^{(h_i)}\chi_j^{(h_j)},\chi_k\rangle_{C(h_i)\cap C(h_j)},

        where `\chi_i^{(h_i)}` is the character `\chi_i` of `C(g_i)` conjugated to a
        character of `C(h_i)`, when `h_i` is a conjugate of the fixed representative `g_i`.
        More exactly, there exists `r_i` such that `r_i g_i r_i^{-1}=h_i`, and then
        `\chi_i^{(h_i)}(x)=\chi_i(r_i^{-1}xr_i)`, and this definition does not
        depend on the choice of `r_i`.
        This formula is due to Christopher Goff when the centralizers are normal, and
        to Wenqi Li in the general case.
        """
        G = self._G
        I = G.conjugacy_class(i.g())
        J = G.conjugacy_class(j.g())
        IJ = Set(I_elem * J_elem for I_elem in I for J_elem in J)
        if k.g() not in IJ:
            return 0

        K = G.conjugacy_class(k.g())
        CI = G.centralizer(i.g())
        CJ = G.centralizer(j.g())
        CK = G.centralizer(k.g())

        c = K.cardinality() / G.order()
        summands = [(I_elem, J_elem) for I_elem in I for J_elem in J if I_elem * J_elem == k.g()]
        res = 0
        for p in summands:
            I_elem, J_elem = p
            for g in G:
                if g.inverse() * i.g() * g == I_elem:
                    i_twist = g
                if g.inverse() * j.g() * g == J_elem:
                    j_twist = g
            A = Set(i_twist.inverse() * zi * i_twist for zi in CI)
            B = Set(j_twist.inverse() * zj * j_twist for zj in CJ)
            inner_summands = A.intersection(B).intersection(Set(CK))
            for x in inner_summands:
                res += i.chi()(i_twist * x * i_twist.inverse()) * j.chi()(j_twist * x * j_twist.inverse()) * k.chi()(x).conjugate()
        return c * res

    def field(self):
        return CyclotomicField(4 * self._cyclotomic_order)

    def root_of_unity(self, r, base_coercion=True):
        r"""
        Return `e^{i\pi r}` as an element of ``self.field()`` if possible.

        INPUT:

        - ``r`` -- a rational number

        EXAMPLES::

            sage: A11 = FusionRing("A1", 1)
            sage: A11.field()
            Cyclotomic Field of order 24 and degree 8
            sage: for n in [1..7]:
            ....:     try:
            ....:         print(n, A11.root_of_unity(2/n))
            ....:     except ValueError as err:
            ....:         print(n, err)
            1 1
            2 -1
            3 zeta24^4 - 1
            4 zeta24^6
            5 not a root of unity in the field
            6 zeta24^4
            7 not a root of unity in the field
        """
        n = 2 * r * self._cyclotomic_order
        if n not in ZZ:
            raise ValueError("not a root of unity in the field")
        return  self.field().gen() ** n

    def r_matrix(self, i, j, k):

        r"""
        """
        if self.Nk_ij(i, j, k) == 0:
            return self.field().zero()
        if i != j:
            ret = self.root_of_unity((k.twist() - i.twist() - j.twist()) / 2)
        else:
            i0 = self.one()
            B = self.basis()
            ret = sum(y.ribbon()**2 / (i.ribbon() * x.ribbon()**2)
                   * self.s_ij(i0, y) * self.s_ij(i, z) * self.s_ijconj(x, z)
                   * self.s_ijconj(k, x) * self.s_ijconj(y, z) / self.s_ij(i0, z)
                   for x in B for y in B for z in B) / (self.total_q_order()**4)
        return ret

    def total_q_order(self):
        r"""
        Return the positive square root of :meth:`self.global_q_dimension()
        <global_q_dimension>` as an element of :meth:`self.field() <field>`.

        For the Drinfeld double of a finite group `G`, this equals the
        cardinality of `G`.
        """
        return self._G.order()

    def is_multiplicity_free(self, verbose=False):
        """
        Returns True if all fusion coefficients are at most 1.
        """
        for i in self.basis():
            for j in self.basis():
                for k in self.basis():
                    if self.N_ijk(i,j,k) > 1:
                        if verbose:
                            print ("N(%s,%s,%s)=%s"%(i,j,k,self.N_ijk(i,j,k)))
                        return False
        return True

    def one(self):
        """
        The unit element of the ring.
        """
        return self.basis()[0]

    def dual(self,i):
        sz = self.one()
        for j in self.basis():
            if self.N_ijk(i,j,sz) > 0:
                return j

    def product_on_basis(self, a, b):
        d = {k.support_of_term() : self.N_ijk(self.monomial(a),self.monomial(b),self.dual(k)) for k in self.basis()}
        return self._from_dict(d)

    def _repr_term(self, t):
        return self._names[t]

    def group(self):
        """
        Returns the name of the underlying group.

        EXAMPLE::
            sage: FusionDouble(DiCyclicGroup(4)).group()
            Diyclic group of order 16 as a permutation group
        """
        return self._G

    def IdGroup(self):
        """
        returns the Gap Group ID

        EXAMPLE::
            sage: FusionDouble(DiCyclicGroup(4)).IdGroup()
            [ 16, 9 ]

        """
        return self._G._libgap_().IdGroup()


    class Element(CombinatorialFreeModule.Element):
        def is_simple_object(self):
            r"""
            Determine whether ``self`` is a simple object of the fusion ring.
            """
            return self in self.parent().basis()

        def g(self):
           """
           Returns the conjugacy class representative of the underlying
           group corresponding to a simple object.
           """
           return self.parent()._elt[self.support_of_term()]

        def chi(self):
           """
           Returns the character of the centralizer of a conjugacy class 
           representative of the underlying group corresponding to a simple object.
           """
           return self.parent()._chi[self.support_of_term()]

        def ribbon(self):
            """
            The twist or ribbon of the simple object.

            EXAMPLE::
                sage: H = FusionDouble(CyclicPermutationGroup(3))
                sage: [i.ribbon() for i in H.basis()]
                [1, 1, 1, 1, zeta3, -zeta3 - 1, 1, -zeta3 - 1, zeta3]
            """
            return self.chi()(self.g()) / self.chi()(self.parent()._G.one())

        def twist(self, reduced=True):
            r"""
            Return a rational number `h` such that `\theta = e^{i \pi h}`
            is the twist of ``self``. The quantity `e^{i \pi h}` is
            also available using :meth:`ribbon`.

            This method is only available for simple objects.
            """
            if not self.is_simple_object():
                raise ValueError("quantum twist is only available for simple objects of a FusionRing")
            zeta = self.parent().field().gen()
            rib = self.ribbon()
            for k in range(4*self.parent()._cyclotomic_order):
                if zeta**k == rib:
                    return k/(2*self.parent()._cyclotomic_order)

        def dual(self):
            """
            Returns the dual of self.

            EXAMPLE::
                sage: G = CyclicPermutationGroup(4)
                sage: H = FusionDouble(G, prefix="b")
                sage: [x for x in H.basis() if x==x.dual()]
                [b0, b1, b8, b9]
            """
            return self.parent().dual(self)
            return

# Code below this line should be removed later. 
# Cloned from sage trac repository, commit 1b99dcc .

class FMatrix1:
    """
    Older version of the FMatrix code, cloned from sage
    trac commit 1b99dcc.

    The F-matrix is only determined up to a gauge. It is possible to make
    the F-matrices unitary, or it is possible to make them
    cyclotomic. We choose the latter.

    A good account of the mechanics of finding the F-matrix
    can be found in [Bond2002]_ Section 2.5.

    Due to the large number of equations we may fail to find a
    Groebner basis if there are too many variables. Therefore
    self.get_solution() will refuse to run if there are more
    that 4 fields.

    A1 level 3:
    CPU times: user 1min 58s, sys: 212 ms, total: 1min 59s
    Wall time: 1min 59s

    ..EXAMPLES::

        sage: FR = FusionRing("A1",2)
        sage: FR.fusion_labels(["i0","p","s"],inject_variables=True)
        sage: fmats = FMatrix(FR)
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
        sage: fmats.get_solution()
        Setting up hexagons and pentagons...
        equations: 14
        equations: 41
        Finding a Groebner basis...
        Solving...
        Fixing the gauge...
        adding equation... x1 - 1
        adding equation... x4 - 1
        Done!
        {(p, p, p, p, i0, i0): 1/2*zeta32^12 - 1/2*zeta32^4,
        (p, p, p, p, i0, s): 1,
        (p, p, p, p, s, i0): -1/2,
        (p, p, p, p, s, s): 1/2*zeta32^12 - 1/2*zeta32^4,
        (p, p, s, i0, s, p): 1,
        (p, p, s, s, i0, p): zeta32^12 - zeta32^4,
        (p, s, p, i0, p, p): -1,
        (p, s, p, s, p, p): 1,
        (p, s, s, p, p, i0): 1/2*zeta32^12 - 1/2*zeta32^4,
        (s, p, p, i0, p, s): -1,
        (s, p, p, s, p, i0): -1/2*zeta32^12 + 1/2*zeta32^4,
        (s, p, s, p, p, p): -1,
        (s, s, p, p, i0, p): zeta32^12 - zeta32^4,
        (s, s, s, s, i0, i0): 1}
        sage: F = fmats.assemble_matrices()
        sage: F[p,p,p,p]
        [1/2*zeta32^12 - 1/2*zeta32^4                         -1/2]
        [                           1 1/2*zeta32^12 - 1/2*zeta32^4]

    """
    def __init__(self, fusion_ring, fusion_label="f", var_prefix='fx', max_fields=4):
        self.FR = fusion_ring
        #Set up F-symbols entry by entry
        n_vars = self.findcases()
        self._poly_ring = PolynomialRing(self.FR.field(),n_vars,var_prefix)
        self._poly_ring.inject_variables()
        self._var_to_sextuple, self._fvars = self.findcases(output=True)
        self._max_fields = max_fields

        #Initialize set of defining equations
        self.ideal_basis = set()

        #Initialize empty set of solved F-symbols
        self.solved = set()

    def total_f_symbols(self):
        """
        Count total number of F-symbols. (Initial number of variables)
        """
        return len(self._poly_ring.gens())

    def remaining_vars(self):
        """
        Return a list of unknown F-symbols (reflects current stage of computation)
        """
        return [var for var in self._poly_ring.gens() if var not in self.solved]

    def sreduce_equations(self, nonzeros=None, output=False):
        if nonzeros is None:
            nonzeros = self._poly_ring.gens()
        reduced = [self.sreduce(eq, nonzeros) for eq in self.ideal_basis]
        if output:
            return reduced
        self.ideal_basis = set(reduced)

    #TODO: implement reduction heuristics described in Bonderson p. 37:
    #subsitute for one variable in an equation containing only two terms (perhaps using eqn.specialization)
    #substitute for a variable that appears as a single linear term in an equation
    def reduce_heuristics(self, eqns):
        raise NotImplementedError

    #TODO: implement method. Produce a list of conditions that would make the
    #F-symbols unitary. Based on output parameter, return it or add it to
    #ideal_basis set
    def unitarity_constraints(self, output=True):
        raise NotImplementedError

    def assemble_matrices(self):
        """
        Construct a dictionary mapping a 4-tuple (a, b, c, d) to an F-matrix F_d^{abc}
        """
        ret = dict()
        for a in self.FR.basis():
            for b in self.FR.basis():
                for c in self.FR.basis():
                    for d in self.FR.basis():
                        supp = [self.fmat(a,b,c,d,x,y) for x in self.FR.basis() for y in self.FR.basis() if self.fmat(a,b,c,d,x,y) != 0]

                        if len(supp) > 0:
                            #Assume relevant F-symbols can be assembled into square matrix
                            ret[(a,b,c,d)] = matrix(round(sqrt(len(supp))),supp).transpose()
        return ret

    def fix_gauge(self):
        """
        Fix the gauge by forcing F-symbols not already fixed to equal 1.
        This method should be used AFTER adding hex and pentagon eqns to ideal_basis
        """
        while len(self.solved) < len(self._poly_ring.gens()):
            #Get a variable that has not been fixed
            #In ascending index order, for consistent results
            for var in self._poly_ring.gens():
                if var not in self.solved:
                    break

            #Fix var = 1, substitute, and solve equations
            self.ideal_basis.add(var-1)
            print("adding equation...", var-1)
            self.ideal_basis = set(Ideal(list(self.ideal_basis)).groebner_basis())
            self.substitute_known_values()
            self.update_equations()

    def substitute_known_values(self, eqns=None):
        if not eqns:
            eqns = self.ideal_basis

        for eq in eqns:
            #Ensure polynomial is degree 1 in a single variable
            if sum(eq.degrees()) == 1 and eq.monomials()[0] not in self.solved:
                var = eq.monomials()[0]
                self._fvars[self._var_to_sextuple[var]] = -eq.constant_coefficient()
                #Add variable to set of known values
                self.solved.add(var)

    def update_equations(self):
        """
        Update ideal_basis equations by plugging in known values
        """
        special_values = { known : self._fvars[self._var_to_sextuple[known]] for known in self.solved }
        self.ideal_basis = set(eq.specialization(special_values) for eq in self.ideal_basis)
        self.ideal_basis.discard(0)

    def get_solution(self, factor=False, pent=True, hexa=True):
        """
        Retrieve a set of F-symbols satisfying the hex and pentagon relations.
        """
        if len(self.FR.basis()) <= self._max_fields:
            print("Setting up hexagons and pentagons...")
            eqns = []
            if pent:
                eqns += self.pentagon(factor=factor)
            if hexa:
                eqns += self.hexagon(factor=factor)
            print("Finding a Groebner basis...")
            self.ideal_basis = set(Ideal(eqns).groebner_basis())
            print("Solving...")
            self.substitute_known_values()
            print("Fixing the gauge...")
            self.fix_gauge()
            print("Done!")
            return self._fvars
        else:
            raise NotImplementedError

    def get_ideal(self):
        return Ideal(list(self.ideal_basis))

    def add_equations(self, eqns):
        #TODO: consider replacing set union by ideal intersection. (study computational cost)
        self.ideal_basis = self.ideal_basis.union(set(eqns))

    def clear_equations(self):
        self.ideal_basis = set()

    def clear_vars(self):
        self._fvars = { self._var_to_sextuple[key] : key for key in self._var_to_sextuple }

    def clear(self):
        self.clear_equations()
        self.clear_vars()

    def fmat(self, a, b, c, d, x, y, data=True):
        """
        Set up an entry of the F matrix
        """
        #Determine if fusion tree is admissible
        admissible = self.FR.Nk_ij(a,b,x) * self.FR.Nk_ij(x,c,d) * self.FR.Nk_ij(b,c,y) * self.FR.Nk_ij(a,y,d)

        if admissible == 0:
            return 0

        #Some known zero F-symbols
        if a == self.FR.one():
            if x == b and y == d:
                return 1
            else:
                return 0
        if b == self.FR.one():
            if x == a and y == c:
                return 1
            else:
                return 0
        if c == self.FR.one():
            if x == d and y == b:
                return 1
            else:
                return 0
        if data:
            return self._fvars.get((a,b,c,d,x,y),0)
        else:
            return (a,b,c,d,x,y)

    def findcases(self,output=False):
        """
        Find the unknown F-matrix entries. If run with output=True,
        this returns two dictionaries; otherwise it just returns the
        number of unknown values.
        """
        i = 0
        if output:
            idx_map = dict()
            ret = dict()
        for a in self.FR.basis():
            for b in self.FR.basis():
                for c in self.FR.basis():
                    for d in self.FR.basis():
                        for x in self.FR.basis():
                            for y in self.FR.basis():
                                fm = self.fmat(a, b, c, d, x, y, data=False)
                                if fm is not None and fm not in [0,1]:
                                    if output:
                                        v = self._poly_ring.gens()[i]
                                        ret[(a,b,c,d,x,y)] = v
                                        idx_map[v] = (a, b, c, d, x, y)
                                    i += 1
        if output:
            return idx_map, ret
        else:
            return i

    def singletons(self):
        """
        Find variables that are automatically nonzero, because their F-matrix is 1x1
        """
        ret = []
        for a in self.FR.basis():
            for b in self.FR.basis():
                for c in self.FR.basis():
                    for d in self.FR.basis():
                        supp = [[x,y] for x in self.FR.basis() for y in self.FR.basis() if self.fmat(a,b,c,d,x,y) != 0]
                        if len(supp) == 1:
                            [x,y] = supp[0]
                            if self.fmat(a,b,c,d,x,y) not in [0,1]:
                                ret.append(self.fmat(a,b,c,d,x,y))
        return ret

    def sreduce(self, expr, nonzeros=None):
        """
        From an equation, discard the leading coefficient and
        any factors that are singletons (variables guaranteed to
        be nonzero).
        """
        if nonzeros is None:
            nonzeros = self._poly_ring.gens()
        ret = 1
        for (a,e) in expr.factor()._Factorization__x:
            if a not in nonzeros:
                ret *= a^e
        return ret

    def feq(self, a, b, c, d, e, f, g, k, l):
        """
        Pentagon axiom following Bonderson (2.77)
        """
        lhs = self.fmat(f,c,d,e,g,l)*self.fmat(a,b,l,e,f,k)
        rhs = sum(self.fmat(a,b,c,g,f,h)*self.fmat(a,h,d,e,g,k)*self.fmat(b,c,d,k,h,l) for h in self.FR.basis())
        return lhs - rhs

    def req(self, a, b, c, d, e, g):
        """
        Hexagon axiom following Bonderson (2.78)
        """
        lhs = self.FR.r_matrix(a,c,e)*self.fmat(a,c,b,d,e,g)*self.FR.r_matrix(b,c,g)
        rhs = sum(self.fmat(c,a,b,d,e,f)*self.FR.r_matrix(f,c,d)*self.fmat(a,b,c,d,f,g) for f in self.FR.basis())
        return lhs-rhs

    def hexagon(self, verbose=False, output=True, factor=False):
        """
        Display the information we get from the hexagon relations,
        Bonderson's version
        """
        ret = []
        for a in self.FR.basis():
            for b in self.FR.basis():
                for c in self.FR.basis():
                    for d in self.FR.basis():
                        for e in self.FR.basis():
                            for g in self.FR.basis():
                                rd = self.req(a,b,c,d,e,g)
                                if rd != 0:
                                    if factor:
                                        rd = self.sreduce(rd)
                                    ret.append(rd)
                                    if verbose:
                                        print ("%s,%s,%s,%s,%s,%s : %s"%(a,b,c,d,e,g,rd.factor()))
        print ("equations: %s"%len(ret))
        if output:
            return ret

    def pentagon(self, verbose=False, output=True, factor=None):
        """
        Display the information from the pentagon relations.
        """
        ret = []
        for a in self.FR.basis():
            for b in self.FR.basis():
                for c in self.FR.basis():
                    for d in self.FR.basis():
                        for e in self.FR.basis():
                            for f in self.FR.basis():
                                for g in self.FR.basis():
                                    for k in self.FR.basis():
                                        for l in self.FR.basis():
                                            pd = self.feq(a,b,c,d,e,f,g,k,l)
                                            if pd != 0:
                                                if factor:
                                                    pd = self.sreduce(pd)
                                                ret.append(pd)
                                                if verbose:
                                                    print ("%s,%s,%s,%s,%s,%s,%s,%s,%s : %s"%(a,b,c,d,e,f,g,k,l,pd))
        print ("equations: %s"%len(ret))
        if output:
            return ret

    def equation_graph(self, equations):
        G = graphs.EmptyGraph()
        for e in equations:
            s = [v for v in e.variables()]
            for x in s:
                for y in s:
                    if y!=x:
                        G.add_edge(x,y)
        return(G)

