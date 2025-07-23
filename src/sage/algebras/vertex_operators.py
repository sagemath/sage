r"""
Currently just a sandbox for experimenting with vertex operator code
before a more permanent implementation

TESTS::
    sage: from sage.algebras.vertex_operators import *
    sage: R = SymmetricFunctions(QQ); p = R.p(); s = R.s()
    sage: P = p.completion()
    sage: F = FermionicFockSpace(s)
    sage: x = 3*F(([3,2,1],0)); x
    3*s[]*|[3, 2, 1], 0>
    sage: J = Current()
    sage: J.act(3, x)
    (-3*s[])*|[1, 1, 1], 0> + (-3*s[])*|[3], 0>
    sage: deg(x)
    6
"""

from sage.algebras.weyl_algebra import DifferentialWeylAlgebra
from sage.categories.cartesian_product import cartesian_product
from sage.categories.action import Action
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partitions, Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.data_structures.stream import Stream_function, Stream_cauchy_compose
from sage.functions.other import factorial
from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod
from sage.rings.infinity import PlusInfinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.non_negative_integers import NonNegativeIntegers

class FermionicFockSpace(CombinatorialFreeModule):
    r"""
    A model of the Fermionic Fock space `\mathcal{H}_F`. As a vector space, `\mathcal{H}_F`
     has basis given by pairs `(\lambda, n)` where `\lambda` is a partition, and `n` is an integer.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: F = FermionicFockSpace(QQ)
        sage: F.an_element()
        |[], 0> + 4*|[], 1> + 2*|[1], 0>
    """
    def __init__(self, R):
        I = cartesian_product((Partitions(), ZZ))
        CombinatorialFreeModule.__init__(self, R, I, prefix='', bracket='')

    def _repr_term(self, m):
        return '|' + str(m[0]) + ', ' + str(m[1]) + '>'

class Current(Action):
    """
    The action of the current operators `J_i` on the Fermionic Fock space.
     For `i \neq 0`, we can view `J_i` as moving a particle by `-i` in all possible ways.
    Equivalently, viewed on the bosonic side, the action of `J_i` can be calculated
     using the Murnaghan-Nakayam rule.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: R = SymmetricFunctions(QQ); s = R.s()
        sage: F = FermionicFockSpace(s)
        sage: x = F(([3,2], 0))
        sage: J = Current()
        sage: J.act(-3, x)
        s[]*|[3, 2, 1, 1, 1], 0> + (-s[])*|[3, 2, 2, 1], 0> + (-s[])*|[4, 4], 0>
         + s[]*|[6, 2], 0>
        sage: J.act(1,x)
        s[]*|[2, 2], 0> + s[]*|[3, 1], 0>

    """
    def __init__(self):
        self._R = SymmetricFunctions(QQ)
        self._s = self._R.s()
        self._p = self._R.p()
        self.fockspace = FermionicFockSpace(self._s)
        super().__init__(ZZ, self.fockspace)
    def _act_(self, g, x):
        res = self.fockspace.zero()
        if g > 0: # \partial p_k
            for (m,c) in x.monomial_coefficients().items():
                par = self._s(m[0])
                for (m2,c2) in (par.skew_by(self._p[g])).monomial_coefficients().items():
                    res += c*c2*self.fockspace((m2, m[1]))
            return res
        elif g < 0: #-k*p_k
            for (m,c) in x.monomial_coefficients().items():
                par = self._s(m[0])
                for (m2, c2) in (self._s(self._p[-g]*par)).monomial_coefficients().items():
                    res += c*c2*self.fockspace((m2, m[1]))     
            return res
        else: #multiply by charge
            return sum(m[1]*c*self.fockspace(m) for (m,c) in x.monomial_coefficients().items())
    def act_by_partition(self,lam, sign, x):
        for i in lam:
            # TODO: Make this more efficient by acting by the whole partition at once
            x = self._act_(i*sign, x) 
        return x



@cached_function
def Hamiltonian():
    r"""

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: H = Hamiltonian()
        sage: H.truncate(3)
        p[] + p[1] + (1/2*p[1,1]+1/2*p[2])
    """
    p = SymmetricFunctions(QQ).p()
    P = p.completion()
    return P(lambda n: p[n]/n, valuation=1).exp()

def act_by_H(x):
    r"""
    Computes the action of H on an element ``x`` of the Fermionic Fock Space

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: R = SymmetricFunctions(QQ); s = R.s()
        sage: F = FermionicFockSpace(s)
        sage: act_by_H(F(([2,1],0)))
        (s[2,1])*|[], 0> + (s[1,1]+s[2])*|[1], 0> + s[1]*|[1, 1], 0> + s[1]*|[2], 0> +
         s[]*|[2, 1], 0>
    """
    J = Current()
    H = Hamiltonian().truncate(deg(x) + 1).symmetric_function()
    p = J._p

    res = x.parent().zero()
    for (hm, hc) in H.monomial_coefficients().items():
        res += hc*p(hm)*J.act_by_partition(hm,1, x)
    return res

def H_mat_coeff(lam, mu):
    r"""
    
    EXAMPLES::
    
        sage: from sage.algebras.vertex_operators import *
        sage: lam = [3,2,1]; s = SymmetricFunctions(QQ).s()
        sage: all(H_mat_coeff(lam, mu) == s(lam).skew_by(s(mu)) for mu in Partitions(outer=lam))
        True
    """
    R = SymmetricFunctions(QQ)
    s = R.s()
    F = FermionicFockSpace(s)
    return act_by_H(F((lam, 0))).coefficient((mu,0))

# degree of element of F
def deg(x):
    return max((m[0].size() for (m,_) in x.monomial_coefficients().items()), default=0)



class HalfVertexOperator():
    def __init__(self,f):
        self._stream = Stream_cauchy_compose(
            Stream_function(lambda n: ZZ(1)/factorial(n), False, 0),
            Stream_function(f, False, 1), 
            False
        )
    def __getitem__(self, i):
        if i in NonNegativeIntegers():
            return self._stream[i]
        raise ValueError("Invalid input")

class VertexOperator(Action):
    """
    The action of a Vertex Operator on the Bosonic Fock space.

    INPUT:

    - ``pos`` -- function taking in nonnegative integers indexing the positive
     half of the vertex operator
    - ``neg`` -- function taking in nonnegative integers indexing the negative
     half of the vertex operator
    - ``cutoff`` -- function taking in symmetric functions which determines how
     far to expand the vertex operator (default: ``lambda x: max(x.degree(), 1))
    - ``dcharge`` -- integer (default: ``1``) indicating how this vertex operator
     should change the charge of an element 
    - ``fockspace``-- the space that the vertex operators are acting on.

    """
    def __init__(self, pos, neg, cutoff = lambda x: max(x.degree(), 1), dcharge=1, fockspace=None):
        self.pos = HalfVertexOperator(pos)
        self.neg = HalfVertexOperator(neg)
        if fockspace is None:
            self.fockspace = LaurentPolynomialRing(SymmetricFunctions(QQ).s(), names = ('w')) 
        else: #TODO: check input is correct type
            self.fockspace = fockspace
        self.cutoff = cutoff
        self.dcharge = dcharge
        super().__init__(ZZ, self.fockspace)
        
    def _act_(self, i,x):
        res = 0
        # for each component of constant charge, compute the action 
        for (charge, fn) in x.monomial_coefficients().items():
            res += self.fockspace.gen()**(charge+self.dcharge)*self._act_on_sym(i, fn, charge)

        return res

    def _act_on_sym(self, i, x, c):
        p = self.fockspace.base_ring().symmetric_function_ring().p()
        op = self._get_operator(i, x, c)
        if op in self.fockspace.base_ring().base_ring():
            return op*x
        res = 0
        for (m, c) in op.monomial_coefficients().items():
             par1 = self._weyl_to_par(m[0].dict())
             par2 = self._weyl_to_par(m[1].dict())
             for (m2, c2) in x.monomial_coefficients().items():
                res += c*c2*(self.fockspace.base_ring()(m2).skew_by(p(par2)) * p(par1)/(prod(par1)))
        return res
    def _get_operator(self, i, x, c):
        raise NotImplementedError("Use a subclass of VertexOperator")
    
    def _weyl_to_par(self,m):
        res = []
        for i in m:
            res += [i]*int(m[i])
        return Partition(sorted(res, reverse=True))
    

class CreationOperator(VertexOperator):
    """
    
    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: B.<w> = LaurentPolynomialRing(SymmetricFunctions(QQ).s())
        sage: Cre = CreationOperator(B)
        sage: Cre.act(-1, B.one())
        0
        sage: Cre.act(0, B.one())
        s[]*w
        sage: Cre.act(1, w)
        s[]*w^2
        sage: Cre.act(0, w^-1)
        s[1]
        sage: t = Cre.act(1, Cre.act(-1, w^-2)); t
        s[2, 1]
        sage: Cre.act(3, t)
        s[3, 2, 1]*w
    """
    def __init__(self, fockspace):
        self.weyl_algebra = DifferentialWeylAlgebra(QQ, names=('x'), n=PlusInfinity())
        self.x, self.dx = self.weyl_algebra.gens()
        super().__init__(lambda n: self.x[n], lambda n: -self.dx[n]/n, dcharge=1, fockspace=fockspace)

    def _get_operator(self, i, x, c):
        op = 0 
        for j in range(self.cutoff(x)+1):
            if i + j - c < 0:
                continue
            op += self.pos[i + j - c]*self.neg[j]
        return op


# TODO: Fix mysterious off by one error causing tests to fail.
class AnnihilationOperator(VertexOperator):
    """
    
    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: B.<w> = LaurentPolynomialRing(SymmetricFunctions(QQ).s())
        sage: Ann = AnnihilationOperator(B)
        sage: Ann.act(0, B.one()) # known bug
        0
        sage: Ann.act(-1, B.one()) # known bug
        s[]*w^-1
        sage: Ann.act(-2, w^-1) # known bug
        s[]*w^-2
        sage: Ann.act(-1, w) # known bug
        s[1]
        sage: Ann.act(-2, Ann.act(0, w^2)) # known bug
        s[2, 1]
    """
    def __init__(self, fockspace):
       self.weyl_algebra = DifferentialWeylAlgebra(QQ, names=('x'), n=PlusInfinity())
       self.x, self.dx = self.weyl_algebra.gens()
       super().__init__(lambda n: -self.x[n], lambda n: self.dx[n]/n, dcharge=-1, fockspace=fockspace)

    def _get_operator(self, i, x, c):
        op = 0 
        for j in range(self.cutoff(x)+1):
            if -i + j + c < 0:
                continue
            op += self.pos[-i + j + c]*self.neg[j]
        return op