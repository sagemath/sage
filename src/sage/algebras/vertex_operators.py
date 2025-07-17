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

from sage.categories.cartesian_product import cartesian_product
from sage.categories.action import Action
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.cachefunc import cached_function
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

class FermionicFockSpace(CombinatorialFreeModule):

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

    TESTS::

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

    TESTS::

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
    
    TESTS::
    
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




