r"""
Double Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.rings import Rings
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.free_module import CombinatorialFreeModule

class DoubleSupersymmetricFunctions(UniqueRepresentation, Parent):
    def __init__(self, R):
        if R not in Rings().Commutative():
            raise ValueError("the base ring must be a commutative ring")
        cat = GradedHopfAlgebras(R).Commutative().Cocommutative()
        if R in PrincipalIdealDomains():
            cat &= UniqueFactorizationDomains()
        Parent.__init__(self, base=R, category=cat.WithRealizations())

class DoubleSupersymmetricFunctions_generic(CombinatorialFreeModule):
    pass