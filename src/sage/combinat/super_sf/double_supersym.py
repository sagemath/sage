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
    def __init__(self, Sym):
        self.Sym = Sym

    def homogeneous(self):
        return

    def elementary(self):
        return

    def schur(self):
        return

