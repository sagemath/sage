"""
Space of pseudomorphisms of free modules

AUTHORS:

- Xavier Caruso, Yossef Musleh (2024-09): initial version
"""
# ****************************************************************************
#  Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#                     Yossef Musleh <specialholonomy@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.homset import HomsetWithBase
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.sequence import Sequence
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.modules.free_module_pseudomorphism import FreeModulePseudoMorphism


class FreeModulePseudoHomspace(UniqueRepresentation, HomsetWithBase):
    r"""
    This class implements the space of Pseudomorphisms with a fixed twist.

    For free modules, the elements of a pseudomorphism correspond to matrices
    which define the mapping on elements of a basis.

    TESTS::

        sage: F = GF(125)
        sage: M = F^2
        sage: Frob = F.frobenius_endomorphism()
        sage: PHS = M.pseudoHom(Frob)
        sage: h = PHS([[1, 2], [1, 1]])
        sage: e = M((4*F.gen()^2 + F.gen() + 2, 4*F.gen()^2 + 4*F.gen() + 4))
        sage: h(e)
        (z3, 2*z3^2 + 3*z3 + 3)
    """
    Element = FreeModulePseudoMorphism

    @staticmethod
    def __classcall_private__(cls, domain, codomain, twist):
        r"""
        Constructs the space of pseudomorphisms with a given twist.

        INPUT:

        -  ``domain`` -- a free module,  the domain of this pseudomorphism

        -  ``codomain`` -- a free module, the codomain of this pseudomorphism

        -  ``twist`` -- a twisting morphism/derivation or a Ore polynomial ring

        TESTS::

            sage: F = GF(125)
            sage: Frob = F.frobenius_endomorphism()
            sage: M = F^2
            sage: H = M.pseudoHom(Frob)
            sage: type(H)
            <class 'sage.modules.free_module_pseudohomspace.FreeModulePseudoHomspace_with_category'>

            sage: TestSuite(H).run()
        """
        ring = domain.base_ring()
        if codomain.base_ring() is not ring:
            raise ValueError("the domain and the codomain must be defined over the same ring")
        if isinstance(twist, OrePolynomialRing):
            ore = twist
            if ore.base_ring() is not ring:
                raise ValueError("base rings do not match")
        else:
            ore = OrePolynomialRing(ring, twist, names='x', polcast=False)
        if ore._derivation is not None:
            if not codomain.has_coerce_map_from(domain):
                raise ValueError("the domain does not coerce into the codomain")
        return cls.__classcall__(cls, domain, codomain, ore)

    def __init__(self, domain, codomain, ore):
        r"""
        Initialize this pseudohom space.

        INPUT:

        -  ``domain`` -- a free module,  the domain of this pseudomorphism

        -  ``codomain`` -- a free module, the codomain of this pseudomorphism

        -  ``ore`` -- the underlying Ore polynomial ring

        TESTS::

            sage: F = GF(125)
            sage: Frob = F.frobenius_endomorphism()
            sage: M = F^2
            sage: M.pseudoHom(Frob)
            Set of Pseudoendomorphisms (twisted by z3 |--> z3^5) of
            Vector space of dimension 2 over Finite Field in z3 of size 5^3
        """
        self._domain = domain
        self._codomain = codomain
        super().__init__(domain, codomain, category=None)
        self._ore = ore
        if isinstance(ore, OrePolynomialRing):
            self._morphism = ore._morphism
            self._derivation = ore._derivation
        else:
            self._morphism = self._derivation = None
        ring = ore.base_ring()
        self._matrix_space = MatrixSpace(ring, domain.dimension(), codomain.dimension())

    def _element_constructor_(self, f, side="left"):
        r"""
        Return the element of this parent constructed from the
        given data.

        TESTS::

            sage: F.<z> = GF(5^3)
            sage: Frob = F.frobenius_endomorphism()
            sage: V = F^2
            sage: H = V.pseudoHom(Frob)

            sage: H([[1, z], [z, z^2]])
            Free module pseudomorphism (twisted by z |--> z^5) defined by the matrix
            [  1   z]
            [  z z^2]
            Domain: Vector space of dimension 2 over Finite Field in z of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 5^3
        """
        return self.element_class(self, f, side)

    def __reduce__(self):
        r"""
        TESTS::

            sage: F = GF(125)
            sage: Frob = F.frobenius_endomorphism()
            sage: M = F^2
            sage: H = M.pseudoHom(Frob)
            sage: loads(dumps(M)) is M
            True
        """
        if self._derivation is None:
            twist = self._morphism
        else:
            twist = self._derivation
        return FreeModulePseudoHomspace, (self.domain(), self.codomain(), twist)

    def _repr_(self):
        r"""
        Returns a string representation of this pseudomorphism space.

        EXAMPLES::

            sage: Fq = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: V.pseudoHom(Frob)  # indirect doctest
            Set of Pseudoendomorphisms (twisted by z3 |--> z3^7) of
            Vector space of dimension 2 over Finite Field in z3 of size 7^3

        ::

            sage: V.pseudoHom(Frob, codomain=Fq^3)  # indirect doctest
            Set of Pseudomorphism (twisted by z3 |--> z3^7)
            from Vector space of dimension 2 over Finite Field in z3 of size 7^3
            to Vector space of dimension 3 over Finite Field in z3 of size 7^3

        ::

            sage: A.<t> = QQ[]
            sage: d = A.derivation()
            sage: M = A^3
            sage: M.pseudoHom(d)
            Set of Pseudoendomorphisms (twisted by d/dt) of Ambient free module of rank 3 over
            the principal ideal domain Univariate Polynomial Ring in t over Rational Field
        """
        twist = self._ore._repr_twist()
        if self.domain() is self.codomain():
            return "Set of Pseudoendomorphisms (%s) of %s" % (twist, self.domain())
        else:
            return "Set of Pseudomorphism (%s) from %s to %s" % (twist, self.domain(), self.codomain())

    def ore_ring(self, var='x'):
        r"""
        Return the underlying Ore polynomial ring.

        INPUT:

        - ``var`` -- string (default: ``x``) the name of
          tha variable

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: H = V.pseudoHom(Frob)

            sage: H.ore_ring()
            Ore Polynomial Ring in x over Finite Field in z of size 7^3 twisted by z |--> z^7

            sage: H.ore_ring('y')
            Ore Polynomial Ring in y over Finite Field in z of size 7^3 twisted by z |--> z^7
        """
        return self._ore.change_var(var)

    def matrix_space(self):
        r"""
        Return the matrix space used for representing the
        pseudomorphism in this space.

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: W = Fq^3
            sage: H = V.pseudoHom(Frob, codomain=W)
            sage: H.matrix_space()
            Full MatrixSpace of 2 by 3 dense matrices over Finite Field in z of size 7^3
        """
        return self._matrix_space

    def basis(self, side="left"):
        r"""
        Return a basis for the underlying matrix space.

        EXAMPLES::

            sage: Fq = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: PHS = V.pseudoHom(Frob)
            sage: PHS.basis()
            [Free module pseudomorphism (twisted by z3 |--> z3^7) defined by the matrix
            [1 0]
            [0 0]
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3,
            Free module pseudomorphism (twisted by z3 |--> z3^7) defined by the matrix
            [0 1]
            [0 0]
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3,
            Free module pseudomorphism (twisted by z3 |--> z3^7) defined by the matrix
            [0 0]
            [1 0]
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3,
            Free module pseudomorphism (twisted by z3 |--> z3^7) defined by the matrix
            [0 0]
            [0 1]
            Domain: Vector space of dimension 2 over Finite Field in z3 of size 7^3
            Codomain: Vector space of dimension 2 over Finite Field in z3 of size 7^3]
        """
        return Sequence(self(mat) for mat in self._matrix_space.basis())

    def _test_additive_associativity(self, tester):
        r"""
        Test associativity for (not necessarily all) elements in this parent.

        This test is not relevant for pseudo-morphisms because they are not
        stable by addition.

        TESTS::

            sage: Fq = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: PHS = V.pseudoHom(Frob)
            sage: TestSuite(PHS).run()  # indirect doctest
        """
        pass

    def _test_distributivity(self, tester):
        r"""
        Test distributivity for (not necessarily all) elements in this parent.

        This test is not relevant for pseudo-morphisms because they are not
        stable by addition.

        TESTS::

            sage: Fq = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: PHS = V.pseudoHom(Frob)
            sage: TestSuite(PHS).run()  # indirect doctest
        """
        pass

    def _test_one(self, tester):
        r"""
        Test properties the identity element.

        This test is not relevant for pseudo-morphisms because the identity
        is not a pseudo-morphism in general.

        TESTS::

            sage: Fq = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: PHS = V.pseudoHom(Frob)
            sage: TestSuite(PHS).run()  # indirect doctest
        """
        pass

    def _test_zero(self, tester):
        r"""
        Test properties of the zero element.

        This test is not relevant for pseudo-morphisms because the zero
        map is not a pseudo-morphism in general.

        TESTS::

            sage: Fq = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: PHS = V.pseudoHom(Frob)
            sage: TestSuite(PHS).run()  # indirect doctest
        """
        pass
