r"""
Morphisms between modular abelian varieties, including Hecke operators acting on modular abelian varieties

Sage can compute with Hecke operators on modular abelian varieties.
A Hecke operator is defined by given a modular abelian variety and
an index. Given a Hecke operator, Sage can compute the
characteristic polynomial, and the action of the Hecke operator on
various homology groups.

AUTHORS:

- William Stein (2007-03)

- Craig Citro (2008-03)

EXAMPLES::

    sage: A = J0(54)
    sage: t5 = A.hecke_operator(5); t5
    Hecke operator T_5 on Abelian variety J0(54) of dimension 4
    sage: t5.charpoly().factor()
    (x - 3) * (x + 3) * x^2
    sage: B = A.new_subvariety(); B
    Abelian subvariety of dimension 2 of J0(54)
    sage: t5 = B.hecke_operator(5); t5
    Hecke operator T_5 on Abelian subvariety of dimension 2 of J0(54)
    sage: t5.charpoly().factor()
    (x - 3) * (x + 3)
    sage: t5.action_on_homology().matrix()
    [ 0  3  3 -3]
    [-3  3  3  0]
    [ 3  3  0 -3]
    [-3  6  3 -3]
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.categories.morphism import Morphism as base_Morphism
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

import sage.modules.matrix_morphism
import sage.matrix.matrix_space as matrix_space

from .finite_subgroup import TorsionPoint

class Morphism_abstract(sage.modules.matrix_morphism.MatrixMorphism_abstract):
    """
    A morphism between modular abelian varieties.

    EXAMPLES::

        sage: t = J0(11).hecke_operator(2)
        sage: from sage.modular.abvar.morphism import Morphism
        sage: isinstance(t, Morphism)
        True
    """

    def _repr_(self):
        r"""
        Return string representation of this morphism.

        EXAMPLES::

            sage: t = J0(11).hecke_operator(2)
            sage: sage.modular.abvar.morphism.Morphism_abstract._repr_(t)
            'Abelian variety endomorphism of Abelian variety J0(11) of dimension 1'
            sage: J0(42).projection(J0(42)[0])._repr_()
            'Abelian variety morphism:\n  From: Abelian variety J0(42) of dimension 5\n  To:   Simple abelian subvariety 14a(1,42) of dimension 1 of J0(42)'
        """
        return base_Morphism._repr_(self)

    def _repr_type(self):
        """
        Return type of morphism.

        EXAMPLES::

            sage: t = J0(11).hecke_operator(2)
            sage: sage.modular.abvar.morphism.Morphism_abstract._repr_type(t)
            'Abelian variety'
        """
        return "Abelian variety"

    def complementary_isogeny(self):
        """
        Returns the complementary isogeny of self.

        EXAMPLES::

            sage: J = J0(43)
            sage: A = J[1]
            sage: T5 = A.hecke_operator(5)
            sage: T5.is_isogeny()
            True
            sage: T5.complementary_isogeny()
            Abelian variety endomorphism of Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: (T5.complementary_isogeny() * T5).matrix()
            [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]
        """
        if not self.is_isogeny():
            raise ValueError("self is not an isogeny")
        M = self.matrix()
        try:
            iM, denom = M._invert_iml()
        except AttributeError:
            iM = M.matrix_over_field().invert()
            iM, denom = iM._clear_denom()
        return Morphism(self.parent().reversed(), iM)

    def is_isogeny(self):
        """
        Return True if this morphism is an isogeny of abelian varieties.

        EXAMPLES::

            sage: J = J0(39)
            sage: Id = J.hecke_operator(1)
            sage: Id.is_isogeny()
            True
            sage: J.hecke_operator(19).is_isogeny()
            False
        """
        M = self.matrix()
        return M.nrows() == M.ncols() == M.rank()

    def cokernel(self):
        """
        Return the cokernel of self.

        OUTPUT:


        -  ``A`` - an abelian variety (the cokernel)

        -  ``phi`` - a quotient map from self.codomain() to the
           cokernel of self


        EXAMPLES::

            sage: t = J0(33).hecke_operator(2)
            sage: (t-1).cokernel()
            (Abelian subvariety of dimension 1 of J0(33),
             Abelian variety morphism:
              From: Abelian variety J0(33) of dimension 3
              To:   Abelian subvariety of dimension 1 of J0(33))

        Projection will always have cokernel zero.

        ::

            sage: J0(37).projection(J0(37)[0]).cokernel()
            (Simple abelian subvariety of dimension 0 of J0(37),
             Abelian variety morphism:
              From: Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37)
              To:   Simple abelian subvariety of dimension 0 of J0(37))

        Here we have a nontrivial cokernel of a Hecke operator, as the
        T_2-eigenvalue for the newform 37b is 0.

        ::

            sage: J0(37).hecke_operator(2).cokernel()
            (Abelian subvariety of dimension 1 of J0(37),
             Abelian variety morphism:
              From: Abelian variety J0(37) of dimension 2
              To:   Abelian subvariety of dimension 1 of J0(37))
            sage: AbelianVariety('37b').newform().q_expansion(5)
            q + q^3 - 2*q^4 + O(q^5)
        """
        try:
            return self.__cokernel
        except AttributeError:
            I = self.image()
            C = self.codomain().quotient(I)
            self.__cokernel = C
            return C


    def kernel(self):
        """
        Return the kernel of this morphism.

        OUTPUT:


        -  ``G`` - a finite group

        -  ``A`` - an abelian variety (identity component of
           the kernel)


        EXAMPLES: We compute the kernel of a projection map. Notice that
        the kernel has a nontrivial abelian variety part.

        ::

            sage: A, B, C = J0(33)
            sage: pi = J0(33).projection(B)
            sage: pi.kernel()
            (Finite subgroup with invariants [20] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 2 of J0(33))

        We compute the kernels of some Hecke operators::

            sage: t2 = J0(33).hecke_operator(2)
            sage: t2
            Hecke operator T_2 on Abelian variety J0(33) of dimension 3
            sage: t2.kernel()
            (Finite subgroup with invariants [2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 0 of J0(33))
            sage: t3 = J0(33).hecke_operator(3)
            sage: t3.kernel()
            (Finite subgroup with invariants [3, 3] over QQ of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 0 of J0(33))
        """
        A = self.matrix()
        L = A.image().change_ring(ZZ)
        # Saturate the image of the matrix corresponding to self.
        Lsat = L.saturation()
        # Now find a matrix whose rows map exactly onto the
        # saturation of L.
        X = A.solve_left(Lsat.basis_matrix())
        D = self.domain()
        V = (A.kernel().basis_matrix() * D.vector_space().basis_matrix()).row_module()
        Lambda = V.intersection(D._ambient_lattice())
        from .abvar import ModularAbelianVariety
        abvar = ModularAbelianVariety(D.groups(), Lambda, D.base_ring())

        if Lambda.rank() == 0:
            field_of_definition = QQ
        else:
            field_of_definition = None

        lattice = (X * self.domain().lattice().basis_matrix()).row_module(ZZ)

        K = D.finite_subgroup(lattice, field_of_definition=field_of_definition)

        return K, abvar


    def factor_out_component_group(self):
        r"""
        View self as a morphism `f:A \to B`. Then `\ker(f)`
        is an extension of an abelian variety `C` by a finite
        component group `G`. This function constructs a morphism
        `g` with domain `A` and codomain Q isogenous to
        `C` such that `\ker(g)` is equal to `C`.

        OUTPUT: a morphism

        EXAMPLES::

            sage: A,B,C = J0(33)
            sage: pi = J0(33).projection(A)
            sage: pi.kernel()
            (Finite subgroup with invariants [5] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 2 of J0(33))
            sage: psi = pi.factor_out_component_group()
            sage: psi.kernel()
            (Finite subgroup with invariants [] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 2 of J0(33))

        ALGORITHM: We compute a subgroup `G` of `B` so that
        the composition `h: A\to B \to B/G` has kernel that
        contains `A[n]` and component group isomorphic to
        `(\ZZ/n\ZZ)^{2d}`, where `d` is the
        dimension of `A`. Then `h` factors through
        multiplication by `n`, so there is a morphism
        `g: A\to B/G` such that `g \circ [n] = h`. Then
        `g` is the desired morphism. We give more details below
        about how to transform this into linear algebra.
        """
        try:
            return self.__factor_out
        except AttributeError:
            A = self.matrix()
            L = A.image()
            # Saturate the image of the matrix corresponding to self.
            Lsat = L.saturation()
            if L == Lsat: # easy case
                self.__factor_out = self
                return self
            # Now find a matrix whose rows map exactly onto the
            # saturation of L.
            X = A.solve_left(Lsat.basis_matrix())

            # Find an integer n such that 1/n times the lattice Lambda of A
            # contains the row span of X.
            n = X.denominator()

            # Let M be the lattice of the codomain B of self.
            # Now 1/n * Lambda contains Lambda and maps
            # via the matrix of self to a lattice L' that
            # contains Lsat.   Consider the lattice
            #                R = M + L'.
            # This is a lattice that contains the lattice M of B.
            # Also 1/n*Lambda maps exactly to L' in R.
            # We have
            #     R/L' = (M+L')/L' = M/(L'/\M) = M/Lsat
            # which is torsion free!

            Q      = self.codomain()
            M      = Q.lattice()
            one_over_n = ZZ(1)/n
            Lprime = (one_over_n * self.matrix() * M.basis_matrix()).row_module(ZZ)

            # This R is a lattice in the ambient space for B.
            R = Lprime + M

            from .abvar import ModularAbelianVariety
            C = ModularAbelianVariety(Q.groups(), R, Q.base_field())

            # We have to change the basis of the representation of A
            # to the basis for R instead of the basis for M.  Each row
            # of A is written in terms of M, but needs to be in terms
            # of R's basis, which contains M with finite index.
            change_basis_from_M_to_R = R.basis_matrix().solve_left(M.basis_matrix())
            matrix = one_over_n * A * change_basis_from_M_to_R

            # Finally
            g = Morphism(self.domain().Hom(C), matrix)
            self.__factor_out = g
            return g

    def image(self):
        """
        Return the image of this morphism.

        OUTPUT: an abelian variety

        EXAMPLES: We compute the image of projection onto a factor of
        `J_0(33)`::

            sage: A,B,C = J0(33)
            sage: A
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: f = J0(33).projection(A)
            sage: f.image()
            Abelian subvariety of dimension 1 of J0(33)
            sage: f.image() == A
            True

        We compute the image of a Hecke operator::

            sage: t2 = J0(33).hecke_operator(2); t2.fcp()
            (x - 1) * (x + 2)^2
            sage: phi = t2 + 2
            sage: phi.image()
            Abelian subvariety of dimension 1 of J0(33)

        The sum of the image and the kernel is the whole space::

            sage: phi.kernel()[1] + phi.image() == J0(33)
            True
        """
        return self(self.domain())

    def __call__(self, X):
        """
        INPUT:


        -  ``X`` - abelian variety, finite group, or torsion
           element


        OUTPUT: abelian variety, finite group, torsion element

        EXAMPLES: We apply morphisms to elements::

            sage: t2 = J0(33).hecke_operator(2)
            sage: G  = J0(33).torsion_subgroup(2); G
            Finite subgroup with invariants [2, 2, 2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2(G.0)
            [(-1/2, 0, 1/2, -1/2, 1/2, -1/2)]
            sage: t2(G.0) in G
            True
            sage: t2(G.1)
            [(0, -1, 1/2, 0, 1/2, -1/2)]
            sage: t2(G.2)
            [(0, 0, 0, 0, 0, 0)]
            sage: K = t2.kernel()[0]; K
            Finite subgroup with invariants [2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2(K.0)
            [(0, 0, 0, 0, 0, 0)]

        We apply morphisms to subgroups::

            sage: t2 = J0(33).hecke_operator(2)
            sage: G  = J0(33).torsion_subgroup(2); G
            Finite subgroup with invariants [2, 2, 2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2(G)
            Finite subgroup with invariants [2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2.fcp()
            (x - 1) * (x + 2)^2

        We apply morphisms to abelian subvarieties::

            sage: E11a0, E11a1, B = J0(33)
            sage: t2 = J0(33).hecke_operator(2)
            sage: t3 = J0(33).hecke_operator(3)
            sage: E11a0
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: t3(E11a0)
            Abelian subvariety of dimension 1 of J0(33)
            sage: t3(E11a0).decomposition()
            [
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            ]
            sage: t3(E11a0) == E11a1
            True
            sage: t2(E11a0) == E11a0
            True
            sage: t3(E11a0) == E11a0
            False
            sage: t3(E11a0 + E11a1) == E11a0 + E11a1
            True

        We apply some Hecke operators to the cuspidal subgroup and split it
        up::

            sage: C = J0(33).cuspidal_subgroup(); C
            Finite subgroup with invariants [10, 10] over QQ of Abelian variety J0(33) of dimension 3
            sage: t2 = J0(33).hecke_operator(2); t2.fcp()
            (x - 1) * (x + 2)^2
            sage: (t2 - 1)(C)
            Finite subgroup with invariants [5, 5] over QQ of Abelian variety J0(33) of dimension 3
            sage: (t2 + 2)(C)
            Finite subgroup with invariants [2, 2] over QQ of Abelian variety J0(33) of dimension 3

        Same but on a simple new factor::

            sage: C = J0(33)[2].cuspidal_subgroup(); C
            Finite subgroup with invariants [2, 2] over QQ of Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            sage: t2 = J0(33)[2].hecke_operator(2); t2.fcp()
            x - 1
            sage: t2(C)
            Finite subgroup with invariants [2, 2] over QQ of Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
        """
        from .abvar import is_ModularAbelianVariety
        from .finite_subgroup import FiniteSubgroup
        if isinstance(X, TorsionPoint):
            return self._image_of_element(X)
        elif is_ModularAbelianVariety(X):
            return self._image_of_abvar(X)
        elif isinstance(X, FiniteSubgroup):
            return self._image_of_finite_subgroup(X)
        else:
            raise TypeError("X must be an abelian variety or finite subgroup")

    def _image_of_element(self, x):
        """
        Return the image of the torsion point `x` under this
        morphism.

        The parent of the image element is always the group of all torsion
        elements of the abelian variety.

        INPUT:


        -  ``x`` - a torsion point on an abelian variety


        OUTPUT: a torsion point

        EXAMPLES::

            sage: A = J0(11); t = A.hecke_operator(2)
            sage: t.matrix()
            [-2  0]
            [ 0 -2]
            sage: P = A.cuspidal_subgroup().0; P
            [(0, 1/5)]
            sage: t._image_of_element(P)
            [(0, -2/5)]
            sage: -2*P
            [(0, -2/5)]

        ::

            sage: J = J0(37) ; phi = J._isogeny_to_product_of_simples()
            sage: phi._image_of_element(J.torsion_subgroup(5).gens()[0])
            [(1/5, -1/5, -1/5, 1/5, 1/5, 1/5, 1/5, -1/5)]

        ::

            sage: K = J[0].intersection(J[1])[0]
            sage: K.list()
            [[(0, 0, 0, 0)],
             [(1/2, -1/2, 1/2, 0)],
             [(0, 0, 1, -1/2)],
             [(1/2, -1/2, 3/2, -1/2)]]
            sage: [ phi.restrict_domain(J[0])._image_of_element(k) for k in K ]
            [[(0, 0, 0, 0, 0, 0, 0, 0)],
             [(0, 0, 0, 0, 0, 0, 0, 0)],
             [(0, 0, 0, 0, 0, 0, 0, 0)],
             [(0, 0, 0, 0, 0, 0, 0, 0)]]
        """
        v = x._relative_element() * self.matrix() * self.codomain().lattice().basis_matrix()
        T = self.codomain().qbar_torsion_subgroup()
        return T(v)

    def _image_of_finite_subgroup(self, G):
        """
        Return the image of the finite group `G` under ``self``.

        INPUT:

        - ``G`` -- a finite subgroup of the domain of ``self``

        OUTPUT:

        A finite subgroup of the codomain.

        EXAMPLES::

            sage: J = J0(33); A = J[0]; B = J[1]
            sage: C = A.intersection(B)[0]; C
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: t = J.hecke_operator(3)
            sage: D = t(C); D
            Finite subgroup with invariants [5] over QQ of Abelian variety J0(33) of dimension 3
            sage: D == C
            True

        Or we directly test this function::

            sage: D = t._image_of_finite_subgroup(C); D
            Finite subgroup with invariants [5] over QQ of Abelian variety J0(33) of dimension 3
            sage: phi = J0(11).degeneracy_map(22,2)
            sage: J0(11).rational_torsion_subgroup().order()
            5
            sage: phi._image_of_finite_subgroup(J0(11).rational_torsion_subgroup())
            Finite subgroup with invariants [5] over QQ of Abelian variety J0(22) of dimension 2
        """
        B = G._relative_basis_matrix() * self.restrict_domain(G.abelian_variety()).matrix() * self.codomain().lattice().basis_matrix()
        lattice = B.row_module(ZZ)
        return self.codomain().finite_subgroup(lattice,
                             field_of_definition = G.field_of_definition())

    def _image_of_abvar(self, A):
        """
        Compute the image of the abelian variety `A` under this
        morphism.

        INPUT:


        -  ``A`` - an abelian variety


        OUTPUT an abelian variety

        EXAMPLES::

            sage: t = J0(33).hecke_operator(2)
            sage: t._image_of_abvar(J0(33).new_subvariety())
            Abelian subvariety of dimension 1 of J0(33)

        ::

            sage: t = J0(33).hecke_operator(3)
            sage: A = J0(33)[0]
            sage: B = t._image_of_abvar(A); B
            Abelian subvariety of dimension 1 of J0(33)
            sage: B == A
            False
            sage: A + B == J0(33).old_subvariety()
            True

        ::

             sage: J = J0(37) ; A, B = J.decomposition()
            sage: J.projection(A)._image_of_abvar(A)
            Abelian subvariety of dimension 1 of J0(37)
            sage: J.projection(A)._image_of_abvar(B)
            Abelian subvariety of dimension 0 of J0(37)
            sage: J.projection(B)._image_of_abvar(A)
            Abelian subvariety of dimension 0 of J0(37)
            sage: J.projection(B)._image_of_abvar(B)
            Abelian subvariety of dimension 1 of J0(37)
            sage: J.projection(B)._image_of_abvar(J)
            Abelian subvariety of dimension 1 of J0(37)
        """
        from .abvar import ModularAbelianVariety
        D = self.domain()
        C = self.codomain()
        if A is D:
            B = self.matrix()
        else:
            if not A.is_subvariety(D):
                raise ValueError("A must be an abelian subvariety of self.")
            # Write the vector space corresponding to A in terms of self's
            # vector space, then take the image under self.
            B = D.vector_space().coordinate_module(A.vector_space()).basis_matrix() * self.matrix()

        V = (B * C.vector_space().basis_matrix()).row_module(QQ)

        lattice = V.intersection(C.lattice())
        base_field = C.base_field()
        return ModularAbelianVariety(C.groups(), lattice, base_field)


class Morphism(Morphism_abstract, sage.modules.matrix_morphism.MatrixMorphism):

    def restrict_domain(self, sub):
        """
        Restrict self to the subvariety sub of self.domain().

        EXAMPLES::

            sage: J = J0(37) ; A, B = J.decomposition()
            sage: A.lattice().matrix()
            [ 1 -1  1  0]
            [ 0  0  2 -1]
            sage: B.lattice().matrix()
            [1 1 1 0]
            [0 0 0 1]
            sage: T = J.hecke_operator(2) ; T.matrix()
            [-1  1  1 -1]
            [ 1 -1  1  0]
            [ 0  0 -2  1]
            [ 0  0  0  0]
            sage: T.restrict_domain(A)
            Abelian variety morphism:
              From: Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37)
              To:   Abelian variety J0(37) of dimension 2
            sage: T.restrict_domain(A).matrix()
            [-2  2 -2  0]
            [ 0  0 -4  2]
            sage: T.restrict_domain(B)
            Abelian variety morphism:
              From: Simple abelian subvariety 37b(1,37) of dimension 1 of J0(37)
              To:   Abelian variety J0(37) of dimension 2
            sage: T.restrict_domain(B).matrix()
            [0 0 0 0]
            [0 0 0 0]
        """
        if not sub.is_subvariety(self.domain()):
            raise ValueError("sub must be a subvariety of self.domain()")

        if sub == self.domain():
            return self

        L = self.domain().lattice()
        B = sub.lattice().basis()
        ims = sum([ (L(b)*self.matrix()).list() for b in B], [])
        MS = matrix_space.MatrixSpace(self.base_ring(), len(B), self.codomain().rank())
        H = sub.Hom(self.codomain(), self.category_for())
        return H(MS(ims))

class DegeneracyMap(Morphism):
    def __init__(self, parent, A, t, side="left"):
        """
        Create the degeneracy map of index t in parent defined by the
        matrix A.

        INPUT:


        -  ``parent`` - a space of homomorphisms of abelian
           varieties

        -  ``A`` - a matrix defining self

        -  ``t`` - a list of indices defining the degeneracy
           map


        EXAMPLES::

            sage: J0(44).degeneracy_map(11,2)
            Degeneracy map from Abelian variety J0(44) of dimension 4 to Abelian variety J0(11) of dimension 1 defined by [2]
            sage: J0(44)[0].degeneracy_map(88,2)
            Degeneracy map from Simple abelian subvariety 11a(1,44) of dimension 1 of J0(44) to Abelian variety J0(88) of dimension 9 defined by [2]
        """
        if not isinstance(t, list):
            t = [t]
        self._t = t
        Morphism.__init__(self, parent, A, side)

    def t(self):
        """
        Return the list of indices defining self.

        EXAMPLES::

            sage: J0(22).degeneracy_map(44).t()
            [1]
            sage: J = J0(22) * J0(11)
            sage: J.degeneracy_map([44,44], [2,1])
            Degeneracy map from Abelian variety J0(22) x J0(11) of dimension 3 to Abelian variety J0(44) x J0(44) of dimension 8 defined by [2, 1]
            sage: J.degeneracy_map([44,44], [2,1]).t()
            [2, 1]
        """
        return self._t

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: J0(22).degeneracy_map(44)._repr_()
            'Degeneracy map from Abelian variety J0(22) of dimension 2 to Abelian variety J0(44) of dimension 4 defined by [1]'
        """
        return "Degeneracy map from %s to %s defined by %s"%(self.domain(), self.codomain(), self._t)

class HeckeOperator(Morphism):
    """
    A Hecke operator acting on a modular abelian variety.
    """
    def __init__(self, abvar, n, side="left"):
        """
        Create the Hecke operator of index `n` acting on the
        abelian variety abvar.

        INPUT:


        -  ``abvar`` - a modular abelian variety

        -  ``n`` - a positive integer


        EXAMPLES::

            sage: J = J0(37)
            sage: T2 = J.hecke_operator(2); T2
            Hecke operator T_2 on Abelian variety J0(37) of dimension 2
            sage: T2.parent()
            Endomorphism ring of Abelian variety J0(37) of dimension 2
        """
        from .abvar import is_ModularAbelianVariety
        n = ZZ(n)
        if n <= 0:
            raise ValueError("n must be positive")
        if not is_ModularAbelianVariety(abvar):
            raise TypeError("abvar must be a modular abelian variety")
        self.__abvar = abvar
        self.__n = n
        sage.modules.matrix_morphism.MatrixMorphism_abstract.__init__(self, abvar.Hom(abvar), side)

    def _repr_(self):
        """
        String representation of this Hecke operator.

        EXAMPLES::

            sage: J = J0(37)
            sage: J.hecke_operator(2)._repr_()
            'Hecke operator T_2 on Abelian variety J0(37) of dimension 2'
        """
        return "Hecke operator T_%s on %s"%(self.__n, self.__abvar)

    def index(self):
        """
        Return the index of this Hecke operator. (For example, if this is
        the operator `T_n`, then the index is the integer
        `n`.)

        OUTPUT:


        -  ``n`` - a (Sage) Integer


        EXAMPLES::

            sage: J = J0(15)
            sage: t = J.hecke_operator(53)
            sage: t
            Hecke operator T_53 on Abelian variety J0(15) of dimension 1
            sage: t.index()
            53
            sage: t = J.hecke_operator(54)
            sage: t
            Hecke operator T_54 on Abelian variety J0(15) of dimension 1
            sage: t.index()
            54

        ::

            sage: J = J1(12345)
            sage: t = J.hecke_operator(997) ; t
            Hecke operator T_997 on Abelian variety J1(12345) of dimension 5405473
            sage: t.index()
            997
            sage: type(t.index())
            <class 'sage.rings.integer.Integer'>
        """
        return self.__n

    def n(self):
        r"""
        Alias for ``self.index()``.

        EXAMPLES::

            sage: J = J0(17)
            sage: J.hecke_operator(5).n()
            5
        """
        return self.index()

    def characteristic_polynomial(self, var='x'):
        """
        Return the characteristic polynomial of this Hecke operator in the
        given variable.

        INPUT:


        -  ``var`` - a string (default: 'x')


        OUTPUT: a polynomial in var over the rational numbers.

        EXAMPLES::

            sage: A = J0(43)[1]; A
            Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: f = t2.characteristic_polynomial(); f
            x^2 - 2
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: f.factor()
            x^2 - 2
            sage: t2.characteristic_polynomial('y')
            y^2 - 2
        """
        return self.__abvar.rational_homology().hecke_polynomial(self.__n, var).change_ring(ZZ)

    def charpoly(self, var='x'):
        r"""
        Synonym for ``self.characteristic_polynomial(var)``.

        INPUT:


        -  ``var`` - string (default: 'x')


        EXAMPLES::

            sage: A = J1(13)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Abelian variety J1(13) of dimension 2
            sage: f = t2.charpoly(); f
            x^2 + 3*x + 3
            sage: f.factor()
            x^2 + 3*x + 3
            sage: t2.charpoly('y')
            y^2 + 3*y + 3
        """
        return self.characteristic_polynomial(var)

    def action_on_homology(self, R=ZZ):
        r"""
        Return the action of this Hecke operator on the homology
        `H_1(A; R)` of this abelian variety with coefficients in
        `R`.

        EXAMPLES::

            sage: A = J0(43)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Abelian variety J0(43) of dimension 3
            sage: h2 = t2.action_on_homology(); h2
            Hecke operator T_2 on Integral Homology of Abelian variety J0(43) of dimension 3
            sage: h2.matrix()
            [-2  1  0  0  0  0]
            [-1  1  1  0 -1  0]
            [-1  0 -1  2 -1  1]
            [-1  0  1  1 -1  1]
            [ 0 -2  0  2 -2  1]
            [ 0 -1  0  1  0 -1]
            sage: h2 = t2.action_on_homology(GF(2)); h2
            Hecke operator T_2 on Homology with coefficients in Finite Field of size 2 of Abelian variety J0(43) of dimension 3
            sage: h2.matrix()
            [0 1 0 0 0 0]
            [1 1 1 0 1 0]
            [1 0 1 0 1 1]
            [1 0 1 1 1 1]
            [0 0 0 0 0 1]
            [0 1 0 1 0 1]
        """
        return self.__abvar.homology(R).hecke_operator(self.index())

    def matrix(self):
        r"""
        Return the matrix of self acting on the homology
        `H_1(A, ZZ)` of this abelian variety with coefficients in
        `\ZZ`.

        EXAMPLES::

            sage: J0(47).hecke_operator(3).matrix()
            [ 0  0  1 -2  1  0 -1  0]
            [ 0  0  1  0 -1  0  0  0]
            [-1  2  0  0  2 -2  1 -1]
            [-2  1  1 -1  3 -1 -1  0]
            [-1 -1  1  0  1  0 -1  1]
            [-1  0  0 -1  2  0 -1  0]
            [-1 -1  2 -2  2  0 -1  0]
            [ 0 -1  0  0  1  0 -1  1]

        ::

            sage: J0(11).hecke_operator(7).matrix()
            [-2  0]
            [ 0 -2]
            sage: (J0(11) * J0(33)).hecke_operator(7).matrix()
            [-2  0  0  0  0  0  0  0]
            [ 0 -2  0  0  0  0  0  0]
            [ 0  0  0  0  2 -2  2 -2]
            [ 0  0  0 -2  2  0  2 -2]
            [ 0  0  0  0  2  0  4 -4]
            [ 0  0 -4  0  2  2  2 -2]
            [ 0  0 -2  0  2  2  0 -2]
            [ 0  0 -2  0  0  2  0 -2]

        ::

            sage: J0(23).hecke_operator(2).matrix()
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
        """
        try:
            return self._matrix
        except AttributeError:
            pass
        self._matrix = self.action_on_homology().matrix()
        return self._matrix
