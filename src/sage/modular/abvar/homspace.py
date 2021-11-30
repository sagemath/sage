"""
Spaces of homomorphisms between modular abelian varieties

EXAMPLES:

First, we consider J0(37). This Jacobian has two simple factors,
corresponding to distinct newforms. These two intersect
nontrivially in J0(37).

::

    sage: J = J0(37)
    sage: D = J.decomposition() ; D
    [
    Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37),
    Simple abelian subvariety 37b(1,37) of dimension 1 of J0(37)
    ]
    sage: D[0].intersection(D[1])
    (Finite subgroup with invariants [2, 2] over QQ of Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37),
     Simple abelian subvariety of dimension 0 of J0(37))

As an abstract product, since these newforms are distinct, the
corresponding simple abelian varieties are not isogenous, and so
there are no maps between them. The endomorphism ring of the
corresponding product is thus isomorphic to the direct sum of the
endomorphism rings for each factor. Since the factors correspond to
abelian varieties of dimension 1, these endomorphism rings are each
isomorphic to ZZ.

::

    sage: Hom(D[0],D[1]).gens()
    ()
    sage: A = D[0] * D[1] ; A
    Abelian subvariety of dimension 2 of J0(37) x J0(37)
    sage: A.endomorphism_ring().gens()
    (Abelian variety endomorphism of Abelian subvariety of dimension 2 of J0(37) x J0(37),
     Abelian variety endomorphism of Abelian subvariety of dimension 2 of J0(37) x J0(37))
    sage: [ x.matrix() for x in A.endomorphism_ring().gens() ]
    [
    [1 0 0 0]  [0 0 0 0]
    [0 1 0 0]  [0 0 0 0]
    [0 0 0 0]  [0 0 1 0]
    [0 0 0 0], [0 0 0 1]
    ]

However, these two newforms have a congruence between them modulo
2, which gives rise to interesting endomorphisms of J0(37).

::

    sage: E = J.endomorphism_ring()
    sage: E.gens()
    (Abelian variety endomorphism of Abelian variety J0(37) of dimension 2,
     Abelian variety endomorphism of Abelian variety J0(37) of dimension 2)
    sage: [ x.matrix() for x in E.gens() ]
    [
    [1 0 0 0]  [ 0  1  1 -1]
    [0 1 0 0]  [ 1  0  1  0]
    [0 0 1 0]  [ 0  0 -1  1]
    [0 0 0 1], [ 0  0  0  1]
    ]
    sage: (-1*E.gens()[0] + E.gens()[1]).matrix()
    [-1  1  1 -1]
    [ 1 -1  1  0]
    [ 0  0 -2  1]
    [ 0  0  0  0]

Of course, these endomorphisms will be reflected in the Hecke
algebra, which is in fact the full endomorphism ring of J0(37) in
this case::

    sage: J.hecke_operator(2).matrix()
    [-1  1  1 -1]
    [ 1 -1  1  0]
    [ 0  0 -2  1]
    [ 0  0  0  0]
    sage: T = E.image_of_hecke_algebra()
    sage: T.gens()
    (Abelian variety endomorphism of Abelian variety J0(37) of dimension 2,
     Abelian variety endomorphism of Abelian variety J0(37) of dimension 2)
    sage: [ x.matrix() for x in T.gens() ]
    [
    [1 0 0 0]  [ 0  1  1 -1]
    [0 1 0 0]  [ 1  0  1  0]
    [0 0 1 0]  [ 0  0 -1  1]
    [0 0 0 1], [ 0  0  0  1]
    ]
    sage: T.index_in(E)
    1

Next, we consider J0(33). In this case, we have both oldforms and
newforms. There are two copies of J0(11), one for each degeneracy
map from J0(11) to J0(33). There is also one newform at level 33.
The images of the two degeneracy maps are, of course, isogenous.

::

    sage: J = J0(33)
    sage: D = J.decomposition()
    sage: D
    [
    Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
    Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
    Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
    ]
    sage: Hom(D[0],D[1]).gens()
    (Abelian variety morphism:
      From: Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
      To:   Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),)
    sage: Hom(D[0],D[1]).gens()[0].matrix()
    [ 0  1]
    [-1  0]

Then this gives that the component corresponding to the sum of the
oldforms will have a rank 4 endomorphism ring. We also have a rank
one endomorphism ring for the newform 33a (since it is again
1-dimensional), which gives a rank 5 endomorphism ring for J0(33).

::

    sage: DD = J.decomposition(simple=False) ; DD
    [
    Abelian subvariety of dimension 2 of J0(33),
    Abelian subvariety of dimension 1 of J0(33)
    ]
    sage: A, B = DD
    sage: A == D[0] + D[1]
    True
    sage: A.endomorphism_ring().gens()
    (Abelian variety endomorphism of Abelian subvariety of dimension 2 of J0(33),
     Abelian variety endomorphism of Abelian subvariety of dimension 2 of J0(33),
     Abelian variety endomorphism of Abelian subvariety of dimension 2 of J0(33),
     Abelian variety endomorphism of Abelian subvariety of dimension 2 of J0(33))
    sage: B.endomorphism_ring().gens()
    (Abelian variety endomorphism of Abelian subvariety of dimension 1 of J0(33),)
    sage: E = J.endomorphism_ring() ; E.gens()  # long time (3s on sage.math, 2011)
    (Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
     Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
     Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
     Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
     Abelian variety endomorphism of Abelian variety J0(33) of dimension 3)

In this case, the image of the Hecke algebra will only have rank 3,
so that it is of infinite index in the full endomorphism ring.
However, if we call this image T, we can still ask about the index
of T in its saturation, which is 1 in this case.

::

    sage: T = E.image_of_hecke_algebra()  # long time
    sage: T.gens()  # long time
    (Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
     Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
     Abelian variety endomorphism of Abelian variety J0(33) of dimension 3)
    sage: T.index_in(E)  # long time
    +Infinity
    sage: T.index_in_saturation()  # long time
    1

AUTHORS:

- William Stein (2007-03)

- Craig Citro, Robert Bradshaw (2008-03): Rewrote with modabvar overhaul
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from copy import copy

from sage.categories.homset import HomsetWithBase
from sage.structure.all import parent
from sage.misc.lazy_attribute import lazy_attribute


from . import morphism

import sage.rings.integer_ring
import sage.rings.all

from sage.rings.ring import Ring
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix, identity_matrix
from sage.structure.element import is_Matrix

ZZ = sage.rings.integer_ring.ZZ

class Homspace(HomsetWithBase):
    """
    A space of homomorphisms between two modular abelian varieties.
    """
    Element = morphism.Morphism
    def __init__(self, domain, codomain, cat):
        """
        Create a homspace.

        INPUT:


        -  ``domain, codomain`` - modular abelian varieties

        -  ``cat`` - category


        EXAMPLES::

            sage: H = Hom(J0(11), J0(22)); H
            Space of homomorphisms from Abelian variety J0(11) of dimension 1 to Abelian variety J0(22) of dimension 2
            sage: Hom(J0(11), J0(11))
            Endomorphism ring of Abelian variety J0(11) of dimension 1
            sage: type(H)
            <class 'sage.modular.abvar.homspace.Homspace_with_category'>
            sage: H.homset_category()
            Category of modular abelian varieties over Rational Field
        """
        from .abvar import is_ModularAbelianVariety
        if not is_ModularAbelianVariety(domain):
            raise TypeError("domain must be a modular abelian variety")
        if not is_ModularAbelianVariety(codomain):
            raise TypeError("codomain must be a modular abelian variety")
        self._gens = None
        HomsetWithBase.__init__(self, domain, codomain, category=cat)

    def identity(self):
        """
        Return the identity endomorphism.

        EXAMPLES::

            sage: E = End(J0(11))
            sage: E.identity()
            Abelian variety endomorphism of Abelian variety J0(11) of dimension 1
            sage: E.one()
            Abelian variety endomorphism of Abelian variety J0(11) of dimension 1

            sage: H = Hom(J0(11), J0(22))
            sage: H.identity()
            Traceback (most recent call last):
            ...
            TypeError: the identity map is only defined for endomorphisms
        """
        if self.domain() is not self.codomain():
            raise TypeError("the identity map is only defined for endomorphisms")
        M = self.matrix_space().one()
        return self.element_class(self, M)

    @lazy_attribute
    def _matrix_space(self):
        """
        Return the matrix space of ``self``.

        .. WARNING::

            During unpickling, the domain and codomain may be unable to
            provide the necessary information. This is why this is a lazy
            attribute. See :trac:`14793`.

        EXAMPLES::

            sage: Hom(J0(11), J0(22))._matrix_space
            Full MatrixSpace of 2 by 4 dense matrices over Integer Ring
        """
        return MatrixSpace(ZZ,2*self.domain().dimension(), 2*self.codomain().dimension())

    def _element_constructor_from_element_class(self, *args, **keywords):
        """
        Used in the coercion framework. Unfortunately, the default method
        would get the order of parent and data different from what is expected
        in ``MatrixMorphism.__init__``.

        EXAMPLES::

            sage: H = Hom(J0(11), J0(22))
            sage: phi = H(matrix(ZZ,2,4,[5..12])); phi # indirect doctest
            Abelian variety morphism:
              From: Abelian variety J0(11) of dimension 1
              To:   Abelian variety J0(22) of dimension 2
        """
        return self.element_class(self, *args, **keywords)

    def __call__(self, M, **kwds):
        r"""
        Create a homomorphism in this space from M. M can be any of the
        following:

        - a Morphism of abelian varieties

        - a matrix of the appropriate size
          (i.e. 2\*self.domain().dimension() x
          2\*self.codomain().dimension()) whose entries are coercible
          into self.base_ring()

        - anything that can be coerced into self.matrix_space()

        EXAMPLES::

            sage: H = Hom(J0(11), J0(22))
            sage: phi = H(matrix(ZZ,2,4,[5..12])) ; phi
            Abelian variety morphism:
              From: Abelian variety J0(11) of dimension 1
              To:   Abelian variety J0(22) of dimension 2
            sage: phi.matrix()
            [ 5  6  7  8]
            [ 9 10 11 12]
            sage: phi.matrix().parent()
            Full MatrixSpace of 2 by 4 dense matrices over Integer Ring

        ::

            sage: H = J0(22).Hom(J0(11)*J0(11))
            sage: m1 = J0(22).degeneracy_map(11,1).matrix() ; m1
            [ 0  1]
            [-1  1]
            [-1  0]
            [ 0 -1]
            sage: m2 = J0(22).degeneracy_map(11,2).matrix() ; m2
            [ 1 -2]
            [ 0 -2]
            [ 1 -1]
            [ 0 -1]
            sage: m = m1.transpose().stack(m2.transpose()).transpose() ; m
            [ 0  1  1 -2]
            [-1  1  0 -2]
            [-1  0  1 -1]
            [ 0 -1  0 -1]
            sage: phi = H(m) ; phi
            Abelian variety morphism:
              From: Abelian variety J0(22) of dimension 2
              To:   Abelian variety J0(11) x J0(11) of dimension 2
            sage: phi.matrix()
            [ 0  1  1 -2]
            [-1  1  0 -2]
            [-1  0  1 -1]
            [ 0 -1  0 -1]
        """
        side = kwds.get("side", "left")
        if isinstance(M, morphism.Morphism):
            if M.parent() is self:
                return M
            elif M.domain() == self.domain() and M.codomain() == self.codomain():
                M = M.matrix()
            else:
                raise ValueError("cannot convert %s into %s" % (M, self)) 
        elif is_Matrix(M):
            if M.base_ring() != ZZ:
                M = M.change_ring(ZZ)
            if side == "left":
                if M.nrows() != 2*self.domain().dimension() or M.ncols() != 2*self.codomain().dimension():
                    raise TypeError("matrix has wrong dimension")
            else:
                if M.ncols() != 2*self.domain().dimension() or M.nrows() != 2*self.codomain().dimension():
                    raise TypeError("matrix has wrong dimension")
        elif self.matrix_space().has_coerce_map_from(parent(M)):
            M = self.matrix_space()(M)
        else:
            raise TypeError("can only coerce in matrices or morphisms")
        return self.element_class(self, M, side)

    def _coerce_impl(self, x):
        """
        Coerce x into self, if possible.

        EXAMPLES::

            sage: J = J0(37) ; J.Hom(J)._coerce_impl(matrix(ZZ,4,[5..20]))
            Abelian variety endomorphism of Abelian variety J0(37) of dimension 2
            sage: K = J0(11) * J0(11) ; J.Hom(K)._coerce_impl(matrix(ZZ,4,[5..20]))
            Abelian variety morphism:
              From: Abelian variety J0(37) of dimension 2
              To:   Abelian variety J0(11) x J0(11) of dimension 2
        """
        if self.matrix_space().has_coerce_map_from(parent(x)):
            return self(x)
        else:
            return HomsetWithBase._coerce_impl(self, x)

    def _repr_(self):
        """
        String representation of a modular abelian variety homspace.

        EXAMPLES::

            sage: J = J0(11)
            sage: End(J)._repr_()
            'Endomorphism ring of Abelian variety J0(11) of dimension 1'
        """
        return "Space of homomorphisms from %s to %s"%\
               (self.domain(), self.codomain())

    def _get_matrix(self, g):
        """
        Given an object g, try to return a matrix corresponding to g with
        dimensions the same as those of self.matrix_space().

        INPUT:


        -  ``g`` - a matrix or morphism or object with a list
           method


        OUTPUT: a matrix

        EXAMPLES::

            sage: E = End(J0(11))
            sage: E._get_matrix(matrix(QQ,2,[1,2,3,4]))
            [1 2]
            [3 4]
            sage: E._get_matrix(J0(11).hecke_operator(2))
            [-2  0]
            [ 0 -2]

        ::

            sage: H = Hom(J0(11) * J0(17), J0(22))
            sage: H._get_matrix(tuple([8..23]))
            [ 8  9 10 11]
            [12 13 14 15]
            [16 17 18 19]
            [20 21 22 23]
            sage: H._get_matrix(tuple([8..23]))
            [ 8  9 10 11]
            [12 13 14 15]
            [16 17 18 19]
            [20 21 22 23]
            sage: H._get_matrix([8..23])
            [ 8  9 10 11]
            [12 13 14 15]
            [16 17 18 19]
            [20 21 22 23]
        """
        try:
            if g.parent() is self.matrix_space():
                return g
        except AttributeError:
            pass

        if isinstance(g, morphism.Morphism):
            return g.matrix()
        elif hasattr(g, 'list'):
            return self.matrix_space()(g.list())
        else:
            return self.matrix_space()(g)

    def free_module(self):
        r"""
        Return this endomorphism ring as a free submodule of a big
        `\ZZ^{4nm}`, where `n` is the dimension of
        the domain abelian variety and `m` the dimension of the
        codomain.

        OUTPUT: free module

        EXAMPLES::

            sage: E = Hom(J0(11), J0(22))
            sage: E.free_module()
            Free module of degree 8 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  1  1  1 -1 -1]
            [ 0  1 -3  1  1  1 -1  0]
        """
        self.calculate_generators()
        V = ZZ**(4*self.domain().dimension() * self.codomain().dimension())
        return V.submodule([ V(m.matrix().list()) for m in self.gens() ])

    def gen(self, i=0):
        """
        Return i-th generator of self.

        INPUT:


        -  ``i`` - an integer


        OUTPUT: a morphism

        EXAMPLES::

            sage: E = End(J0(22))
            sage: E.gen(0).matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        self.calculate_generators()
        if i > self.ngens():
            raise ValueError("self only has %s generators"%self.ngens())
        return self.element_class(self, self._gens[i])

    def ngens(self):
        """
        Return number of generators of self.

        OUTPUT: integer

        EXAMPLES::

            sage: E = End(J0(22))
            sage: E.ngens()
            4
        """
        self.calculate_generators()
        return len(self._gens)

    def gens(self):
        """
        Return tuple of generators for this endomorphism ring.

        EXAMPLES::

            sage: E = End(J0(22))
            sage: E.gens()
            (Abelian variety endomorphism of Abelian variety J0(22) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(22) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(22) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(22) of dimension 2)
        """
        try:
            return self._gen_morphisms
        except AttributeError:
            self.calculate_generators()
            self._gen_morphisms = tuple([self.gen(i) for i in range(self.ngens())])
            return self._gen_morphisms

    def matrix_space(self):
        """
        Return the underlying matrix space that we view this endomorphism
        ring as being embedded into.

        EXAMPLES::

            sage: E = End(J0(22))
            sage: E.matrix_space()
            Full MatrixSpace of 4 by 4 dense matrices over Integer Ring
        """
        return self._matrix_space

    def calculate_generators(self):
        """
        If generators haven't already been computed, calculate generators
        for this homspace. If they have been computed, do nothing.

        EXAMPLES::

            sage: E = End(J0(11))
            sage: E.calculate_generators()
        """

        if self._gens is not None:
            return

        if (self.domain() == self.codomain()) and (self.domain().dimension() == 1):
            self._gens = tuple([ identity_matrix(ZZ,2) ])
            return

        phi = self.domain()._isogeny_to_product_of_powers()
        psi = self.codomain()._isogeny_to_product_of_powers()

        H_simple = phi.codomain().Hom(psi.codomain())
        im_gens = H_simple._calculate_product_gens()

        M = phi.matrix()
        Mt = psi.complementary_isogeny().matrix()

        R = ZZ**(4*self.domain().dimension()*self.codomain().dimension())
        gens = R.submodule([ (M*self._get_matrix(g)*Mt).list()
                             for g in im_gens ]).saturation().basis()
        self._gens = tuple([ self._get_matrix(g) for g in gens ])

    def _calculate_product_gens(self):
        """
        For internal use.

        Calculate generators for self, assuming that self is a product of
        simple factors.

        EXAMPLES::

            sage: E = End(J0(37))
            sage: E.gens()
            (Abelian variety endomorphism of Abelian variety J0(37) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(37) of dimension 2)
            sage: [ x.matrix() for x in E.gens() ]
            [
            [1 0 0 0]  [ 0  1  1 -1]
            [0 1 0 0]  [ 1  0  1  0]
            [0 0 1 0]  [ 0  0 -1  1]
            [0 0 0 1], [ 0  0  0  1]
            ]
            sage: E._calculate_product_gens()
            [
            [1 0 0 0]  [0 0 0 0]
            [0 1 0 0]  [0 0 0 0]
            [0 0 0 0]  [0 0 1 0]
            [0 0 0 0], [0 0 0 1]
            ]
        """
        Afactors = self.domain().decomposition(simple=False)
        Bfactors = self.codomain().decomposition(simple=False)
        if len(Afactors) == 1 and len(Bfactors) == 1:
            Asimples = Afactors[0].decomposition()
            Bsimples = Bfactors[0].decomposition()
            if len(Asimples) == 1 and len(Bsimples) == 1:
                # Handle the base case of A, B simple
                gens = self._calculate_simple_gens()

            else:
                # Handle the case of A, B simple powers
                gens = []
                phi_matrix = Afactors[0]._isogeny_to_product_of_simples().matrix()
                psi_t_matrix = Bfactors[0]._isogeny_to_product_of_simples().complementary_isogeny().matrix()
                for i in range(len(Asimples)):
                    for j in range(len(Bsimples)):
                        hom_gens = Asimples[i].Hom(Bsimples[j]).gens()
                        for sub_gen in hom_gens:
                            sub_mat = sub_gen.matrix()
                            M = copy(self.matrix_space().zero_matrix())
                            M.set_block(sub_mat.nrows()*i, sub_mat.ncols()*j, sub_mat)
                            gens.append(phi_matrix * M * psi_t_matrix)

        else:
            # Handle the case of A, B generic
            gens = []
            cur_row = 0
            for Afactor in Afactors:
                cur_row += Afactor.dimension() * 2
                cur_col = 0
                for Bfactor in Bfactors:
                    cur_col += Bfactor.dimension() * 2
                    Asimple = Afactor[0]
                    Bsimple = Bfactor[0]
                    if Asimple.newform_label() == Bsimple.newform_label():
                        for sub_gen in Afactor.Hom(Bfactor).gens():
                            sub_mat = sub_gen.matrix()
                            M = copy(self.matrix_space().zero_matrix())
                            M.set_block(cur_row - sub_mat.nrows(),
                                        cur_col - sub_mat.ncols(),
                                        sub_mat)
                            gens.append(M)

        return gens

    def _calculate_simple_gens(self):
        """
        Calculate generators for self, where both the domain and codomain
        for self are assumed to be simple abelian varieties. The saturation
        of the span of these generators in self will be the full space of
        homomorphisms from the domain of self to its codomain.

        EXAMPLES::

            sage: H = Hom(J0(11), J0(22)[0])
            sage: H._calculate_simple_gens()
            [
            [1 0]
            [1 1]
            ]
            sage: J = J0(11) * J0(33) ; J.decomposition()
            [
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J0(33),
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(11) x J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(11) x J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(11) x J0(33)
            ]
            sage: J[0].Hom(J[1])._calculate_simple_gens()
            [
            [ 0 -1]
            [ 1 -1]
            ]
            sage: J[0].Hom(J[2])._calculate_simple_gens()
            [
            [-1  0]
            [-1 -1]
            ]
            sage: J[0].Hom(J[0])._calculate_simple_gens()
            [
            [1 0]
            [0 1]
            ]
            sage: J[1].Hom(J[2])._calculate_simple_gens()
            [
            [ 0 -4]
            [ 4  0]
            ]

        ::

            sage: J = J0(23) ; J.decomposition()
            [
            Simple abelian variety J0(23) of dimension 2
            ]
            sage: J[0].Hom(J[0])._calculate_simple_gens()
            [
            [1 0 0 0]  [ 0  1 -1  0]
            [0 1 0 0]  [ 0  1 -1  1]
            [0 0 1 0]  [-1  2 -2  1]
            [0 0 0 1], [-1  1  0 -1]
            ]
            sage: J.hecke_operator(2).matrix()
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]

        ::

            sage: H = Hom(J0(11), J0(22)[0])
            sage: H._calculate_simple_gens()
            [
            [1 0]
            [1 1]
            ]
        """
        A = self.domain()
        B = self.codomain()

        if A.newform_label() != B.newform_label():
            return []

        f = A._isogeny_to_newform_abelian_variety()
        g = B._isogeny_to_newform_abelian_variety().complementary_isogeny()

        Af = f.codomain()
        ls = Af._calculate_endomorphism_generators()

        Mf = f.matrix()
        Mg = g.matrix()

        return [ Mf * self._get_matrix(e) * Mg for e in ls ]

# NOTE/WARNING/TODO:  Below in the __init__, etc. we do *not* check
# that the input gens are give something that spans a sub*ring*, as apposed
# to just a subgroup.
class EndomorphismSubring(Homspace, Ring):

    def __init__(self, A, gens=None, category=None):
        """
        A subring of the endomorphism ring.

        INPUT:


        -  ``A`` - an abelian variety

        -  ``gens`` - (default: None); optional; if given
           should be a tuple of the generators as matrices


        EXAMPLES::

            sage: J0(23).endomorphism_ring()
            Endomorphism ring of Abelian variety J0(23) of dimension 2
            sage: sage.modular.abvar.homspace.EndomorphismSubring(J0(25))
            Endomorphism ring of Abelian variety J0(25) of dimension 0
            sage: E = J0(11).endomorphism_ring()
            sage: type(E)
            <class 'sage.modular.abvar.homspace.EndomorphismSubring_with_category'>
            sage: E.homset_category()
            Category of modular abelian varieties over Rational Field
            sage: E.category()
            Category of endsets of modular abelian varieties over Rational Field
            sage: E in Rings()
            True
            sage: TestSuite(E).run(skip=["_test_prod"])

        TESTS:

        The following tests against a problem on 32 bit machines that
        occurred while working on :trac:`9944`::

            sage: sage.modular.abvar.homspace.EndomorphismSubring(J1(12345))
            Endomorphism ring of Abelian variety J1(12345) of dimension 5405473

        :trac:`16275` removed the custom ``__reduce__`` method, since
        :meth:`Homset.__reduce__` already implements appropriate
        unpickling by construction::

            sage: E.__reduce__.__module__
            'sage.categories.homset'
            sage: E.__reduce__()
            (<function Hom at ...>,
             (Abelian variety J0(11) of dimension 1,
              Abelian variety J0(11) of dimension 1,
              Category of modular abelian varieties over Rational Field,
             False))
        """
        self._J = A.ambient_variety()
        self._A = A

        # Initialise self with the correct category.
        # We need to initialise it as a ring first
        if category is None:
            homset_cat = A.category()
        else:
            homset_cat = category
        # Remark: Ring.__init__ will automatically form the join
        # of the category of rings and of homset_cat
        Ring.__init__(self, A.base_ring(), category=homset_cat.Endsets())
        Homspace.__init__(self, A, A, cat=homset_cat)
        if gens is None:
            self._gens = None
        else:
            self._gens = tuple([ self._get_matrix(g) for g in gens ])
        self._is_full_ring = gens is None

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: J0(31).endomorphism_ring()._repr_()
            'Endomorphism ring of Abelian variety J0(31) of dimension 2'
            sage: J0(31).endomorphism_ring().image_of_hecke_algebra()._repr_()
            'Subring of endomorphism ring of Abelian variety J0(31) of dimension 2'
        """
        if self._is_full_ring:
            return "Endomorphism ring of %s" % self._A
        else:
            return "Subring of endomorphism ring of %s" % self._A

    def abelian_variety(self):
        """
        Return the abelian variety that this endomorphism ring is attached
        to.

        EXAMPLES::

            sage: J0(11).endomorphism_ring().abelian_variety()
            Abelian variety J0(11) of dimension 1
        """
        return self._A

    def index_in(self, other, check=True):
        """
        Return the index of self in other.

        INPUT:


        -  ``other`` - another endomorphism subring of the
           same abelian variety

        -  ``check`` - bool (default: True); whether to do some
           type and other consistency checks


        EXAMPLES::

            sage: R = J0(33).endomorphism_ring()
            sage: R.index_in(R)
            1
            sage: J = J0(37) ; E = J.endomorphism_ring() ; T = E.image_of_hecke_algebra()
            sage: T.index_in(E)
            1
            sage: J = J0(22) ; E = J.endomorphism_ring() ; T = E.image_of_hecke_algebra()
            sage: T.index_in(E)
            +Infinity
        """
        if check:
            if not isinstance(other, EndomorphismSubring):
                raise ValueError("other must be a subring of an endomorphism ring of an abelian variety.")
            if not (self.abelian_variety() == other.abelian_variety()):
                raise ValueError("self and other must be endomorphisms of the same abelian variety")

        M = self.free_module()
        N = other.free_module()
        if M.rank() < N.rank():
            return sage.rings.all.Infinity
        return M.index_in(N)

    def index_in_saturation(self):
        """
        Given a Hecke algebra T, compute its index in its saturation.

        EXAMPLES::

            sage: End(J0(23)).image_of_hecke_algebra().index_in_saturation()
            1
            sage: End(J0(44)).image_of_hecke_algebra().index_in_saturation()
            2
        """
        A = self.abelian_variety()
        d = A.dimension()
        M = ZZ**(4*d**2)
        gens = [ x.matrix().list() for x in self.gens() ]
        R = M.submodule(gens)
        return R.index_in_saturation()

    def discriminant(self):
        """
        Return the discriminant of this ring, which is the discriminant of
        the trace pairing.

        .. note::

           One knows that for modular abelian varieties, the
           endomorphism ring should be isomorphic to an order in a
           number field. However, the discriminant returned by this
           function will be `2^n` ( `n =`
           self.dimension()) times the discriminant of that order,
           since the elements are represented as 2d x 2d
           matrices. Notice, for example, that the case of a one
           dimensional abelian variety, whose endomorphism ring must
           be ZZ, has discriminant 2, as in the example below.

        EXAMPLES::

            sage: J0(33).endomorphism_ring().discriminant()
            -64800
            sage: J0(46).endomorphism_ring().discriminant()  # long time (6s on sage.math, 2011)
            24200000000
            sage: J0(11).endomorphism_ring().discriminant()
            2
        """
        g = self.gens()
        M = Matrix(ZZ,len(g), [ (g[i]*g[j]).trace()
                                for i in range(len(g)) for j in range(len(g)) ])
        return M.determinant()

    def image_of_hecke_algebra(self, check_every=1):
        """
        Compute the image of the Hecke algebra inside this endomorphism
        subring.

        We simply calculate Hecke operators up to the Sturm bound, and look
        at the submodule spanned by them. While computing, we can check to
        see if the submodule spanned so far is saturated and of maximal
        dimension, in which case we may be done. The optional argument
        check_every determines how many Hecke operators we add in before
        checking to see if this condition is met.

        INPUT:

        - ``check_every`` -- integer (default: 1) If this integer is positive,
          this integer determines how many Hecke operators we add in before
          checking to see if the submodule spanned so far is maximal and
          saturated.

        OUTPUT:

        - The image of the Hecke algebra as an subring of ``self``.

        EXAMPLES::

            sage: E = J0(33).endomorphism_ring()
            sage: E.image_of_hecke_algebra()
            Subring of endomorphism ring of Abelian variety J0(33) of dimension 3
            sage: E.image_of_hecke_algebra().gens()
            (Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
             Abelian variety endomorphism of Abelian variety J0(33) of dimension 3,
             Abelian variety endomorphism of Abelian variety J0(33) of dimension 3)
            sage: [ x.matrix() for x in E.image_of_hecke_algebra().gens() ]
            [
            [1 0 0 0 0 0]  [ 0  2  0 -1  1 -1]  [ 0  0  1 -1  1 -1]
            [0 1 0 0 0 0]  [-1 -2  2 -1  2 -1]  [ 0 -1  1  0  1 -1]
            [0 0 1 0 0 0]  [ 0  0  1 -1  3 -1]  [ 0  0  1  0  2 -2]
            [0 0 0 1 0 0]  [-2  2  0  1  1 -1]  [-2  0  1  1  1 -1]
            [0 0 0 0 1 0]  [-1  1  0  2  0 -3]  [-1  0  1  1  0 -1]
            [0 0 0 0 0 1], [-1  1 -1  1  1 -2], [-1  0  0  1  0 -1]
            ]
            sage: J0(33).hecke_operator(2).matrix()
            [-1  0  1 -1  1 -1]
            [ 0 -2  1  0  1 -1]
            [ 0  0  0  0  2 -2]
            [-2  0  1  0  1 -1]
            [-1  0  1  1 -1 -1]
            [-1  0  0  1  0 -2]
        """
        try:
            return self.__hecke_algebra_image
        except AttributeError:
            pass

        A = self.abelian_variety()
        if not A.is_hecke_stable():
            raise ValueError("ambient variety is not Hecke stable")

        M = A.modular_symbols()

        d = A.dimension()
        EndVecZ = ZZ**(4*d**2)

        if d == 1:
            self.__hecke_algebra_image = EndomorphismSubring(A, [[1,0,0,1]])
            return self.__hecke_algebra_image

        V = EndVecZ.submodule([A.hecke_operator(1).matrix().list()])

        for n in range(2,M.sturm_bound()+1):
            if (check_every > 0 and
                    n % check_every == 0 and
                    V.dimension() == d and
                    V.index_in_saturation() == 1):
                break
            V += EndVecZ.submodule([ A.hecke_operator(n).matrix().list() ])

        self.__hecke_algebra_image = EndomorphismSubring(A, V.basis())
        return self.__hecke_algebra_image



