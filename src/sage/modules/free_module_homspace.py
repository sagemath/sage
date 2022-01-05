r"""
Homspaces between free modules

EXAMPLES: We create `\mathrm{End}(\ZZ^2)` and compute a
basis.

::

    sage: M = FreeModule(IntegerRing(),2)
    sage: E = End(M)
    sage: B = E.basis()
    sage: len(B)
    4
    sage: B[0]
    Free module morphism defined by the matrix
    [1 0]
    [0 0]
    Domain: Ambient free module of rank 2 over the principal ideal domain ...
    Codomain: Ambient free module of rank 2 over the principal ideal domain ...

We create `\mathrm{Hom}(\ZZ^3, \ZZ^2)` and
compute a basis.

::

    sage: V3 = FreeModule(IntegerRing(),3)
    sage: V2 = FreeModule(IntegerRing(),2)
    sage: H = Hom(V3,V2)
    sage: H
    Set of Morphisms from Ambient free module of rank 3 over
     the principal ideal domain Integer Ring
     to Ambient free module of rank 2
     over the principal ideal domain Integer Ring
     in Category of finite dimensional modules with basis over
     (euclidean domains and infinite enumerated sets and metric spaces)
    sage: B = H.basis()
    sage: len(B)
    6
    sage: B[0]
    Free module morphism defined by the matrix
    [1 0]
    [0 0]
    [0 0]...

TESTS::

    sage: H = Hom(QQ^2, QQ^1)
    sage: loads(dumps(H)) == H
    True

See :trac:`5886`::

    sage: V = (ZZ^2).span_of_basis([[1,2],[3,4]])
    sage: V.hom([V.0, V.1])
    Free module morphism defined by the matrix
    [1 0]
    [0 1]...

"""

#*****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import sage.categories.homset
from sage.structure.element import is_Matrix
from sage.matrix.constructor import matrix, identity_matrix
from sage.matrix.matrix_space import MatrixSpace
from inspect import isfunction
from sage.misc.cachefunc import cached_method


def is_FreeModuleHomspace(x):
    r"""
    Return ``True`` if ``x`` is a free module homspace.

    EXAMPLES:

    Notice that every vector space is a free module, but when we construct
    a set of morphisms between two vector spaces, it is a
    ``VectorSpaceHomspace``, which qualifies as a ``FreeModuleHomspace``,
    since the former is special case of the latter.

        sage: H = Hom(ZZ^3, ZZ^2)
        sage: type(H)
        <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(H)
        True

        sage: K = Hom(QQ^3, ZZ^2)
        sage: type(K)
        <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(K)
        True

        sage: L = Hom(ZZ^3, QQ^2)
        sage: type(L)
        <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(L)
        True

        sage: P = Hom(QQ^3, QQ^2)
        sage: type(P)
        <class 'sage.modules.vector_space_homspace.VectorSpaceHomspace_with_category'>
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(P)
        True

        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace('junk')
        False
    """
    return isinstance(x, FreeModuleHomspace)

class FreeModuleHomspace(sage.categories.homset.HomsetWithBase):
    def __call__(self, A, check=True):
        r"""
        INPUT:

        - A -- either a matrix or a list/tuple of images of generators,
          or a function returning elements of the codomain for elements
          of the domain.
        - check -- bool (default: True)

        If A is a matrix, then it is the matrix of this linear
        transformation, with respect to the basis for the domain and
        codomain.  Thus the identity matrix always defines the
        identity morphism.

        EXAMPLES::

            sage: V = (ZZ^3).span_of_basis([[1,1,0],[1,0,2]])
            sage: H = V.Hom(V); H
            Set of Morphisms from ...
            sage: H([V.0,V.1])                    # indirect doctest
            Free module morphism defined by the matrix
            [1 0]
            [0 1]...
            sage: phi = H([V.1,V.0]); phi
            Free module morphism defined by the matrix
            [0 1]
            [1 0]...
            sage: phi(V.1) == V.0
            True
            sage: phi(V.0) == V.1
            True

        The following tests against a bug that was fixed in
        :trac:`9944`. The method ``zero()`` calls this hom space with
        a function, not with a matrix, and that case had previously
        not been taken care of::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ)
            sage: V.Hom(V).zero()   # indirect doctest
            Free module morphism defined by the matrix
            [0 0 0]
            [0 0 0]
            [0 0 0]
            Domain: Free module of degree 3 and rank 3 over Integer Ring
            Echelon ...
            Codomain: Free module of degree 3 and rank 3 over Integer Ring
            Echelon ...

        The following tests the bug fixed in :trac:`31818`. If there is no 
        coercion between base rings, one can only define the zero morphism, 
        as morphism of additive groups. Before one could for example use an 
        integer matrix to define a morphism from the rational numbers to the
        integers. ::

            sage: V = QQ^2; W = ZZ^2; m = identity_matrix(2)                                
            sage: H = V.Hom(W); H(m)                                                        
            Traceback (most recent call last):
            ...
            TypeError: nontrivial morphisms require a coercion map from the base ring of the domain to the base ring of the codomain
            sage: n = zero_matrix(2); 
            sage: h = H(n); h                                                                                                                                                                                           
            Free module morphism defined by the matrix
            [0 0]
            [0 0]
            Domain: Vector space of dimension 2 over Rational Field
            Codomain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: [h(v) for v in V.gens()]                                                                                                                                                                              
            [(0, 0), (0, 0)]
        """
        from . import free_module_morphism
        if not is_Matrix(A):
            # Compute the matrix of the morphism that sends the
            # generators of the domain to the elements of A.
            C = self.codomain()
            try:
                if isfunction(A):
                    v = [C(A(g)) for g in self.domain().gens()]
                else:
                    v = [C(a) for a in A]
                A = matrix([C.coordinates(a) for a in v], ncols=C.rank())
            except TypeError:
                # Let us hope that FreeModuleMorphism knows to handle
                # that case
                pass
        if not self.codomain().base_ring().has_coerce_map_from(self.domain().base_ring()) and not A.is_zero():
            raise TypeError("nontrivial morphisms require a coercion map from the base ring of the domain to the base ring of the codomain")
        return free_module_morphism.FreeModuleMorphism(self, A)

    @cached_method
    def zero(self):
        """
        EXAMPLES::

            sage: E = ZZ^2
            sage: F = ZZ^3
            sage: H = Hom(E, F)
            sage: f = H.zero()
            sage: f
            Free module morphism defined by the matrix
            [0 0 0]
            [0 0 0]
            Domain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: f(E.an_element())
            (0, 0, 0)
            sage: f(E.an_element()) == F.zero()
            True

        TESTS:

        We check that ``H.zero()`` is picklable::

            sage: loads(dumps(f.parent().zero()))
            Free module morphism defined by the matrix
            [0 0 0]
            [0 0 0]
            Domain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return self(lambda x: self.codomain().zero())

    def _matrix_space(self):
        """
        Return underlying matrix space that contains the matrices that define
        the homomorphisms in this free module homspace.

        OUTPUT:

        - matrix space

        EXAMPLES::

            sage: H = Hom(QQ^3, QQ^2)
            sage: H._matrix_space()
            Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        """
        try:
            return self.__matrix_space
        except AttributeError:
            R = self.codomain().base_ring()
            M = MatrixSpace(R, self.domain().rank(), self.codomain().rank())
            self.__matrix_space = M
            return M

    def basis(self):
        """
        Return a basis for this space of free module homomorphisms.

        OUTPUT:

        - tuple

        EXAMPLES::

            sage: H = Hom(ZZ^2, ZZ^1)
            sage: H.basis()
            (Free module morphism defined by the matrix
            [1]
            [0]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 1 over the principal ideal domain ..., Free module morphism defined by the matrix
            [0]
            [1]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 1 over the principal ideal domain ...)
        """
        try:
            return self.__basis
        except AttributeError:
            M = self._matrix_space()
            B = M.basis()
            self.__basis = tuple([self(x) for x in B])
            return self.__basis

    def identity(self):
        r"""
        Return identity morphism in an endomorphism ring.

        EXAMPLES::

            sage: V=FreeModule(ZZ,5)
            sage: H=V.Hom(V)
            sage: H.identity()
            Free module morphism defined by the matrix
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            Domain: Ambient free module of rank 5 over the principal ideal domain ...
            Codomain: Ambient free module of rank 5 over the principal ideal domain ...
        """
        if self.is_endomorphism_set():
            return self(identity_matrix(self.base_ring(),self.domain().rank()))
        else:
            raise TypeError("Identity map only defined for endomorphisms. Try natural_map() instead.")

