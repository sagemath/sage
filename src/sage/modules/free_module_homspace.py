# sage_setup: distribution = sagemath-modules
r"""
Homspaces between free modules

EXAMPLES:

We create `\mathrm{End}(\ZZ^2)` and compute a basis. ::

    sage: M = FreeModule(IntegerRing(), 2)
    sage: E = End(M)
    sage: B = E.basis()
    sage: len(B)
    4
    sage: B[0]
    Free module morphism defined by the matrix
    [1 0]
    [0 0]
    Domain:   Ambient free module of rank 2 over the principal ideal domain ...
    Codomain: Ambient free module of rank 2 over the principal ideal domain ...

We create `\mathrm{Hom}(\ZZ^3, \ZZ^2)` and compute a basis. ::

    sage: V3 = FreeModule(IntegerRing(), 3)
    sage: V2 = FreeModule(IntegerRing(), 2)
    sage: H = Hom(V3, V2)
    sage: H
    Set of Morphisms
     from Ambient free module of rank 3 over the principal ideal domain Integer Ring
       to Ambient free module of rank 2 over the principal ideal domain Integer Ring
       in Category of finite dimensional modules with basis over
          (Dedekind domains and euclidean domains
           and noetherian rings
           and infinite enumerated sets and metric spaces)
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

See :issue:`5886`::

    sage: V = (ZZ^2).span_of_basis([[1,2],[3,4]])
    sage: V.hom([V.0, V.1])
    Free module morphism defined by the matrix
    [1 0]
    [0 1]...

See :issue:`13321`::

    sage: (GF(7)^2).hom([[20, 0], [0, 21]], ZZ^2)
    Traceback (most recent call last):
    ...
    TypeError: nontrivial morphisms require a coercion map from the base ring
    of the domain to the base ring of the codomain
"""
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sage.categories.homset
from sage.structure.element import Matrix
from sage.matrix.constructor import matrix, identity_matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method


def is_FreeModuleHomspace(x):
    r"""
    Return ``True`` if ``x`` is a free module homspace.

    Notice that every vector space is a free module, but when we construct
    a set of morphisms between two vector spaces, it is a
    ``VectorSpaceHomspace``, which qualifies as a ``FreeModuleHomspace``,
    since the former is special case of the latter.

    EXAMPLES::

        sage: H = Hom(ZZ^3, ZZ^2)
        sage: type(H)
        <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(H)
        doctest:warning...
        DeprecationWarning: the function is_FreeModuleHomspace is deprecated;
        use 'isinstance(..., FreeModuleHomspace)' instead
        See https://github.com/sagemath/sage/issues/37924 for details.
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
    from sage.misc.superseded import deprecation
    deprecation(37924, "the function is_FreeModuleHomspace is deprecated; use 'isinstance(..., FreeModuleHomspace)' instead")
    return isinstance(x, FreeModuleHomspace)


class FreeModuleHomspace(sage.categories.homset.HomsetWithBase):
    def __call__(self, A, **kwds):
        r"""
        INPUT:

        - ``A`` -- either a matrix or a list/tuple of images of generators,
          or a function returning elements of the codomain for elements of the domain
        - ``check`` -- boolean (default: ``True``)
        - the keyword ``side`` can be assigned the values ``'left'`` or
          ``'right'``. It corresponds to the side of vectors relative to the matrix.

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
        :issue:`9944`. The method ``zero()`` calls this hom space with
        a function, not with a matrix, and that case had previously
        not been taken care of::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ)
            sage: V.Hom(V).zero()   # indirect doctest
            Free module morphism defined by the matrix
            [0 0 0]
            [0 0 0]
            [0 0 0]
            Domain:   Free module of degree 3 and rank 3 over Integer Ring
                      Echelon ...
            Codomain: Free module of degree 3 and rank 3 over Integer Ring
                      Echelon ...

        The following tests the bug fixed in :issue:`31818`. If there is no
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
            Domain:   Vector space of dimension 2 over Rational Field
            Codomain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: [h(v) for v in V.gens()]
            [(0, 0), (0, 0)]
        """
        from . import free_module_morphism
        side = kwds.get("side", "left")
        if not isinstance(A, Matrix):
            # Compute the matrix of the morphism that sends the
            # generators of the domain to the elements of A.
            C = self.codomain()
            try:
                if callable(A):
                    v = [C(A(g)) for g in self.domain().gens()]
                    A = matrix([C.coordinates(a) for a in v], ncols=C.rank())
                    if side == "right":
                        A = A.transpose()
                else:
                    v = [C(a) for a in A]
                    if side == "right":
                        A = matrix([C.coordinates(a) for a in v], ncols=C.rank()).transpose()
                    else:
                        A = matrix([C.coordinates(a) for a in v], ncols=C.rank())
            except TypeError:
                # Let us hope that FreeModuleMorphism knows to handle
                # that case
                pass
        if not self.codomain().base_ring().has_coerce_map_from(self.domain().base_ring()) and not A.is_zero():
            raise TypeError("nontrivial morphisms require a coercion map from the base ring of the domain to the base ring of the codomain")
        return free_module_morphism.FreeModuleMorphism(self, A, side)

    @cached_method
    def zero(self, side='left'):
        """
        INPUT:

        - ``side`` -- side of the vectors acted on by the matrix
          (default: ``'left'``)

        EXAMPLES::

            sage: E = ZZ^2
            sage: F = ZZ^3
            sage: H = Hom(E, F)
            sage: f = H.zero()
            sage: f
            Free module morphism defined by the matrix
            [0 0 0]
            [0 0 0]
            Domain:   Ambient free module of rank 2 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: f(E.an_element())
            (0, 0, 0)
            sage: f(E.an_element()) == F.zero()
            True
            sage: H.zero("right")
            Free module morphism defined as left-multiplication by the matrix
            [0 0]
            [0 0]
            [0 0]
            Domain:   Ambient free module of rank 2 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring


        TESTS:

        We check that ``H.zero()`` is picklable::

            sage: loads(dumps(f.parent().zero()))
            Free module morphism defined by the matrix
            [0 0 0]
            [0 0 0]
            Domain:   Ambient free module of rank 2 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return self(lambda x: self.codomain().zero(), side=side)

    @cached_method
    def _matrix_space(self, side='left'):
        """
        INPUT:

        - ``side`` -- side of the vectors acted on by the matrix
          (default: ``'left'``)

        Return underlying matrix space that contains the matrices that define
        the homomorphisms in this free module homspace.

        OUTPUT: matrix space

        EXAMPLES::

            sage: H = Hom(QQ^3, QQ^2)
            sage: H._matrix_space()
            Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        """
        if side not in ["left", "right"]:
            raise ValueError("the side must be either 'left' or 'right'")
        R = self.codomain().base_ring()
        if side == "left":
            return MatrixSpace(R, self.domain().rank(), self.codomain().rank())
        elif side == "right":
            return MatrixSpace(R, self.codomain().rank(), self.domain().rank())

    @cached_method
    def basis(self, side='left'):
        """
        Return a basis for this space of free module homomorphisms.

        INPUT:

        - ``side`` -- side of the vectors acted on by the matrix
          (default: ``'left'``)

        OUTPUT: tuple

        EXAMPLES::

            sage: H = Hom(ZZ^2, ZZ^1)
            sage: H.basis()
            (Free module morphism defined by the matrix
              [1]
              [0]
              Domain:   Ambient free module of rank 2 over the principal ideal domain ...
              Codomain: Ambient free module of rank 1 over the principal ideal domain ...,
             Free module morphism defined by the matrix
              [0]
              [1]
              Domain:   Ambient free module of rank 2 over the principal ideal domain ...
              Codomain: Ambient free module of rank 1 over the principal ideal domain ...)
            sage: H.basis("right")
            (Free module morphism defined as left-multiplication by the matrix
              [1 0]
              Domain:   Ambient free module of rank 2 over the principal ideal domain ...
              Codomain: Ambient free module of rank 1 over the principal ideal domain ...,
             Free module morphism defined as left-multiplication by the matrix
              [0 1]
              Domain:   Ambient free module of rank 2 over the principal ideal domain ...
              Codomain: Ambient free module of rank 1 over the principal ideal domain ...)
        """
        M = self._matrix_space(side)
        B = M.basis()
        return tuple([self(x, side=side) for x in B])

    def identity(self, side='left'):
        r"""
        Return identity morphism in an endomorphism ring.

        INPUT:

        - ``side`` -- side of the vectors acted on by the matrix
          (default: ``'left'``)

        EXAMPLES::

            sage: V = FreeModule(ZZ,5)
            sage: H = V.Hom(V)
            sage: H.identity()
            Free module morphism defined by the matrix
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            Domain:   Ambient free module of rank 5 over the principal ideal domain ...
            Codomain: Ambient free module of rank 5 over the principal ideal domain ...
        """
        if self.is_endomorphism_set():
            return self(identity_matrix(self.base_ring(), self.domain().rank()),
                        side=side)
        raise TypeError("Identity map only defined for endomorphisms. Try natural_map() instead.")
