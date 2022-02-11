r"""
Space of Morphisms of Vector Spaces (Linear Transformations)

AUTHOR:

    - Rob Beezer: (2011-06-29)

A :class:`VectorSpaceHomspace` object represents the set of all
possible homomorphisms from one vector space to another.
These mappings are usually known as linear transformations.

For more information on the use of linear transformations,
consult the documentation for vector space morphisms at
:mod:`sage.modules.vector_space_morphism`. Also, this is
an extremely thin veneer on free module homspaces
(:mod:`sage.modules.free_module_homspace`) and free module
morphisms (:mod:`sage.modules.free_module_morphism`) -
objects which might also be useful, and places
where much of the documentation resides.

EXAMPLES:

Creation and basic examination is simple. ::

    sage: V = QQ^3
    sage: W = QQ^2
    sage: H = Hom(V, W)
    sage: H
    Set of Morphisms (Linear Transformations) from
    Vector space of dimension 3 over Rational Field to
    Vector space of dimension 2 over Rational Field
    sage: H.domain()
    Vector space of dimension 3 over Rational Field
    sage: H.codomain()
    Vector space of dimension 2 over Rational Field

Homspaces have a few useful properties.  A basis is provided by
a list of matrix representations, where these matrix representatives
are relative to the bases of the domain and codomain.  ::

    sage: K = Hom(GF(3)^2, GF(3)^2)
    sage: B = K.basis()
    sage: for f in B:
    ....:     print(f)
    ....:     print("\n")
    Vector space morphism represented by the matrix:
    [1 0]
    [0 0]
    Domain: Vector space of dimension 2 over Finite Field of size 3
    Codomain: Vector space of dimension 2 over Finite Field of size 3
    <BLANKLINE>
    Vector space morphism represented by the matrix:
    [0 1]
    [0 0]
    Domain: Vector space of dimension 2 over Finite Field of size 3
    Codomain: Vector space of dimension 2 over Finite Field of size 3
    <BLANKLINE>
    Vector space morphism represented by the matrix:
    [0 0]
    [1 0]
    Domain: Vector space of dimension 2 over Finite Field of size 3
    Codomain: Vector space of dimension 2 over Finite Field of size 3
    <BLANKLINE>
    Vector space morphism represented by the matrix:
    [0 0]
    [0 1]
    Domain: Vector space of dimension 2 over Finite Field of size 3
    Codomain: Vector space of dimension 2 over Finite Field of size 3
    <BLANKLINE>

The zero and identity mappings are properties of the space.
The identity mapping will only be available if the domain and codomain
allow for endomorphisms (equal vector spaces with equal bases).  ::

    sage: H = Hom(QQ^3, QQ^3)
    sage: g = H.zero()
    sage: g([1, 1/2, -3])
    (0, 0, 0)
    sage: f = H.identity()
    sage: f([1, 1/2, -3])
    (1, 1/2, -3)

The homspace may be used with various representations of a
morphism in the space to create the morphism.  We demonstrate
three ways to create the same linear transformation between
two two-dimensional subspaces of ``QQ^3``.  The ``V.n`` notation
is a shortcut to the generators of each vector space, better
known as the basis elements.  Note that the matrix representations
are relative to the bases, which are purposely fixed when the
subspaces are created ("user bases").  ::

    sage: U = QQ^3
    sage: V = U.subspace_with_basis([U.0+U.1, U.1-U.2])
    sage: W = U.subspace_with_basis([U.0, U.1+U.2])
    sage: H = Hom(V, W)

First, with a matrix.  Note that the matrix representation
acts by matrix multiplication with the vector on the left.
The input to the linear transformation, ``(3, 1, 2)``,
is converted to the coordinate vector ``(3, -2)``, then
matrix multiplication yields the vector ``(-3, -2)``,
which represents the vector ``(-3, -2, -2)`` in the codomain.  ::

    sage: m = matrix(QQ, [[1, 2], [3, 4]])
    sage: f1 = H(m)
    sage: f1
    Vector space morphism represented by the matrix:
    [1 2]
    [3 4]
    Domain: Vector space of degree 3 and dimension 2 over Rational Field
    User basis matrix:
    [ 1  1  0]
    [ 0  1 -1]
    Codomain: Vector space of degree 3 and dimension 2 over Rational Field
    User basis matrix:
    [1 0 0]
    [0 1 1]
    sage: f1([3,1,2])
    (-3, -2, -2)

Second, with a list of images of the domain's basis elements.  ::

    sage: img = [1*(U.0) + 2*(U.1+U.2), 3*U.0 + 4*(U.1+U.2)]
    sage: f2 = H(img)
    sage: f2
    Vector space morphism represented by the matrix:
    [1 2]
    [3 4]
    Domain: Vector space of degree 3 and dimension 2 over Rational Field
    User basis matrix:
    [ 1  1  0]
    [ 0  1 -1]
    Codomain: Vector space of degree 3 and dimension 2 over Rational Field
    User basis matrix:
    [1 0 0]
    [0 1 1]
    sage: f2([3,1,2])
    (-3, -2, -2)

Third, with a linear function taking the domain to the codomain.  ::

    sage: g = lambda x: vector(QQ, [-2*x[0]+3*x[1], -2*x[0]+4*x[1], -2*x[0]+4*x[1]])
    sage: f3 = H(g)
    sage: f3
    Vector space morphism represented by the matrix:
    [1 2]
    [3 4]
    Domain: Vector space of degree 3 and dimension 2 over Rational Field
    User basis matrix:
    [ 1  1  0]
    [ 0  1 -1]
    Codomain: Vector space of degree 3 and dimension 2 over Rational Field
    User basis matrix:
    [1 0 0]
    [0 1 1]
    sage: f3([3,1,2])
    (-3, -2, -2)

The three linear transformations look the same, and are the same.  ::

    sage: f1 == f2
    True
    sage: f2 == f3
    True

TESTS::

    sage: V = QQ^2
    sage: W = QQ^3
    sage: H = Hom(QQ^2, QQ^3)
    sage: loads(dumps(H))
    Set of Morphisms (Linear Transformations) from
    Vector space of dimension 2 over Rational Field to
    Vector space of dimension 3 over Rational Field
    sage: loads(dumps(H)) == H
    True
"""

####################################################################################
#       Copyright (C) 2011 Rob Beezer <beezer@ups.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
####################################################################################

import sage.matrix.all as matrix
import sage.modules.free_module_homspace

# This module initially overrides just the minimum functionality necessary
# from  sage.modules.free_module_homspace.FreeModuleHomSpace.
# If additional methods here override the free module homspace methods,
# consider adjusting the free module doctests, since many are written with
# examples that are actually vector spaces and not so many use "pure" modules
# for the examples.


def is_VectorSpaceHomspace(x):
    r"""
    Return ``True`` if ``x`` is a vector space homspace.

    INPUT:

    ``x`` - anything

    EXAMPLES:

    To be a vector space morphism, the domain and codomain must both be
    vector spaces, in other words, modules over fields.  If either
    set is just a module, then the ``Hom()`` constructor will build a
    space of free module morphisms.  ::

        sage: H = Hom(QQ^3, QQ^2)
        sage: type(H)
        <class 'sage.modules.vector_space_homspace.VectorSpaceHomspace_with_category'>
        sage: sage.modules.vector_space_homspace.is_VectorSpaceHomspace(H)
        True

        sage: K = Hom(QQ^3, ZZ^2)
        sage: type(K)
        <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        sage: sage.modules.vector_space_homspace.is_VectorSpaceHomspace(K)
        False

        sage: L = Hom(ZZ^3, QQ^2)
        sage: type(L)
        <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        sage: sage.modules.vector_space_homspace.is_VectorSpaceHomspace(L)
        False

        sage: sage.modules.vector_space_homspace.is_VectorSpaceHomspace('junk')
        False
    """
    return isinstance(x, VectorSpaceHomspace)


class VectorSpaceHomspace(sage.modules.free_module_homspace.FreeModuleHomspace):

    def __call__(self, A, check=True, **kwds):
        r"""
        INPUT:

        - ``A`` - one of several possible inputs representing
          a morphism from this vector space homspace.
          - a vector space morphism in this homspace
          - a matrix representation relative to the bases of the vector spaces,
            which acts on a vector placed to the left of the matrix
          - a list or tuple containing images of the domain's basis vectors
          - a function from the domain to the codomain
        - ``check`` (default: True) - ``True`` or ``False``, required for
          compatibility with calls from
          :meth:`sage.structure.parent.Parent.hom`.
        - the keyword ``side`` can be assigned the values ``"left"`` or
          ``"right"``. It corresponds to the side of vectors relative to the
          matrix.

        EXAMPLES::

            sage: V = (QQ^3).span_of_basis([[1,1,0],[1,0,2]])
            sage: H = V.Hom(V)
            sage: H
            Set of Morphisms (Linear Transformations) from
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]
            to
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]

        Coercing a matrix::

            sage: A = matrix(QQ, [[0, 1], [1, 0]])
            sage: rho = H(A)          # indirect doctest
            sage: rho
            Vector space morphism represented by the matrix:
            [0 1]
            [1 0]
            Domain: Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]
            Codomain: Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]

        Coercing a list of images::

            sage: phi = H([V.1, V.0])
            sage: phi(V.1) == V.0
            True
            sage: phi(V.0) == V.1
            True
            sage: phi
            Vector space morphism represented by the matrix:
            [0 1]
            [1 0]
            Domain: Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]
            Codomain: Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]

        Coercing a lambda function::

            sage: f = lambda x: vector(QQ, [x[0], (1/2)*x[2], 2*x[1]])
            sage: zeta = H(f)
            sage: zeta
            Vector space morphism represented by the matrix:
            [0 1]
            [1 0]
            Domain: Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]
            Codomain: Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 1 0]
            [1 0 2]

        Coercing a vector space morphism into the parent of a second vector
        space morphism will unify their parents::

            sage: U = FreeModule(QQ,3, sparse=True ); V = QQ^4
            sage: W = FreeModule(QQ,3, sparse=False); X = QQ^4
            sage: H = Hom(U, V)
            sage: K = Hom(W, X)
            sage: H is K, H == K
            (False, True)

            sage: A = matrix(QQ, 3, 4, [0]*12)
            sage: f = H(A)
            sage: B = matrix(QQ, 3, 4, range(12))
            sage: g = K(B)
            sage: f.parent() is H and g.parent() is K
            True

            sage: h = H(g)
            sage: f.parent() is h.parent()
            True

        See other examples in the module-level documentation.

        TESTS::

            sage: V = GF(3)^0
            sage: W = GF(3)^1
            sage: H = V.Hom(W)
            sage: H.zero().is_zero()
            True

        Previously the above code resulted in a TypeError because the
        dimensions of the matrix were incorrect.
        """
        from .vector_space_morphism import is_VectorSpaceMorphism, VectorSpaceMorphism
        D = self.domain()
        C = self.codomain()
        side = kwds.get("side", "left")
        from sage.structure.element import is_Matrix
        if is_Matrix(A):
            pass
        elif is_VectorSpaceMorphism(A):
            A = A.matrix()
        elif callable(A):
            try:
                images = [A(g) for g in D.basis()]
            except (ValueError, TypeError, IndexError) as e:
                msg = 'function cannot be applied properly to some basis element because\n' + e.args[0]
                raise ValueError(msg)
            try:
                A = matrix.matrix(D.dimension(), C.dimension(), [C.coordinates(C(a)) for a in images])
            except (ArithmeticError, TypeError) as e:
                msg = 'some image of the function is not in the codomain, because\n' + e.args[0]
                raise ArithmeticError(msg)
            if side == "right":
                A = A.transpose()
        elif isinstance(A, (list, tuple)):
            if len(A) != len(D.basis()):
                msg = "number of images should equal the size of the domain's basis (={0}), not {1}"
                raise ValueError(msg.format(len(D.basis()), len(A)))
            try:
                v = [C(a) for a in A]
                A = matrix.matrix(D.dimension(), C.dimension(), [C.coordinates(a) for a in v])
            except (ArithmeticError, TypeError) as e:
                msg = 'some proposed image is not in the codomain, because\n' + e.args[0]
                raise ArithmeticError(msg)
            if side == "right":
                A = A.transpose()
        else:
            msg = 'vector space homspace can only coerce matrices, vector space morphisms, functions or lists, not {0}'
            raise TypeError(msg.format(A))
        return VectorSpaceMorphism(self, A, side=side)

    def _repr_(self):
        r"""
        Text representation of a space of vector space morphisms.

        EXAMPLES::

            sage: H = Hom(QQ^2, QQ^3)
            sage: H._repr_().split(' ')
            ['Set', 'of', 'Morphisms', '(Linear', 'Transformations)',
            'from', 'Vector', 'space', 'of', 'dimension', '2', 'over',
            'Rational', 'Field', 'to', 'Vector', 'space', 'of',
            'dimension', '3', 'over', 'Rational', 'Field']
        """
        msg = 'Set of Morphisms (Linear Transformations) from {0} to {1}'
        return msg.format(self.domain(), self.codomain())
