"""
Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.indexed_generators import standardize_names_index_set
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras, LiftMorphism
from sage.categories.rings import Rings
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom

from sage.algebras.lie_algebras.lie_algebra_element import (LieAlgebraElementWrapper,
                                                            LieAlgebraMatrixWrapper)
from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.sets.family import Family, AbstractFamily


class LieAlgebra(Parent, UniqueRepresentation):  # IndexedGenerators):
    r"""
    A Lie algebra `L` over a base ring `R`.

    A Lie algebra is an `R`-module `L` with a bilinear operation called
    Lie bracket `[\cdot, \cdot] : L \times L \to L` such that
    `[x, x] = 0` and the following relation holds:

    .. MATH::

        \bigl[ x, [y, z] \bigr] + \bigl[ y, [z, x] \bigr]
        + \bigl[ z, [x, y] \bigr] = 0.

    This relation is known as the *Jacobi identity* (or sometimes the Jacobi
    relation). We note that from `[x, x] = 0`, we have `[x + y, x + y] = 0`.
    Next from bilinearity, we see that

    .. MATH::

        0 = [x + y, x + y] = [x, x] + [x, y] + [y, x] + [y, y]
        = [x, y] + [y, x],

    thus `[x, y] = -[y, x]` and the Lie bracket is antisymmetric.

    Lie algebras are closely related to Lie groups. Let `G` be a Lie group
    and fix some `g \in G`. We can construct the Lie algebra `L` of `G` by
    considering the tangent space at `g`. We can also (partially) recover `G`
    from `L` by using what is known as the exponential map.

    Given any associative algebra `A`, we can construct a Lie algebra `L`
    on the `R`-module `A` by defining the Lie bracket to be the commutator
    `[a, b] = ab - ba`. We call an associative algebra `A` which contains
    `L` in this fashion an *enveloping algebra* of `L`. The embedding
    `L \to A` which sends the Lie bracket to the commutator will be called
    a Lie embedding. Now if we are given a Lie algebra `L`, we
    can construct an enveloping algebra `U_L` with Lie embedding `h : L \to
    U_L` which has the following universal property: for any enveloping
    algebra `A` with Lie embedding `f : L \to A`, there exists a unique unital
    algebra homomorphism `g : U_L \to A` such that `f = g \circ h`. The
    algebra `U_L` is known as the *universal enveloping algebra* of `L`.

    INPUT:

    See examples below for various input options.

    EXAMPLES:

    **1.** The simplest examples of Lie algebras are *abelian Lie
    algebras*. These are Lie algebras whose Lie bracket is (identically)
    zero. We can create them using the ``abelian`` keyword::

        sage: L.<x,y,z> = LieAlgebra(QQ, abelian=True); L
        Abelian Lie algebra on 3 generators (x, y, z) over Rational Field

    **2.** A Lie algebra can be built from any associative algebra by
    defining the Lie bracket to be the commutator. For example, we can
    start with the descent algebra::

        sage: D = DescentAlgebra(QQ, 4).D()
        sage: L = LieAlgebra(associative=D); L
        Lie algebra of Descent algebra of 4 over Rational Field
         in the standard basis
        sage: L(D[2]).bracket(L(D[3]))
        D{1, 2} - D{1, 3} + D{2} - D{3}

    Next we use a free algebra and do some simple computations::

        sage: R.<a,b,c> = FreeAlgebra(QQ, 3)
        sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
        sage: x-y+z
        a - b + c
        sage: L.bracket(x-y, x-z)
        a*b - a*c - b*a + b*c + c*a - c*b
        sage: L.bracket(x-y, L.bracket(x,y))
        a^2*b - 2*a*b*a + a*b^2 + b*a^2 - 2*b*a*b + b^2*a

    We can also use a subset of the elements as a generating set
    of the Lie algebra::

        sage: R.<a,b,c> = FreeAlgebra(QQ, 3)
        sage: L.<x,y> = LieAlgebra(associative=[a,b+c])
        sage: L.bracket(x, y)
        a*b + a*c - b*a - c*a

    Now for a more complicated example using the group ring of `S_3` as our
    base algebra::

        sage: G = SymmetricGroup(3)
        sage: S = GroupAlgebra(G, QQ)
        sage: L.<x,y> = LieAlgebra(associative=S.gens())
        sage: L.bracket(x, y)
        (2,3) - (1,3)
        sage: L.bracket(x, y-x)
        (2,3) - (1,3)
        sage: L.bracket(L.bracket(x, y), y)
        2*(1,2,3) - 2*(1,3,2)
        sage: L.bracket(x, L.bracket(x, y))
        (2,3) - 2*(1,2) + (1,3)
        sage: L.bracket(x, L.bracket(L.bracket(x, y), y))
        0

    Here is an example using matrices::

        sage: MS = MatrixSpace(QQ,2)
        sage: m1 = MS([[0, -1], [1, 0]])
        sage: m2 = MS([[-1, 4], [3, 2]])
        sage: L.<x,y> = LieAlgebra(associative=[m1, m2])
        sage: x
        [ 0 -1]
        [ 1  0]
        sage: y
        [-1  4]
        [ 3  2]
        sage: L.bracket(x,y)
        [-7 -3]
        [-3  7]
        sage: L.bracket(y,y)
        [0 0]
        [0 0]
        sage: L.bracket(y,x)
        [ 7  3]
        [ 3 -7]
        sage: L.bracket(x, L.bracket(y,x))
        [-6 14]
        [14  6]

    (See :class:`LieAlgebraFromAssociative` for other examples.)

    **3.** We can also creating a Lie algebra by inputting a set of
    structure coefficients. For example, we can create the Lie algebra
    of `\QQ^3` under the Lie bracket `\times` (cross-product)::

        sage: d = {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}}
        sage: L.<x,y,z> = LieAlgebra(QQ, d)
        sage: L
        Lie algebra on 3 generators (x, y, z) over Rational Field

    To compute the Lie bracket of two elements, you cannot use the ``*``
    operator. Indeed, this automatically lifts up to the universal
    enveloping algebra and takes the (associative) product there.
    To get elements in the Lie algebra, you must use :meth:`bracket`::

        sage: L = LieAlgebra(QQ, {('e','h'): {'e':-2}, ('f','h'): {'f':2},
        ....:                     ('e','f'): {'h':1}}, names='e,f,h')
        sage: e,f,h = L.lie_algebra_generators()
        sage: L.bracket(h, e)
        2*e
        sage: elt = h*e; elt
        e*h + 2*e
        sage: P = elt.parent(); P
        Noncommutative Multivariate Polynomial Ring in e, f, h over Rational Field,
         nc-relations: {...}
        sage: R = P.relations()
        sage: for rhs in sorted(R, key=str): print("{} = {}".format(rhs, R[rhs]))
        f*e = e*f - h
        h*e = e*h + 2*e
        h*f = f*h - 2*f

    For convenience, there are two shorthand notations for computing
    Lie brackets::

        sage: L([h,e])
        2*e
        sage: L([h,[e,f]])
        0
        sage: L([[h,e],[e,f]])
        -4*e
        sage: L[h, e]
        2*e
        sage: L[h, L[e, f]]
        0

    .. WARNING::

        Because this is a modified (abused) version of python syntax, it
        does **NOT** work with addition. For example ``L([e + [h, f], h])``
        and ``L[e + [h, f], h]`` will both raise errors. Instead you must
        use ``L[e + L[h, f], h]``.

    **4.** We can construct a Lie algebra from a Cartan type by using
    the ``cartan_type`` option::

        sage: L = LieAlgebra(ZZ, cartan_type=['C', 3])
        sage: L.inject_variables()
        Defining e1, e2, e3, f1, f2, f3, h1, h2, h3
        sage: e1.bracket(e2)
        -E[alpha[1] + alpha[2]]
        sage: L([[e1, e2], e2])
        0
        sage: L([[e2, e3], e3])
        0
        sage: L([e2, [e2, e3]])
        2*E[2*alpha[2] + alpha[3]]

        sage: L = LieAlgebra(ZZ, cartan_type=['E', 6])
        sage: L
        Lie algebra of ['E', 6] in the Chevalley basis

    When the Cartan type is finite type and simply-laced, we can also
    specify an asymmetry function from [Ka1990]_ using a Dynkin diagram
    orientation with the ``epsilon`` option::

        sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(1, 2)])
        sage: e1, e2 = L.e()
        sage: L[e1, e2]
        -E[alpha[1] + alpha[2]]

        sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(2, 1)])
        sage: e1, e2 = L.e()
        sage: L[e1, e2]
        E[alpha[1] + alpha[2]]

    We also have matrix versions of the classical Lie algebras::

        sage: L = LieAlgebra(ZZ, cartan_type=['A', 2], representation='matrix')
        sage: L.gens()
        (
        [0 1 0]  [0 0 0]  [0 0 0]  [0 0 0]  [ 1  0  0]  [ 0  0  0]
        [0 0 0]  [0 0 1]  [1 0 0]  [0 0 0]  [ 0 -1  0]  [ 0  1  0]
        [0 0 0], [0 0 0], [0 0 0], [0 1 0], [ 0  0  0], [ 0  0 -1]
        )

    There is also the compact real form of matrix Lie algebras
    implemented (the base ring must currently be a field)::

        sage: L = LieAlgebra(QQ, cartan_type=['A', 2], representation="compact real")
        sage: list(L.basis())
        [
        [ 0  1  0]  [ 0  0  1]  [ 0  0  0]  [ i  0  0]  [0 i 0]  [0 0 i]
        [-1  0  0]  [ 0  0  0]  [ 0  0  1]  [ 0  0  0]  [i 0 0]  [0 0 0]
        [ 0  0  0], [-1  0  0], [ 0 -1  0], [ 0  0 -i], [0 0 0], [i 0 0],
        <BLANKLINE>
        [ 0  0  0]  [0 0 0]
        [ 0  i  0]  [0 0 i]
        [ 0  0 -i], [0 i 0]
        ]

    **5.** We construct a free Lie algebra in a few different ways. There are
    two primary representations, as brackets and as polynomials::

        sage: L = LieAlgebra(QQ, 'x,y,z'); L
        Free Lie algebra generated by (x, y, z) over Rational Field
        sage: P.<a,b,c> = LieAlgebra(QQ, representation='polynomial'); P
        Lie algebra generated by (a, b, c) in
         Free Algebra on 3 generators (a, b, c) over Rational Field

    This has the basis given by Hall and the one indexed by Lyndon words.
    We do some computations and convert between the bases::

        sage: H = L.Hall()
        doctest:warning...:
        FutureWarning: The Hall basis has not been fully proven correct, but currently no bugs are known
        See https://github.com/sagemath/sage/issues/16823 for details.
        sage: H
        Free Lie algebra generated by (x, y, z) over Rational Field in the Hall basis
        sage: Lyn = L.Lyndon()
        sage: Lyn
        Free Lie algebra generated by (x, y, z) over Rational Field in the Lyndon basis
        sage: x,y,z = Lyn.lie_algebra_generators()
        sage: a = Lyn([x, [[z, [x, y]], [y, x]]]); a
        -[x, [[x, y], [x, [y, z]]]] - [x, [[x, y], [[x, z], y]]]
        sage: H(a)
        [[x, y], [z, [x, [x, y]]]] - [[x, y], [[x, y], [x, z]]]
         + [[x, [x, y]], [z, [x, y]]]

    We also have the free Lie algebra given in the polynomial
    representation, which is the canonical embedding of the free
    Lie algebra into the free algebra (i.e., the ring of
    noncommutative polynomials).
    So the generators of the free Lie algebra are the generators of the
    free algebra and the Lie bracket is the commutator::

        sage: P.<a,b,c> = LieAlgebra(QQ, representation='polynomial'); P
        Lie algebra generated by (a, b, c) in
         Free Algebra on 3 generators (a, b, c) over Rational Field
        sage: P.bracket(a, b) + P.bracket(a - c, b + 3*c)
        2*a*b + 3*a*c - 2*b*a + b*c - 3*c*a - c*b

    **6.** Nilpotent Lie algebras are Lie algebras such that there exists an
    integer `s` such that all iterated brackets of length longer than `s`
    are zero. They can be constructed from structural coefficients using the
    ``nilpotent`` keyword::

        sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
        sage: L
        Nilpotent Lie algebra on 3 generators (X, Y, Z) over Rational Field
        sage: L.category()
        Category of finite dimensional nilpotent Lie algebras with basis over Rational Field

    A second example defining the Engel Lie algebra::

        sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
        sage: E.<X,Y,Z,W> = LieAlgebra(QQ, sc, nilpotent=True); E
        Nilpotent Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: E.step()
        3
        sage: E[X, Y + Z]
        Z + W
        sage: E[X, [X, Y + Z]]
        W
        sage: E[X, [X, [X, Y + Z]]]
        0

    A nilpotent Lie algebra will also be constructed if given a ``category``
    of a nilpotent Lie algebra::

        sage: C = LieAlgebras(QQ).Nilpotent().FiniteDimensional().WithBasis()
        sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, category=C); L
        Nilpotent Lie algebra on 3 generators (X, Y, Z) over Rational Field

    **7.** Free nilpotent Lie algebras are the truncated versions of the free
    Lie algebras. That is, the only relations other than anticommutativity
    and the Jacobi identity among the Lie brackets are that brackets of
    length higher than the nilpotency step vanish. They can be created by
    using the ``step`` keyword::

        sage: L = LieAlgebra(ZZ, 2, step=3); L
        Free Nilpotent Lie algebra on 5 generators (X_1, X_2, X_12, X_112, X_122) over Integer Ring
        sage: L.step()
        3

    REFERENCES:

    - [deG2000]_ Willem A. de Graaf. *Lie Algebras: Theory and Algorithms*.
    - [Ka1990]_ Victor Kac, *Infinite dimensional Lie algebras*.
    - :wikipedia:`Lie_algebra`
    """
    # This works because it is an abstract base class and this
    #    __classcall_private__ will only be called when calling LieAlgebra
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, arg1=None, names=None,
                              index_set=None, abelian=False, nilpotent=False,
                              category=None, **kwds):
        """
        Select the correct parent based upon input.

        TESTS::

            sage: LieAlgebra(QQ, abelian=True, names='x,y,z')
            Abelian Lie algebra on 3 generators (x, y, z) over Rational Field
            sage: LieAlgebra(QQ, {('e','h'): {'e':-2}, ('f','h'): {'f':2},
            ....:                 ('e','f'): {'h':1}}, names='e,f,h')
            Lie algebra on 3 generators (e, f, h) over Rational Field
        """
        # Parse associative algebra input
        # -----

        assoc = kwds.get("associative", None)
        if assoc is not None:
            return LieAlgebraFromAssociative(assoc, names=names, index_set=index_set,
                                             category=category)

        # Parse input as a Cartan type
        # -----

        ct = kwds.pop("cartan_type", None)
        if ct is not None:
            from sage.combinat.root_system.cartan_type import CartanType
            ct = CartanType(ct)
            if ct.is_affine():
                from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
                return AffineLieAlgebra(R, cartan_type=ct,
                                        kac_moody=kwds.get("kac_moody", True))
            if not ct.is_finite():
                raise NotImplementedError("non-finite types are not implemented yet, see trac #14901 for details")
            rep = kwds.pop("representation", "bracket")
            if rep == 'bracket':
                from sage.algebras.lie_algebras.classical_lie_algebra import LieAlgebraChevalleyBasis
                return LieAlgebraChevalleyBasis(R, ct, **kwds)
            if rep == 'matrix':
                from sage.algebras.lie_algebras.classical_lie_algebra import ClassicalMatrixLieAlgebra
                return ClassicalMatrixLieAlgebra(R, ct, **kwds)
            if rep == 'compact real':
                from sage.algebras.lie_algebras.classical_lie_algebra import MatrixCompactRealForm
                return MatrixCompactRealForm(R, ct, **kwds)
            raise ValueError("invalid representation")

        # Parse the remaining arguments
        # -----

        if R is None:
            raise ValueError("invalid arguments")

        def check_assoc(A):
            return (isinstance(A, MatrixSpace)
                    or A in Rings()
                    or A in Algebras(R).Associative())
        if arg0 in ZZ or check_assoc(arg1):
            # Check if we need to swap the arguments
            arg0, arg1 = arg1, arg0

        # Parse the first argument
        # -----

        if isinstance(arg0, dict):
            if not arg0:
                from sage.algebras.lie_algebras.abelian import AbelianLieAlgebra
                return AbelianLieAlgebra(R, names, index_set)
            elif isinstance(next(iter(arg0.keys())), (list, tuple)):
                # We assume it is some structure coefficients
                arg1, arg0 = arg0, arg1

        if isinstance(arg0, (list, tuple)):
            if all(isinstance(x, str) for x in arg0):
                # If they are all strings, then it is a list of variables
                names = tuple(arg0)

        if isinstance(arg0, str):
            names = tuple(arg0.split(','))
        elif isinstance(names, str):
            names = tuple(names.split(','))

        # Parse the second argument

        if isinstance(arg1, dict):
            # Assume it is some structure coefficients
            if nilpotent or (category is not None and category.is_subcategory(LieAlgebras(R).Nilpotent())):
                from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra_dense
                return NilpotentLieAlgebra_dense(R, arg1, names, index_set,
                                                 category=category, **kwds)

            from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
            return LieAlgebraWithStructureCoefficients(R, arg1, names, index_set,
                                                       category=category, **kwds)

        # Otherwise it must be either a free (nilpotent) or abelian Lie algebra

        if arg1 in ZZ:
            step = kwds.get("step", None)
            if step:
                # Parse input as a free nilpotent Lie algebra
                from sage.algebras.lie_algebras.nilpotent_lie_algebra import FreeNilpotentLieAlgebra
                del kwds["step"]
                return FreeNilpotentLieAlgebra(R, arg1, step, names=names, **kwds)
            elif nilpotent:
                raise ValueError("free nilpotent Lie algebras must have a"
                                 " 'step' parameter given")

            if isinstance(arg0, str):
                names = arg0
            if names is None:
                index_set = list(range(arg1))
            else:
                if isinstance(names, str):
                    names = tuple(names.split(','))
                    if arg1 != 1 and len(names) == 1:
                        names = tuple('{}{}'.format(names[0], i)
                                      for i in range(arg1))
                if arg1 != len(names):
                    raise ValueError("the number of names must equal the"
                                     " number of generators")

        if "step" in kwds or nilpotent:
            raise ValueError("free nilpotent Lie algebras must have both"
                             " a number of generators and step parameters"
                             " specified")

        if abelian:
            from sage.algebras.lie_algebras.abelian import AbelianLieAlgebra
            return AbelianLieAlgebra(R, names, index_set)

        # Otherwise it is the free Lie algebra
        rep = kwds.get("representation", "bracket")
        if rep == "polynomial":
            # Construct the free Lie algebra from polynomials in the
            #   free (associative unital) algebra
            # TODO: Change this to accept an index set once FreeAlgebra accepts one
            from sage.algebras.free_algebra import FreeAlgebra
            F = FreeAlgebra(R, names)
            if index_set is None:
                index_set = F.variable_names()
            # TODO: As part of #16823, this should instead construct a
            #   subclass with specialized methods for the free Lie algebra
            return LieAlgebraFromAssociative(F, F.gens(), names=names, index_set=index_set)

        from sage.algebras.lie_algebras.free_lie_algebra import FreeLieAlgebra
        return FreeLieAlgebra(R, names, index_set)

    def __init__(self, R, names=None, category=None):
        """
        The Lie algebra.

        INPUT:

        - ``R`` -- the base ring

        - ``names`` -- (optional) the names of the generators

        - ``category`` -- the category of the Lie algebra; the default is the
          category of Lie algebras over ``R``

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.category()
            Category of finite dimensional nilpotent Lie algebras with basis over Rational Field
        """
        category = LieAlgebras(R).or_subcategory(category)
        Parent.__init__(self, base=R, names=names, category=category)

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation='polynomial')
            sage: elt = L([x, y]); elt
            x*y - y*x
            sage: elt.parent() is L
            True

        TESTS:

        Check that `0` gives the zero element::

            sage: L = lie_algebras.pwitt(GF(5), 5)
            sage: L(0)
            0
        """
        if isinstance(x, list) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))

        if x == 0:
            return self.zero()

        try:
            if x in self.module():
                return self.from_vector(x)
        except AttributeError:
            pass

        if x in self.base_ring():
            # We have already handled the case when x == 0
            raise ValueError("can only convert the scalar 0 into a Lie algebra element")

        return self.element_class(self, x)

    def __getitem__(self, x):
        """
        If ``x`` is a pair `(a, b)`, return the Lie bracket `[a, b]
        (including if `a` or `b` are Lie (sub)algebras, in which case the
        corresponding ideal is constructed).
        Otherwise try to return the `x`-th element of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation='polynomial')
            sage: L[x, [y, x]]
            -x^2*y + 2*x*y*x - y*x^2

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L[L, L]
            Ideal () of Abelian Lie algebra on 2 generators (x, y) over Rational Field

            sage: L = lie_algebras.Heisenberg(QQ, 1)
            sage: Z = L[L, L]; Z
            Ideal (z) of Heisenberg algebra of rank 1 over Rational Field
            sage: L[Z, L]
            Ideal () of Heisenberg algebra of rank 1 over Rational Field

            sage: p,q,z = L.basis(); (p, q, z)
            (p1, q1, z)
            sage: L[p, L]
            Ideal (p1) of Heisenberg algebra of rank 1 over Rational Field
            sage: L[L, p+q]
            Ideal (p1 + q1) of Heisenberg algebra of rank 1 over Rational Field
        """
        if isinstance(x, tuple) and len(x) == 2:
            # Check if we need to construct an ideal
            if x[0] in LieAlgebras:
                if x[1] in LieAlgebras:
                    return x[0].product_space(x[1])
                return x[0].ideal(x[1])
            elif x[1] in LieAlgebras:
                return x[1].ideal(x[0])
            # Otherwise it is the bracket of two elements
            return self(x[0])._bracket_(self(x[1]))
        return super().__getitem__(x)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self`` and
        ``False`` otherwise.

        The things that coerce into ``self`` are:

        - Lie algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - A module which coerces into the base vector space of ``self``.

        TESTS::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L._coerce_map_from_(L.module())
            True
            sage: L._coerce_map_from_(FreeModule(ZZ, 2))
            True
        """
        if not isinstance(R, LieAlgebra):
            # Should be moved to LieAlgebrasWithBasis somehow since it is a generic coercion
            if self.module is not NotImplemented:
                return self.module().has_coerce_map_from(R)
            return False

        # We check if it is a subalgebra of something that can coerce into ``self``
        # from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
        # if isinstance(R, LieSubalgebra) and self.has_coerce_map_from(R._ambient):
        #    return R.ambient_lift

        # Lie algebras in the same indices over any base that coerces in
        if R._indices != self._indices:
            return False

        return self.base_ring().has_coerce_map_from(R.base_ring())

    def _Hom_(self, Y, category):
        """
        Return the homset from ``self`` to ``Y`` in the category ``category``.

        INPUT:

        - ``Y`` -- a Lie algebra
        - ``category`` -- a subcategory of :class:`LieAlgebras` or ``None``

        The sole purpose of this method is to construct the homset
        as a :class:`~sage.algebras.lie_algebras.morphism.LieAlgebraHomset`.

        This method is not meant to be called directly. Please use
        :func:`sage.categories.homset.Hom` instead.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ)
            sage: H = L.Hall()
            sage: Hom(H, H)
            Set of Lie algebra morphisms from
             Free Lie algebra generated by (x, y) over Rational Field in the Hall basis
             to Free Lie algebra generated by (x, y) over Rational Field in the Hall basis
        """
        cat = LieAlgebras(self.base_ring())
        if category is not None and not category.is_subcategory(cat):
            raise TypeError(f"{category} is not a subcategory of Lie algebras")
        if Y not in cat:
            raise TypeError(f"{Y} is not a Lie algebra")
        from sage.algebras.lie_algebras.morphism import LieAlgebraHomset
        return LieAlgebraHomset(self, Y, category=category)

    @cached_method
    def zero(self):
        """
        Return the element `0`.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation='polynomial')
            sage: L.zero()
            0
        """
        return self.element_class(self, {})

    # The following methods should belong to ModulesWithBasis?
    def _from_dict(self, d, coerce=False, remove_zeros=True):
        """
        Construct an element of ``self`` from an ``{index: coefficient}``
        dictionary.

        INPUT:

        - ``d`` -- dictionary ``{index: coeff}`` where each ``index`` is the
          index of a basis element and each ``coeff`` belongs to the
          coefficient ring ``self.base_ring()``

        - ``coerce`` -- boolean (default: ``False``); whether to coerce the
          ``coeff`` to the coefficient ring

        - ``remove_zeros`` -- boolean (default: ``True``); if some
          ``coeff`` may be zero and should therefore be removed

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: d = {'p1': 4, 'q3': 1/2, 'z': -2}
            sage: L._from_dict(d)
            4*p1 + 1/2*q3 - 2*z
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            d = {key: R(coeff) for key, coeff in d.items()}
        if remove_zeros:
            d = {key: coeff for key, coeff in d.items() if coeff}
        return self.element_class(self, d)

    def monomial(self, i):
        """
        Return the monomial indexed by ``i``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.monomial('p1')
            p1
        """
        return self.element_class(self, {i: self.base_ring().one()})

    def term(self, i, c=None):
        """
        Return the term indexed by ``i`` with coefficient ``c``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.term('p1', 4)
            4*p1
        """
        if c is None:
            c = self.base_ring().one()
        else:
            c = self.base_ring()(c)
        return self.element_class(self, {i: c})

    def get_order(self):
        """
        Return an ordering of the basis indices.

        .. TODO::

            Remove this method and in :class:`CombinatorialFreeModule`
            in favor of a method in the category of (finite dimensional)
            modules with basis.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {})
            sage: L.get_order()
            ('x', 'y')
        """
        try:
            return self._basis_ordering
        except AttributeError:
            raise ValueError("the Lie algebra is not finite dimensional with a basis")

    # Element = LieAlgebraElement # Default for all Lie algebras


class LieAlgebraWithGenerators(LieAlgebra):
    """
    A Lie algebra with distinguished generators.
    """
    def __init__(self, R, names=None, index_set=None, category=None, prefix='L', **kwds):
        """
        The Lie algebra.

        INPUT:

        - ``R`` -- the base ring
        - ``names`` -- (optional) the names of the generators
        - ``index_set`` -- (optional) the indexing set
        - ``category`` -- the category of the Lie algebra; the default is the
          category of Lie algebras over ``R``
        - ``prefix`` -- (optional) the prefix for the generator representation
        - any keyword accepted by
          :class:`~sage.structure.indexed_generators.IndexedGenerators`

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.category()
            Category of finite dimensional nilpotent Lie algebras with basis over Rational Field
        """
        self._indices = index_set
        LieAlgebra.__init__(self, R, names, category)

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation='polynomial')
            sage: L.lie_algebra_generators()
            Finite family {'x': x, 'y': y}
        """
        return Family(self._indices, self.monomial, name="monomial map")

    @cached_method
    def gens(self):
        """
        Return a tuple whose entries are the generators for this
        object, in some order.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.gens()
            (x, y)
        """
        G = self.lie_algebra_generators()
        try:
            return tuple(G[i] for i in self.variable_names())
        except (KeyError, IndexError):
            return tuple(G[i] for i in self.indices())
        except ValueError:
            return tuple(G)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.gen(0)
            x
        """
        return self.gens()[i]

    def indices(self):
        """
        Return the indices of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation='polynomial')
            sage: L.indices()
            {'x', 'y'}
        """
        return self._indices


class FinitelyGeneratedLieAlgebra(LieAlgebraWithGenerators):
    r"""
    A finitely generated Lie algebra.
    """
    def __init__(self, R, names=None, index_set=None, category=None):
        """
        Initialize ``self``.

        INPUT:

        - ``R`` -- the base ring

        - ``names`` -- the names of the generators

        - ``index_set`` -- the index set of the generators

        - ``category`` -- the category of the Lie algebra

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.category()
            Category of finite dimensional nilpotent Lie algebras with basis over Rational Field
        """
        LieAlgebraWithGenerators.__init__(self, R, names, index_set, category)
        self.__ngens = len(self._indices)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x,y> = LieAlgebra(QQ, {('x','y'): {'x': 1}})
            sage: F
            Lie algebra on 2 generators (x, y) over Rational Field
        """
        if self.__ngens == 1:
            return "Lie algebra on the generator {} over {}".format(
                self.gen(0), self.base_ring())
        return "Lie algebra on {} generators {} over {}".format(
            self.__ngens, self.gens(), self.base_ring())

    @lazy_attribute
    def _ordered_indices(self):
        """
        Return the index set of the basis of ``self`` in (some) order.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L._ordered_indices
            ('x', 'y')
        """
        return tuple(self.basis().keys())

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.an_element()
            x + y
        """
        return self.sum(self.lie_algebra_generators())


class InfinitelyGeneratedLieAlgebra(LieAlgebraWithGenerators):
    r"""
    An infinitely generated Lie algebra.
    """
    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L._an_element_()
            p2 + q2 - 1/2*q3 + z
        """
        return self.lie_algebra_generators()[self._indices.an_element()]

# Do we want this to return lie_algebra_generators()? Perhaps in the category?
#    def gens(self):
#        """
#        Return a tuple whose entries are the generators for this
#        object, in some order.
#
#        EXAMPLES::
#
#            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
#            sage: L.gens()
#            (x, y)
#        """
#        return self.lie_algebra_generators()


class LieAlgebraFromAssociative(LieAlgebraWithGenerators):
    """
    A Lie algebra whose elements are from an associative algebra and whose
    bracket is the commutator.

    .. TODO::

        Split this class into 2 classes, the base class for the Lie
        algebra corresponding to the full associative algebra and a
        subclass for the Lie subalgebra (of the full algebra)
        generated by a generating set?

    .. TODO::

        Return the subalgebra generated by the basis
        elements of ``self`` for the universal enveloping algebra.

    EXAMPLES:

    For the first example, we start with a commutative algebra.
    Note that the bracket of everything will be 0::

        sage: R = SymmetricGroupAlgebra(QQ, 2)
        sage: L = LieAlgebra(associative=R)
        sage: x, y = L.basis()
        sage: L.bracket(x, y)
        0

    Next we use a free algebra and do some simple computations::

        sage: R.<a,b> = FreeAlgebra(QQ, 2)
        sage: L = LieAlgebra(associative=R)
        sage: x,y = L(a), L(b)
        sage: x-y
        a - b
        sage: L.bracket(x-y, x)
        a*b - b*a
        sage: L.bracket(x-y, L.bracket(x,y))
        a^2*b - 2*a*b*a + a*b^2 + b*a^2 - 2*b*a*b + b^2*a

    We can also use a subset of the generators as a generating set
    of the Lie algebra::

        sage: R.<a,b,c> = FreeAlgebra(QQ, 3)
        sage: L.<x,y> = LieAlgebra(associative=[a,b])

    Now for a more complicated example using the group ring of `S_3`
    as our base algebra::

        sage: G = SymmetricGroup(3)
        sage: S = GroupAlgebra(G, QQ)
        sage: L.<x,y> = LieAlgebra(associative=S.gens())
        sage: L.bracket(x, y)
        (2,3) - (1,3)
        sage: L.bracket(x, y-x)
        (2,3) - (1,3)
        sage: L.bracket(L.bracket(x, y), y)
        2*(1,2,3) - 2*(1,3,2)
        sage: L.bracket(x, L.bracket(x, y))
        (2,3) - 2*(1,2) + (1,3)
        sage: L.bracket(x, L.bracket(L.bracket(x, y), y))
        0

    Here is an example using matrices::

        sage: MS = MatrixSpace(QQ,2)
        sage: m1 = MS([[0, -1], [1, 0]])
        sage: m2 = MS([[-1, 4], [3, 2]])
        sage: L.<x,y> = LieAlgebra(associative=[m1, m2])
        sage: x
        [ 0 -1]
        [ 1  0]
        sage: y
        [-1  4]
        [ 3  2]
        sage: L.bracket(x,y)
        [-7 -3]
        [-3  7]
        sage: L.bracket(y,y)
        [0 0]
        [0 0]
        sage: L.bracket(y,x)
        [ 7  3]
        [ 3 -7]
        sage: L.bracket(x, L.bracket(y,x))
        [-6 14]
        [14  6]
    """
    @staticmethod
    def __classcall_private__(cls, A, gens=None, names=None, index_set=None,
                              free_lie_algebra=False, category=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L1 = LieAlgebra(associative=tuple(S.gens()), names=['x','y'])
            sage: L2 = LieAlgebra(associative=[ S(G((1,2,3))), S(G((1,2))) ], names='x,y')
            sage: L1 is L2
            True

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: L1 = LieAlgebra(associative=F.algebra_generators(), names='x,y,z')
            sage: L2.<x,y,z> = LieAlgebra(associative=F.gens())
            sage: L1 is L2
            True
        """
        # If A is not a ring, then we treat it as a set of generators
        if isinstance(A, Parent) and A.category().is_subcategory(Rings()):
            if gens is None and index_set is None:
                # Use the indexing set of the basis
                try:
                    index_set = A.basis().keys()
                except (AttributeError, NotImplementedError):
                    pass
        else:
            gens = A
            A = None

        if index_set is None:
            # See if we can get an index set from the generators
            try:
                index_set = gens.keys()
            except (AttributeError, ValueError):
                pass

        ngens = None
        if isinstance(gens, AbstractFamily):
            if index_set is None and names is None:
                index_set = gens.keys()
            if gens.cardinality() < float('inf'):
                # TODO: This makes the generators of a finitely generated
                #    Lie algebra into an ordered list for uniqueness and then
                #    reconstructs the family. Instead create a key for the
                #    cache this way and then pass the family.
                try:
                    gens = tuple([gens[i] for i in index_set])
                except KeyError:
                    gens = tuple(gens)
                ngens = len(gens)

        elif isinstance(gens, dict):
            if index_set is None and names is None:
                index_set = gens.keys()
            gens = gens.values()
            ngens = len(gens)
        elif gens is not None:  # Assume it is list-like
            gens = tuple(gens)
            ngens = len(gens)
            if index_set is None and names is None:
                index_set = list(range(ngens))

        if ngens is not None:
            if A is None:
                A = gens[0].parent()
            # Make sure all the generators have the same parent of A
            gens = tuple([A(g) for g in gens])

        names, index_set = standardize_names_index_set(names, index_set, ngens)

        # We strip the following axioms from the category of the assoc. algebra:
        #   FiniteDimensional and WithBasis
        category = LieAlgebras(A.base_ring()).or_subcategory(category)
        if 'FiniteDimensional' in A.category().axioms():
            category = category.FiniteDimensional()
        if 'WithBasis' in A.category().axioms() and gens is None:
            category = category.WithBasis()

        if isinstance(A, MatrixSpace):
            if gens is not None:
                for g in gens:
                    g.set_immutable()
            return MatrixLieAlgebraFromAssociative(A, gens, names=names,
                                                   index_set=index_set,
                                                   category=category)

        return super().__classcall__(cls, A, gens, names=names,
                                     index_set=index_set, category=category)

    def __init__(self, A, gens=None, names=None, index_set=None, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: TestSuite(L).run()

        TESTS::

            sage: from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative as LAFA
            sage: LAFA(MatrixSpace(QQ, 0, sparse=True), [], names=())
            Lie algebra generated by () in Full MatrixSpace of 0 by 0 sparse matrices over Rational Field
        """
        self._assoc = A
        R = self._assoc.base_ring()

        LieAlgebraWithGenerators.__init__(self, R, names, index_set, category)

        if isinstance(gens, tuple):
            # This guarantees that the generators have a specified ordering
            d = {self._indices[i]: self.element_class(self, v)
                 for i, v in enumerate(gens)}
            gens = Family(self._indices, lambda i: d[i])
        elif gens is not None:  # It is a family
            gens = Family(self._indices,
                          lambda i: self.element_class(self, gens[i]),
                          name="generator map")
        self._gens = gens

        # We don't need to store the original generators because we can
        #   get them from lifting this object's generators

        # We construct the lift map to the ambient associative algebra
        LiftMorphismToAssociative(self, self._assoc).register_as_coercion()

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: MS = MatrixSpace(QQ,2)
            sage: L.<x> = LieAlgebra(associative=[MS.one()])
            sage: L._repr_option('element_ascii_art')
            True
        """
        return self._assoc._repr_option(key)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: LieAlgebra(associative=S)
            Lie algebra of Symmetric group algebra of order 3
            over Rational Field
            sage: LieAlgebra(associative=S.gens())
            Lie algebra generated by ((1,2,3), (1,2))
            in Symmetric group algebra of order 3 over Rational Field
        """
        if self._gens is not None:
            return "Lie algebra generated by {} in {}".format(tuple(self._gens), self._assoc)
        return "Lie algebra of {}".format(self._assoc)

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: L = LieAlgebra(associative=S)
            sage: x,y = S.algebra_generators()
            sage: elt = L(x - y); elt
            [2, 1, 3] - [2, 3, 1]
            sage: elt.parent() is L
            True
            sage: elt == L(x) - L(y)
            True
            sage: L([x, y])
            -[1, 3, 2] + [3, 2, 1]
            sage: L(2)
            2*[1, 2, 3]
        """
        if isinstance(x, list) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))
        return self.element_class(self, self._assoc(x))

    def associative_algebra(self):
        """
        Return the associative algebra used to construct ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.associative_algebra() is S
            True
        """
        return self._assoc

    def lie_algebra_generators(self):
        """
        Return the Lie algebra generators of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.lie_algebra_generators()
            Finite family {(): (), (1,3,2): (1,3,2), (1,2,3): (1,2,3),
                           (2,3): (2,3), (1,3): (1,3), (1,2): (1,2)}
        """
        if self._gens is not None:
            return self._gens
        try:
            ngens = self._indices.cardinality()
        except AttributeError:
            ngens = len(self._indices)
        if ngens < float('inf'):
            return Family(list(self._indices), self.monomial)
        return Family(self._indices, self.monomial, name="generator map")

    def monomial(self, i):
        """
        Return the monomial indexed by ``i``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ)
            sage: L = LieAlgebra(associative=F)
            sage: L.monomial(x.leading_support())
            x
        """
        if i not in self._assoc.basis().keys():
            # return self(self._assoc.monomial(i))
            raise ValueError("not an index")
        return self.element_class(self, self._assoc.monomial(i))

    def term(self, i, c=None):
        """
        Return the term indexed by ``i`` with coefficient ``c``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ)
            sage: L = LieAlgebra(associative=F)
            sage: L.term(x.leading_support(), 4)
            4*x
        """
        if i not in self._assoc.basis().keys():
            # return self(self._assoc.term(i, c))
            raise ValueError("not an index")
        return self.element_class(self, self._assoc.term(i, c))

    @cached_method
    def zero(self):
        """
        Return the element `0` in ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.zero()
            0
        """
        return self.element_class(self, self._assoc.zero())

    def is_abelian(self):
        """
        Return ``True`` if ``self`` is abelian.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 2, 'x,y')
            sage: L = LieAlgebra(associative=R.gens())
            sage: L.is_abelian()
            False

            sage: R = PolynomialRing(QQ, 'x,y')
            sage: L = LieAlgebra(associative=R.gens())
            sage: L.is_abelian()
            True

        An example with a Lie algebra from the group algebra::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.is_abelian()
            False

        Now we construct a Lie algebra from commuting elements in the group
        algebra::

            sage: G = SymmetricGroup(5)
            sage: S = GroupAlgebra(G, QQ)
            sage: gens = map(S, [G((1, 2)), G((3, 4))])
            sage: L.<x,y> = LieAlgebra(associative=gens)
            sage: L.is_abelian()
            True
        """
        if self._assoc.is_commutative():
            return True
        return super().is_abelian()

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ)

        An infinitely generated example::

            sage: L = LieAlgebra(associative=F)
            sage: L.an_element()
            1

        A finitely generated example::

            sage: L = LieAlgebra(associative=F.gens())
            sage: L.an_element()
            x + y
        """
        G = self.lie_algebra_generators()
        if G.cardinality() < float('inf'):
            return self.sum(G)
        return G[self._indices.an_element()]

    class Element(LieAlgebraElementWrapper):
        def _bracket_(self, rhs):
            """
            Return the bracket ``[self, rhs]``.

            EXAMPLES::

                sage: L.<x,y,z> = LieAlgebra(QQ, representation='polynomial')
                sage: L.bracket(x, y)
                x*y - y*x

                sage: G = SymmetricGroup(3)
                sage: S = GroupAlgebra(G, QQ)
                sage: L.<x,y> = LieAlgebra(associative=S.gens())
                sage: L.bracket(x, y)
                (2,3) - (1,3)

                sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
                sage: L.bracket(L.gen(0), L.gen(1))
                [ 1  0]
                [ 0 -1]
            """
            ret = self.value * rhs.value - rhs.value * self.value
            return self.__class__(self.parent(), ret)

        def lift_associative(self):
            """
            Lift ``self`` to the ambient associative algebra (which
            might be smaller than the universal enveloping algebra).

            EXAMPLES::

                sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
                sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
                sage: x.lift_associative()
                x
                sage: x.lift_associative().parent()
                Free Algebra on 3 generators (x, y, z) over Rational Field
            """
            return self.value

        def monomial_coefficients(self, copy=True):
            """
            Return the monomial coefficients of ``self`` (if this
            notion makes sense for ``self.parent()``).

            EXAMPLES::

                sage: R.<x,y,z> = FreeAlgebra(QQ)
                sage: L = LieAlgebra(associative=R)
                sage: elt = L(x) + 2*L(y) - L(z)
                sage: sorted(elt.monomial_coefficients().items())
                [(x, 1), (y, 2), (z, -1)]

                sage: L = LieAlgebra(associative=[x,y])
                sage: elt = L(x) + 2*L(y)
                sage: elt.monomial_coefficients()
                Traceback (most recent call last):
                ...
                NotImplementedError: the basis is not defined
            """
            if self.parent()._gens is not None:
                raise NotImplementedError("the basis is not defined")
            return self.value.monomial_coefficients(copy)


class LiftMorphismToAssociative(LiftMorphism):
    """
    The natural lifting morphism from a Lie algebra constructed from
    an associative algebra `A` to `A`.
    """
    def preimage(self, x):
        """
        Return the preimage of ``x`` under ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'a,b,c')
            sage: L = LieAlgebra(associative=R)
            sage: x,y,z = R.gens()
            sage: f = R.coerce_map_from(L)
            sage: p = f.preimage(x*y - z); p
            -c + a*b
            sage: p.parent() is L
            True
        """
        return self.domain().element_class(self.domain(), x)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: f = R.coerce_map_from(L)
            sage: a = f(L([x,y]) + z); a
            z + x*y - y*x
            sage: a.parent() is R
            True
        """
        return x.value

    def section(self):
        """
        Return the section map of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: f = R.coerce_map_from(L)
            sage: f.section()
            Generic morphism:
              From: Free Algebra on 3 generators (x, y, z) over Rational Field
              To:   Lie algebra generated by (x, y, z) in Free Algebra on 3 generators (x, y, z) over Rational Field
        """
        return SetMorphism(Hom(self.codomain(), self.domain()),
                           self.preimage)


class MatrixLieAlgebraFromAssociative(LieAlgebraFromAssociative):
    """
    A Lie algebra constructed from a matrix algebra.

    This means a Lie algebra consisting of matrices,
    with commutator as Lie bracket.
    """
    class Element(LieAlgebraMatrixWrapper, LieAlgebraFromAssociative.Element):
        def matrix(self):
            r"""
            Return ``self`` as element of the underlying matrix algebra.

            OUTPUT: an instance of the element class of MatrixSpace

            EXAMPLES::

                sage: sl3m = lie_algebras.sl(ZZ, 3, representation='matrix')
                sage: e1,e2, f1, f2, h1, h2 = sl3m.gens()
                sage: h1m = h1.matrix(); h1m
                [ 1  0  0]
                [ 0 -1  0]
                [ 0  0  0]
                sage: h1m.parent()
                Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
                sage: matrix(h2)
                [ 0  0  0]
                [ 0  1  0]
                [ 0  0 -1]
                sage: L = lie_algebras.so(QQ['z'], 5, representation='matrix')
                sage: matrix(L.an_element())
                [ 1  1  0  0  0]
                [ 1  1  0  0  2]
                [ 0  0 -1 -1  0]
                [ 0  0 -1 -1 -1]
                [ 0  1  0 -2  0]

                sage: gl2 = lie_algebras.gl(QQ, 2)
                sage: matrix(gl2.an_element())
                [1 1]
                [1 1]
            """
            return self.value

        _matrix_ = matrix
