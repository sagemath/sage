r"""
Representations of Lie algebras

AUTHORS:

- Travis Scrimshaw (2023-08-31): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.family import Family, AbstractFamily
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.modules import Modules
from copy import copy


class Representation_abstract:
    r"""
    Mixin class for (left) representations of Lie algebras.

    INPUT:

    - ``lie_algebra`` -- a Lie algebra
    """
    def __init__(self, lie_algebra):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: R = L.trivial_representation()
            sage: TestSuite(R).run()
        """
        self._lie_algebra = lie_algebra

    def lie_algebra(self):
        r"""
        Return the Lie algebra whose representation ``self`` is.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 4)
            sage: R = L.trivial_representation()
            sage: R.lie_algebra() is L
            True
        """
        return self._lie_algebra

    def side(self):
        r"""
        Return that ``self`` is a left representation.

        OUTPUT:

        - the string ``"left"``

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 4)
            sage: R = L.trivial_representation()
            sage: R.side()
            'left'
        """
        return 'left'

    def _test_representation(self, **options):
        r"""
        Check (on some elements) that ``self`` is a representation of the
        given Lie algebra using the basis of the Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 3)
            sage: f = {b: b.adjoint_matrix() for b in L.basis()}
            sage: R = L.representation(f)
            sage: R._test_representation()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        elts = self._lie_algebra.basis()
        if elts.cardinality() == float('inf'):
            elts = list(elts.some_elements())
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(elts, 2, tester._max_runs):
            for v in S:
                tester.assertEqual(x.bracket(y) * v, x * (y * v) - y * (x * v))


class RepresentationByMorphism(CombinatorialFreeModule, Representation_abstract):
    r"""
    Representation of a Lie algebra defined by a Lie algebra morphism.

    INPUT:

    - ``lie_algebra`` -- a Lie algebra
    - ``f`` -- the Lie algebra morphism defining the action of the basis
      elements of ``lie_algebra``
    - ``index_set`` -- (optional) the index set of the module basis
    - ``on_basis`` -- (default: ``False``) the function ``f`` defines a
      map from the basis elements or from a generic element of ``lie_algebra``

    If ``f`` is encoded as a ``dict`` or ``Family``, then the keys must
    be indices of the basis of ``lie_algebra`` and the values being the
    corresponding matrix defining the action. This sets ``on_basis=True``.

    EXAMPLES::

        sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
        sage: f = {x: Matrix([[1,0],[0,0]]), y: Matrix([[0,1],[0,0]])}
        sage: L.representation(f)
        Representation of Lie algebra on 2 generators (x, y) over Rational Field defined by:
               [1 0]
        x |--> [0 0]
               [0 1]
        y |--> [0 0]

    We construct the direct sum of two copies of the trivial representation
    for an infinite dimensional Lie algebra::

        sage: L = lie_algebras.Affine(QQ, ['E',6,1])
        sage: R = L.representation(lambda b: matrix.zero(QQ, 2), index_set=['a','b'])
        sage: x = L.an_element()
        sage: v = R.an_element(); v
        2*R['a'] + 2*R['b']
        sage: x * v
        0

    We construct a finite dimensional representation of the affline Lie algebra
    of type `A_2^{(1)}`::

        sage: L = lie_algebras.Affine(QQ, ['A',2,1]).derived_subalgebra()
        sage: Phi_plus = list(RootSystem(['A',2]).root_lattice().positive_roots())
        sage: def aff_action(key):
        ....:     mat = matrix.zero(QQ, 3)
        ....:     if key == 'c':  # central element
        ....:         return mat
        ....:     b, ell = key
        ....:     if b in Phi_plus:  # positive root
        ....:         ind = tuple(sorted(b.to_ambient().support()))
        ....:         mat[ind] = 1
        ....:         if ind[0] + 1 != ind[1]:
        ....:             mat[ind] = -1
        ....:     elif -b in Phi_plus:  # negative root
        ....:         ind = tuple(sorted(b.to_ambient().support(), reverse=True))
        ....:         mat[ind] = 1
        ....:         if ind[0] - 1 != ind[1]:
        ....:             mat[ind] = -1
        ....:     else:  # must be in the Cartan
        ....:         i = b.leading_support()
        ....:         mat[i,i] = -1
        ....:         mat[i-1,i-1] = 1
        ....:     return mat
        sage: F = Family(L.basis(), aff_action, name="lifted natural repr")
        sage: R = L.representation(index_set=range(1,4), on_basis=F)
        sage: x = L.an_element(); x
        (E[alpha[2]] + E[alpha[1]] + h1 + h2 + E[-alpha[2]] + E[-alpha[1]])#t^0
         + (E[-alpha[1] - alpha[2]])#t^1 + (E[alpha[1] + alpha[2]])#t^-1 + c
        sage: v = R.an_element(); v
        2*R[1] + 2*R[2] + 3*R[3]
        sage: x * v
        R[1] + 5*R[2] - 3*R[3]
        sage: R._test_representation()  # verify that it is a representation
    """
    @staticmethod
    def __classcall_private__(cls, lie_algebra, f=None, index_set=None, on_basis=False, **kwargs):
        r"""
        Normalize inpute to ensure a unique representation.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f1 = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: R1 = L.representation(f1)
            sage: f2 = Family({x: Matrix([[1,0],[0,0]]), y: Matrix(QQ, [[0,1],[0,0]])})
            sage: R2 = L.representation(f2)
            sage: R1 is R2
            True

        TESTS::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0]]), 'y': Matrix([[0,1]])}
            sage: L.representation(f)
            Traceback (most recent call last):
            ...
            ValueError: all matrices must be square

            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0]])}
            sage: L.representation(f)
            Traceback (most recent call last):
            ...
            ValueError: all matrices must be square of size 2

            sage: L.representation(index_set=[1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: either 'f' or 'on_basis' must be specified
            sage: L.representation(on_basis=lambda x: QQ.zero())
            Traceback (most recent call last):
            ...
            ValueError: the index set needs to be specified
        """
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        base = lie_algebra.base_ring()
        C = Modules(base).WithBasis().FiniteDimensional()
        C = C.or_subcategory(kwargs.pop('category', C))
        B = lie_algebra.basis()
        if not isinstance(on_basis, bool):
            f = on_basis
            on_basis = True
        if isinstance(f, AbstractFamily):
            if f.cardinality() < float('inf'):
                f = dict(f)
            on_basis = True
        if isinstance(f, dict):
            data = {}
            dim = None
            for k, mat in f.items():
                if k in B:
                    k = k.leading_support()
                if not mat.is_square():
                    raise ValueError("all matrices must be square")
                if dim is None:
                    dim = mat.nrows()
                elif mat.nrows() != dim or mat.ncols() != dim:
                    raise ValueError("all matrices must be square of size {}".format(dim))
                data[k] = mat.change_ring(base)
                data[k].set_immutable()

            if index_set is None:
                index_set = FiniteEnumeratedSet(range(dim))
            f = Family(data)
            on_basis = True

        if f is None:
            raise ValueError("either 'f' or 'on_basis' must be specified")
        if index_set is None:
            raise ValueError("the index set needs to be specified")

        index_set = FiniteEnumeratedSet(index_set)

        return super(cls, RepresentationByMorphism).__classcall__(cls, lie_algebra,
             f, index_set, on_basis, category=C, **kwargs)

    def __init__(self, lie_algebra, f, index_set, on_basis, category, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: R = L.representation(f)
            sage: TestSuite(R).run()
        """
        if on_basis:
            self._family = f
            self._f = f.__getitem__
        else:
            self._f = f
        prefix = kwargs.pop("prefix", 'R')
        self._on_basis = on_basis

        Representation_abstract.__init__(self, lie_algebra)
        CombinatorialFreeModule.__init__(self, lie_algebra.base_ring(), index_set,
                                         category=category, prefix=prefix, **kwargs)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: L.representation(f)
            Representation of Lie algebra on 2 generators (x, y) over Rational Field defined by:
                   [1 0]
            x |--> [0 0]
                   [0 1]
            y |--> [0 0]

            sage: L = lie_algebras.Affine(QQ, ['E',6,1])
            sage: F = Family(L.basis(), lambda b: matrix.zero(QQ, 2), name="zero map")
            sage: L.representation(F, index_set=['a','b'], on_basis=True)
            Representation of Affine Kac-Moody algebra of ['E', 6] in the Chevalley basis defined by:
            Lazy family (zero map(i))_{i in Lazy family...}

            sage: L.representation(lambda b: matrix.zero(QQ, 2), index_set=['a','b'])
            Representation of Affine Kac-Moody algebra of ['E', 6] in the Chevalley basis defined by:
            <function <lambda> at 0x...>
        """
        ret = "Representation of {} defined by:".format(self._lie_algebra)
        from sage.typeset.ascii_art import ascii_art
        if self._on_basis:
            B = self._lie_algebra.basis()
            if B.cardinality() < float('inf'):
                for k in B.keys():
                    ret += '\n' + repr(ascii_art(B[k], self._f(k), sep=" |--> ", sep_baseline=0))
            else:
                ret += '\n' + repr(self._family)
        else:
            ret += '\n' + repr(self._f)
        return ret

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
                sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
                sage: R = L.representation(f)
                sage: v = R.an_element(); v
                2*R[0] + 2*R[1]
                sage: x * v
                2*R[0]
                sage: y * v
                2*R[0]
                sage: (2*x + 5*y) * v
                14*R[0]
                sage: v * x
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: ...

                sage: v = sum((i+4) * b for i, b in enumerate(R.basis())); v
                4*R[0] + 5*R[1]
                sage: (1/3*x - 5*y) * v
                -71/3*R[0]

                sage: L = lie_algebras.Affine(QQ, ['E',6,1])
                sage: F = Family(L.basis(), lambda b: matrix.zero(QQ, 2), name="zero map")
                sage: R = L.representation(F, index_set=['a','b'], on_basis=True)
                sage: R.an_element()
                2*R['a'] + 2*R['b']
                sage: L.an_element() * R.an_element()
                0
            """
            P = self.parent()
            if scalar in P._lie_algebra:
                if self_on_left:
                    return None
                if not self:  # we are (already) the zero vector
                    return self
                scalar = P._lie_algebra(scalar)
                if not scalar:  # we are acting by zero
                    return P.zero()
                if P._on_basis:
                    mat = sum(c * P._f(k) for k, c in scalar.monomial_coefficients(copy=False).items())
                else:
                    mat = P._f(scalar)
                return P.from_vector(mat * self.to_vector())

            return super()._acted_upon_(scalar, self_on_left)


class TrivialRepresentation(CombinatorialFreeModule, Representation_abstract):
    r"""
    The trivial representation of a Lie algebra.

    The trivial representation of a Lie algebra `L` over a commutative ring
    `R` is the `1`-dimensional `R`-module on which every element of `L`
    acts by zero.

    INPUT:

    - ``lie_algebra`` -- a Lie algebra

    REFERENCES:

    - :wikipedia:`Trivial_representation`
    """
    def __init__(self, lie_algebra, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: R = L.trivial_representation()
            sage: TestSuite(R).run()
        """
        R = lie_algebra.base_ring()
        cat = Modules(R).WithBasis().FiniteDimensional()
        Representation_abstract.__init__(self, lie_algebra)
        CombinatorialFreeModule.__init__(self, R, ['v'], prefix='T', category=cat, **kwargs)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: L.trivial_representation()
            Trivial representation of The Virasoro algebra over Rational Field
        """
        return "Trivial representation of {}".format(self._lie_algebra)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: L = lie_algebras.VirasoroAlgebra(QQ)
                sage: R = L.trivial_representation()
                sage: L.an_element() * R.an_element()
                0
                sage: R.an_element() * L.an_element()
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: ...
                sage: 3 / 5 * R.an_element()
                6/5*T['v']
            """
            P = self.parent()
            if scalar in P._lie_algebra:
                if self_on_left:
                    return None
                return P.zero()
            return super()._acted_upon_(scalar, self_on_left)
