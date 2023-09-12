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

from sage.sets.family import Family
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
      elements of ``lie_algebra`` encoded as a ``dict`` with keys being
      indices of the basis of ``lie_algebra`` and the values being the
      corresponding matrix defining the action

    EXAMPLES::

        sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
        sage: f = ({x: Matrix([[1,0],[0,0]]), y: Matrix([[0,1],[0,0]])})
        sage: L.representation(f)
        Representation of Lie algebra on 2 generators (x, y) over Rational Field defined by:
               [1 0]
        x |--> [0 0]
               [0 1]
        y |--> [0 0]
    """
    @staticmethod
    def __classcall_private__(cls, lie_algebra, f, **kwargs):
        r"""
        Normalize inpute to ensure a unique representation.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f1 = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: R1 = L.representation(f1)
            sage: f2 = Family({x: Matrix([[1,0],[0,0]]), y: Matrix([[0,1],[0,0]])})
            sage: R2 = L.representation(f2)
            sage: R1 is R2
            True
        """
        C = Modules(lie_algebra.base_ring()).WithBasis().FiniteDimensional()
        C = C.or_subcategory(kwargs.pop('category', C))
        B = lie_algebra.basis()
        data = {}
        for k, mat in f.items():
            if k in B:
                k = k.leading_support()
            data[k] = copy(mat)
            data[k].set_immutable()
        f = Family(data)
        return super(cls, RepresentationByMorphism).__classcall__(cls, lie_algebra, f, category=C, **kwargs)

    def __init__(self, lie_algebra, f, category, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: R = L.representation(f)
            sage: TestSuite(R).run()
        """
        it = iter(f)
        mat = next(it)
        if not mat.is_square():
            raise ValueError("all matrices must be square")
        dim = mat.nrows()
        self._mat_space = mat.parent()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        if any(mat.nrows() != dim or mat.ncols() != dim for mat in it):
            raise ValueError("all matrices must be square of size {}".format(dim))
        self._f = dict(f)

        I = FiniteEnumeratedSet(range(dim))
        Representation_abstract.__init__(self, lie_algebra)
        CombinatorialFreeModule.__init__(self, lie_algebra.base_ring(), I, prefix='R', category=category)

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
        """
        ret = "Representation of {} defined by:".format(self._lie_algebra)
        from sage.typeset.ascii_art import ascii_art
        B = self._lie_algebra.basis()
        for k in self._f:
            ret += '\n' + repr(ascii_art(B[k], self._f[k], sep=" |--> ", sep_baseline=0))
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
            """
            P = self.parent()
            if scalar in P._lie_algebra:
                if self_on_left:
                    return None
                scalar = P._lie_algebra(scalar)
                f = P._f
                mat = P._mat_space.sum(scalar[k] * f[k] for k in f if scalar[k])
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
    def __init__(self, lie_algebra):
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
        CombinatorialFreeModule.__init__(self, R, ['v'], prefix='T', category=cat)

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
