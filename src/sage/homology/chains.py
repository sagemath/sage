# -*- coding: utf-8 -*-
r"""
Chains and cochains

This module implements formal linear combinations of cells of a given
cell complex (:class:`Chains`) and their dual (:class:`Cochains`). It
is closely related to the :mod:`sage.topology.chain_complex`
module. The main differences are that chains and cochains here are of
homogeneous dimension only, and that they reference their cell
complex.
"""

#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.structure.element import coercion_model


class CellComplexReference(object):

    def __init__(self, cell_complex, degree, cells=None):
        """
        Auxiliary base class for chains and cochains

        INPUT:

        - ``cell_complex`` -- The cell complex to reference

        - ``degree`` -- integer. The degree of the (co)chains

        - ``cells`` -- tuple of cells or ``None``. Does not necessarily have to
          be the cells in the given degree, for computational purposes this
          could also be any collection that is in one-to-one correspondence with
          the cells. If ``None``, the cells of the complex in the given degree
          are used.

        EXAMPLES::

            sage: X = simplicial_complexes.Simplex(2)
            sage: from sage.homology.chains import CellComplexReference
            sage: c = CellComplexReference(X, 1)
            sage: c.cell_complex() is X
            True
        """
        self._cell_complex = cell_complex
        self._degree = degree
        if cells is not None:
            self._cells = cells
        else:
            self._cells = tuple(cell_complex.n_cells(degree))

    def cell_complex(self):
        """
        Return the underlying cell complex

        OUTPUT:

        A cell complex.

        EXAMPLES::

            sage: X = simplicial_complexes.Simplex(2)
            sage: X.n_chains(1).cell_complex() is X
            True
        """
        return self._cell_complex

    def degree(self):
        """
        Return the dimension of the cells

        OUTPUT:

        Integer. The dimension of the cells.

        EXAMPLES::

            sage: X = simplicial_complexes.Simplex(2)
            sage: X.n_chains(1).degree()
            1
        """
        return self._degree


class Chains(CellComplexReference, CombinatorialFreeModule):
    r"""
    Class for the free module of chains in a given degree.

    INPUT:

    - ``n_cells`` -- tuple of `n`-cells, which thus forms a basis for
      this module
    - ``base_ring`` -- optional (default `\ZZ`)

    One difference between chains and cochains is notation. In a
    simplicial complex, for example, a simplex ``(0,1,2)`` is written
    as "(0,1,2)" in the group of chains but as "\\chi_(0,1,2)" in the
    group of cochains.

    Also, since the free modules of chains and cochains are dual,
    there is a pairing `\langle c, z \rangle`, sending a cochain `c`
    and a chain `z` to a scalar.

    EXAMPLES::

        sage: S2 = simplicial_complexes.Sphere(2)
        sage: C_2 = S2.n_chains(1)
        sage: C_2_co = S2.n_chains(1, cochains=True)
        sage: x = C_2.basis()[Simplex((0,2))]
        sage: y = C_2.basis()[Simplex((1,3))]
        sage: z = x+2*y
        sage: a = C_2_co.basis()[Simplex((1,3))]
        sage: b = C_2_co.basis()[Simplex((0,3))]
        sage: c = 3*a-2*b
        sage: z
        (0, 2) + 2*(1, 3)
        sage: c
        -2*\chi_(0, 3) + 3*\chi_(1, 3)
        sage: c.eval(z)
        6
    """
    def __init__(self, cell_complex, degree, cells=None, base_ring=None):
        """
        EXAMPLES::

            sage: T = cubical_complexes.Torus()
            sage: C = T.n_chains(2, QQ)
            sage: C.dimension()
            16
            sage: sorted(C.gens())
            [[0,0] x [0,1] x [0,0] x [0,1],
             [0,0] x [0,1] x [0,1] x [0,0],
             [0,0] x [0,1] x [0,1] x [1,1],
             [0,0] x [0,1] x [1,1] x [0,1],
             [0,1] x [0,0] x [0,0] x [0,1],
             [0,1] x [0,0] x [0,1] x [0,0],
             [0,1] x [0,0] x [0,1] x [1,1],
             [0,1] x [0,0] x [1,1] x [0,1],
             [0,1] x [1,1] x [0,0] x [0,1],
             [0,1] x [1,1] x [0,1] x [0,0],
             [0,1] x [1,1] x [0,1] x [1,1],
             [0,1] x [1,1] x [1,1] x [0,1],
             [1,1] x [0,1] x [0,0] x [0,1],
             [1,1] x [0,1] x [0,1] x [0,0],
             [1,1] x [0,1] x [0,1] x [1,1],
             [1,1] x [0,1] x [1,1] x [0,1]]

        TESTS::

            sage: T.n_chains(2).base_ring()
            Integer Ring
            sage: T.n_chains(8).dimension()
            0
            sage: T.n_chains(-3).dimension()
            0
        """
        if base_ring is None:
            base_ring = ZZ
        CellComplexReference.__init__(self, cell_complex, degree, cells=cells)
        CombinatorialFreeModule.__init__(
            self, base_ring, self._cells,
            prefix='', bracket=False
        )

    def dual(self):
        """
        Return the cochains.

        OUTPUT:

        The cochains of the same cells with the same base ring.

        EXAMPLES::

            sage: square = cubical_complexes.Cube(2)
            sage: chains = square.n_chains(1, ZZ);  chains
            Free module generated by {[0,0] x [0,1], [0,1] x [0,0], [0,1] x [1,1], [1,1] x [0,1]} over Integer Ring
            sage: chains.dual()
            Free module generated by {[0,0] x [0,1], [0,1] x [0,0], [0,1] x [1,1], [1,1] x [0,1]} over Integer Ring
            sage: type(chains)
            <class 'sage.homology.chains.Chains_with_category'>
            sage: type(chains.dual())
            <class 'sage.homology.chains.Cochains_with_category'>
        """
        return Cochains(
            self.cell_complex, self.degree,
            cells=self._cells, base_ring=self.base_ring()
        )

    def chain_complex(self):
        """
        Return the chain complex.

        OUTPUT:

        Chain complex, see :mod:`sage.homology.chain_complex`.

        EXAMPLES::

            sage: square = cubical_complexes.Cube(2)
            sage: CC = square.n_chains(2, QQ).chain_complex(); CC
            Chain complex with at most 3 nonzero terms over Rational Field
            sage: ascii_art(CC)
                        [-1 -1  0  0]       [-1]
                        [ 1  0 -1  0]       [ 1]
                        [ 0  1  0 -1]       [-1]
                        [ 0  0  1  1]       [ 1]
             0 <-- C_0 <-------------- C_1 <----- C_2 <-- 0
        """
        return self.cell_complex().chain_complex(
            base_ring=self.base_ring(),
            cochain=False,
        )

    class Element(CombinatorialFreeModule.Element):

        def to_complex(self):
            """
            Return the corresponding chain complex element

            OUTPUT:

            An element of the chain complex, see :mod:`sage.homology.chain_complex`.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ)
                sage: from sage.topology.cubical_complex import Cube
                sage: chain = C1(Cube([[1, 1], [0, 1]]))
                sage: chain.to_complex()
                Chain(1:(0, 0, 0, 1))
                sage: ascii_art(_)
                   d_0  [0]  d_1  [0]  d_2       d_3
                0 <---- [0] <---- [0] <---- [0] <---- 0
                        [0]       [0]
                        [0]       [1]
            """
            return self.parent().chain_complex()({
                self.parent().degree(): self.to_vector()
            })

        def boundary(self):
            """
            Return the boundary of the chain

            OUTPUT:

            The boundary as a chain in one degree lower.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ)
                sage: from sage.topology.cubical_complex import Cube
                sage: chain = C1(Cube([[1, 1], [0, 1]])) - 2 * C1(Cube([[0, 1], [0, 0]]))
                sage: chain
                -2*[0,1] x [0,0] + [1,1] x [0,1]
                sage: chain.boundary()
                2*[0,0] x [0,0] - 3*[1,1] x [0,0] + [1,1] x [1,1]
            """
            chains = self.parent()
            degree = chains.degree()
            d = chains.chain_complex().differential(degree)
            codomain = chains.cell_complex().n_chains(
                degree - 1, base_ring=self.base_ring()
            )
            return codomain.from_vector(d * self.to_vector())

        def is_cycle(self):
            """
            Test whether the chain is a cycle

            OUTPUT:

            Boolean. Whether the :meth:`boundary` vanishes.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ)
                sage: from sage.topology.cubical_complex import Cube
                sage: chain = C1(Cube([[1, 1], [0, 1]])) - C1(Cube([[0, 1], [0, 0]]))
                sage: chain.is_cycle()
                False

            TESTS::

                sage: chain.is_cocycle()
                Traceback (most recent call last):
                ...
                AttributeError: 'Chains_with_category.element_class' object has no attribute 'is_cocycle'
            """
            return self.to_complex().is_cycle()

        def is_boundary(self):
            """
            Test whether the chain is a boundary

            OUTPUT:

            Boolean. Whether the chain is the :meth:`boundary` of a chain in one
            degree higher.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ)
                sage: from sage.topology.cubical_complex import Cube
                sage: chain = C1(Cube([[1, 1], [0, 1]])) - C1(Cube([[0, 1], [0, 0]]))
                sage: chain.is_boundary()
                False

            TESTS::

                sage: chain.is_coboundary()
                Traceback (most recent call last):
                ...
                AttributeError: 'Chains_with_category.element_class' object has no attribute 'is_coboundary'
            """
            return self.to_complex().is_boundary()


class Cochains(CellComplexReference, CombinatorialFreeModule):
    r"""
    Class for the free module of cochains in a given degree.

    INPUT:

    - ``n_cells`` -- tuple of `n`-cells, which thus forms a basis for
      this module
    - ``base_ring`` -- optional (default `\ZZ`)

    One difference between chains and cochains is notation. In a
    simplicial complex, for example, a simplex ``(0,1,2)`` is written
    as "(0,1,2)" in the group of chains but as "\\chi_(0,1,2)" in the
    group of cochains.

    Also, since the free modules of chains and cochains are dual,
    there is a pairing `\langle c, z \rangle`, sending a cochain `c`
    and a chain `z` to a scalar.

    EXAMPLES::

        sage: S2 = simplicial_complexes.Sphere(2)
        sage: C_2 = S2.n_chains(1)
        sage: C_2_co = S2.n_chains(1, cochains=True)
        sage: x = C_2.basis()[Simplex((0,2))]
        sage: y = C_2.basis()[Simplex((1,3))]
        sage: z = x+2*y
        sage: a = C_2_co.basis()[Simplex((1,3))]
        sage: b = C_2_co.basis()[Simplex((0,3))]
        sage: c = 3*a-2*b
        sage: z
        (0, 2) + 2*(1, 3)
        sage: c
        -2*\chi_(0, 3) + 3*\chi_(1, 3)
        sage: c.eval(z)
        6
    """
    def __init__(self, cell_complex, degree, cells=None, base_ring=None):
        """
        EXAMPLES::

            sage: T = cubical_complexes.Torus()
            sage: C = T.n_chains(2, QQ)
            sage: C.dimension()
            16
            sage: sorted(C.gens())
            [[0,0] x [0,1] x [0,0] x [0,1],
             [0,0] x [0,1] x [0,1] x [0,0],
             [0,0] x [0,1] x [0,1] x [1,1],
             [0,0] x [0,1] x [1,1] x [0,1],
             [0,1] x [0,0] x [0,0] x [0,1],
             [0,1] x [0,0] x [0,1] x [0,0],
             [0,1] x [0,0] x [0,1] x [1,1],
             [0,1] x [0,0] x [1,1] x [0,1],
             [0,1] x [1,1] x [0,0] x [0,1],
             [0,1] x [1,1] x [0,1] x [0,0],
             [0,1] x [1,1] x [0,1] x [1,1],
             [0,1] x [1,1] x [1,1] x [0,1],
             [1,1] x [0,1] x [0,0] x [0,1],
             [1,1] x [0,1] x [0,1] x [0,0],
             [1,1] x [0,1] x [0,1] x [1,1],
             [1,1] x [0,1] x [1,1] x [0,1]]

        TESTS::

            sage: T.n_chains(2, cochains=True).base_ring()
            Integer Ring
            sage: T.n_chains(8, cochains=True).dimension()
            0
            sage: T.n_chains(-3, cochains=True).dimension()
            0
        """
        if base_ring is None:
            base_ring = ZZ
        CellComplexReference.__init__(self, cell_complex, degree, cells=cells)
        CombinatorialFreeModule.__init__(
            self, base_ring, self._cells,
            prefix='\\chi',
            bracket=['_', '']
        )

    def dual(self):
        """
        Return the chains

        OUTPUT:

        The chains of the same cells with the same base ring.

        EXAMPLES::

            sage: square = cubical_complexes.Cube(2)
            sage: cochains = square.n_chains(1, ZZ, cochains=True);  cochains
            Free module generated by {[0,0] x [0,1], [0,1] x [0,0], [0,1] x [1,1], [1,1] x [0,1]} over Integer Ring
            sage: cochains.dual()
            Free module generated by {[0,0] x [0,1], [0,1] x [0,0], [0,1] x [1,1], [1,1] x [0,1]} over Integer Ring
            sage: type(cochains)
            <class 'sage.homology.chains.Cochains_with_category'>
            sage: type(cochains.dual())
            <class 'sage.homology.chains.Chains_with_category'>
        """
        return Chains(
            self.cell_complex, self.degree,
            cells=self._cells, base_ring=self.base_ring()
        )

    def cochain_complex(self):
        """
        Return the cochain complex.

        OUTPUT:

        Cochain complex, see :mod:`sage.homology.chain_complex`.

        EXAMPLES::

            sage: square = cubical_complexes.Cube(2)
            sage: C2 = square.n_chains(2, QQ, cochains=True)
            sage: C2.cochain_complex()
            Chain complex with at most 3 nonzero terms over Rational Field
            sage: ascii_art(C2.cochain_complex())
                                            [-1  1  0  0]
                                            [-1  0  1  0]
                                            [ 0 -1  0  1]
                        [-1  1 -1  1]       [ 0  0 -1  1]
             0 <-- C_2 <-------------- C_1 <-------------- C_0 <-- 0
        """
        return self.cell_complex().chain_complex(
            base_ring=self.base_ring(),
            cochain=True,
        )

    class Element(CombinatorialFreeModule.Element):

        def to_complex(self):
            """
            Return the corresponding cochain complex element

            OUTPUT:

            An element of the cochain complex, see :mod:`sage.homology.chain_complex`.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ, cochains=True)
                sage: from sage.topology.cubical_complex import Cube
                sage: cochain = C1(Cube([[1, 1], [0, 1]]))
                sage: cochain.to_complex()
                Chain(1:(0, 0, 0, 1))
                sage: ascii_art(_)
                   d_2       d_1  [0]  d_0  [0]  d_-1
                0 <---- [0] <---- [0] <---- [0] <----- 0
                                  [0]       [0]
                                  [1]       [0]
            """
            return self.parent().cochain_complex()({
                self.parent().degree(): self.to_vector()
            })

        def coboundary(self):
            r"""
            Return the coboundary of this cochain

            OUTPUT:

            The coboundary as a cochain in one degree higher.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ, cochains=True)
                sage: from sage.topology.cubical_complex import Cube
                sage: cochain = C1(Cube([[1, 1], [0, 1]])) - 2 * C1(Cube([[0, 1], [0, 0]]))
                sage: cochain
                -2*\chi_[0,1] x [0,0] + \chi_[1,1] x [0,1]
                sage: cochain.coboundary()
                -\chi_[0,1] x [0,1]
            """
            cochains = self.parent()
            degree = cochains.degree()
            d = cochains.cochain_complex().differential(degree)
            codomain = cochains.cell_complex().n_chains(
                degree + 1, base_ring=self.base_ring(), cochains=True
            )
            return codomain.from_vector(d * self.to_vector())

        def is_cocycle(self):
            """
            Test whether the cochain is a cocycle

            OUTPUT:

            Boolean. Whether the :meth:`coboundary` vanishes.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ, cochains=True)
                sage: from sage.topology.cubical_complex import Cube
                sage: cochain = C1(Cube([[1, 1], [0, 1]])) - C1(Cube([[0, 1], [0, 0]]))
                sage: cochain.is_cocycle()
                True

            TESTS::

                sage: cochain.is_cycle()
                Traceback (most recent call last):
                ...
                AttributeError: 'Cochains_with_category.element_class' object has no attribute 'is_cycle'
            """
            return self.to_complex().is_cycle()

        def is_coboundary(self):
            """
            Test whether the cochain is a coboundary

            OUTPUT:

            Boolean. Whether the cochain is the :meth:`coboundary` of a cochain
            in one degree lower.

            EXAMPLES::

                sage: square = cubical_complexes.Cube(2)
                sage: C1 = square.n_chains(1, QQ, cochains=True)
                sage: from sage.topology.cubical_complex import Cube
                sage: cochain = C1(Cube([[1, 1], [0, 1]])) - C1(Cube([[0, 1], [0, 0]]))
                sage: cochain.is_coboundary()
                True

            TESTS::

                sage: cochain.is_boundary()
                Traceback (most recent call last):
                ...
                AttributeError: 'Cochains_with_category.element_class' object has no attribute 'is_boundary'
            """
            return self.to_complex().is_boundary()

        def eval(self, other):
            r"""
            Evaluate this cochain on the chain ``other``.

            INPUT:

            - ``other`` -- a chain for the same cell complex in the
              same dimension with the same base ring

            OUTPUT: scalar

            EXAMPLES::

                sage: S2 = simplicial_complexes.Sphere(2)
                sage: C_2 = S2.n_chains(1)
                sage: C_2_co = S2.n_chains(1, cochains=True)
                sage: x = C_2.basis()[Simplex((0,2))]
                sage: y = C_2.basis()[Simplex((1,3))]
                sage: z = x+2*y
                sage: a = C_2_co.basis()[Simplex((1,3))]
                sage: b = C_2_co.basis()[Simplex((0,3))]
                sage: c = 3*a-2*b
                sage: z
                (0, 2) + 2*(1, 3)
                sage: c
                -2*\chi_(0, 3) + 3*\chi_(1, 3)
                sage: c.eval(z)
                6

            TESTS::

                sage: z.eval(c) # z is not a cochain
                Traceback (most recent call last):
                ...
                AttributeError: 'Chains_with_category.element_class' object has no attribute 'eval'
                sage: c.eval(c) # can't evaluate a cochain on a cochain
                Traceback (most recent call last):
                ...
                ValueError: argument is not a chain
            """
            if not isinstance(other.parent(), Chains):
                raise ValueError('argument is not a chain')
            if other.parent().indices() != self.parent().indices():
                raise ValueError('the cells are not compatible')
            result = sum(coeff * other.coefficient(cell)
                         for cell, coeff in self)
            R = self.base_ring()
            if R != other.base_ring():
                R = coercion_model.common_parent(R, other.base_ring())
            return R(result)

        def cup_product(self, cochain):
            r"""
            Return the cup product with another cochain.

            INPUT:

            - ``cochain`` -- cochain over the same cell complex

            EXAMPLES::

                sage: T2 = simplicial_complexes.Torus()
                sage: C1 = T2.n_chains(1, base_ring=ZZ, cochains=True)
                sage: def l(i, j):
                ....:      return C1(Simplex([i, j]))
                sage: l1 = l(1, 3) + l(1, 4) + l(1, 6) + l(2, 4) - l(4, 5) + l(5, 6)
                sage: l2 = l(1, 6) - l(2, 3) - l(2, 5) + l(3, 6) - l(4, 5) + l(5, 6)

            The two one-cocycles are cohomology generators::

                sage: l1.is_cocycle(), l1.is_coboundary()
                (True, False)
                sage: l2.is_cocycle(), l2.is_coboundary()
                (True, False)

            Their cup product is a two-cocycle that is again non-trivial in
            cohomology::

                sage: l12 = l1.cup_product(l2)
                sage: l12
                \chi_(1, 3, 6) - \chi_(2, 4, 5) - \chi_(4, 5, 6)
                sage: l1.parent().degree(), l2.parent().degree(), l12.parent().degree()
                (1, 1, 2)
                sage: l12.is_cocycle(), l12.is_coboundary()
                (True, False)
            """
            if not isinstance(cochain.parent(), Cochains):
                raise ValueError('argument must be a cochain')
            if cochain.parent().cell_complex() != self.parent().cell_complex():
                raise ValueError('cochain must be over the same cell complex')
            left_deg = self.parent().degree()
            right_deg = cochain.parent().degree()
            left_chains = self.parent().dual()
            right_chains = cochain.parent().dual()
            base_ring = coercion_model.common_parent(
                left_chains.base_ring(), right_chains.base_ring())
            cx = self.parent().cell_complex()
            codomain = cx.n_chains(
                left_deg + right_deg, base_ring=base_ring, cochains=True)
            accumulator = codomain.zero()
            for cell in codomain.indices():
                for (coeff, left_cell, right_cell) in cx.alexander_whitney(cell, left_deg):
                        if not coeff:
                            continue
                        left = left_chains(left_cell)
                        right = right_chains(right_cell)
                        accumulator += codomain(cell) * coeff * self.eval(left) * cochain.eval(right)
            return accumulator
