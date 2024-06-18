r"""
Symplectic vector spaces

AUTHORS:

- Tobias Diez (2021): initial version
"""

# *****************************************************************************
#       Copyright (C) 2020 Tobias Diez
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import annotations

from typing import Optional, Tuple

from sage.categories.manifolds import Manifolds
from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
from sage.manifolds.differentiable.symplectic_form import (SymplecticForm,
                                                           SymplecticFormParal)
from sage.rings.real_mpfr import RR


class StandardSymplecticSpace(EuclideanSpace):
    r"""
    The vector space `\RR^{2n}` equipped with its standard symplectic form.
    """

    _symplectic_form: SymplecticForm

    def __init__(
        self,
        dimension: int,
        name: Optional[str] = None,
        latex_name: Optional[str] = None,
        coordinates: str = "Cartesian",
        symbols: Optional[str] = None,
        symplectic_name: Optional[str] = "omega",
        symplectic_latex_name: Optional[str] = None,
        start_index: int = 1,
        base_manifold: Optional[StandardSymplecticSpace] = None,
        names: Optional[Tuple[str]] = None,
    ):
        r"""
        INPUT:

        - ``dimension`` -- dimension of the space over the real field (has to be even)
        - ``name`` -- name (symbol) given to the underlying vector space;
          if ``None``, the name will be set to ``'Rn'``, where ``n`` is the ``dimension``
        - ``latex_name`` -- LaTeX symbol to denote the underlying vector space;
          if ``None``, it is set to ``name``
        - ``coordinates`` -- (default: ``'Cartesian'``) the
          type of coordinates to be initialized at the Euclidean space
          creation; allowed values are

            - ``'Cartesian'`` (canonical coordinates on `\RR^{2n}`)
            - ``'polar'`` for ``dimension=2`` only (see
                :meth:`~sage.manifolds.differentiable.examples.euclidean.EuclideanPlane.polar_coordinates`)

        - ``symbols`` -- the coordinate text symbols and LaTeX symbols, with the same conventions as the
            argument ``coordinates`` in :class:`~sage.manifolds.differentiable.chart.RealDiffChart`, namely
            ``symbols`` is a string of coordinate fields separated by a blank
            space, where each field contains the coordinate's text symbol and
            possibly the coordinate's LaTeX symbol (when the latter is different
            from the text symbol), both symbols being separated by a colon
            (``:``); if ``None``, the symbols will be automatically generated
            according to the value of ``coordinates``
        - ``symplectic_name`` -- name (symbol) given to the symplectic form
        - ``symplectic_latex_name`` -- LaTeX symbol to denote the symplectic form;
            if none is provided, it is set to ``symplectic_name``
        - ``start_index`` -- lower value of the range of
            indices used for "indexed objects" in the vector space, e.g.
            coordinates of a chart
        - ``base_manifold`` -- if not ``None``, the created object is then an open subset
            of ``base_manifold``
        - ``names`` -- (default: ``None``) unused argument, except if
            ``symbols`` is not provided; it must then be a tuple containing
            the coordinate symbols (this is guaranteed if the shortcut operator
            ``<,>`` is used)
            If ``names`` is specified, then ``dimension`` does not have to be specified.

        EXAMPLES:

        Standard symplectic form on `\RR^2`::

            sage: M.<q, p> = manifolds.StandardSymplecticSpace(2, symplectic_name='omega')
            sage: omega = M.symplectic_form()
            sage: omega.display()
            omega = -dq∧dp

        An isomomorphism of its tangent space (at any point) with an indefinite inner product space
        with distinguished basis::

            sage: Q_M_qp = omega[:]; Q_M_qp
            [ 0 -1]
            [ 1  0]
            sage: W_M_qp = VectorSpace(RR, 2, inner_product_matrix=Q_M_qp); W_M_qp
            Ambient quadratic space of dimension 2 over Real Field with 53 bits of precision
            Inner product matrix:
            [0.000000000000000 -1.00000000000000]
            [ 1.00000000000000 0.000000000000000]
            sage: T = M.tangent_space(M.point(), base_ring=RR); T
            Tangent space at Point on the Standard symplectic space R2
            sage: phi_M_qp = T.isomorphism_with_fixed_basis(T.default_basis(), codomain=W_M_qp); phi_M_qp
            Generic morphism:
            From: Tangent space at Point on the Standard symplectic space R2
            To:   Ambient quadratic space of dimension 2 over Real Field with 53 bits of precision
            Inner product matrix:
            [0.000000000000000 -1.00000000000000]
            [ 1.00000000000000 0.000000000000000]
        """
        # Check that manifold is even dimensional
        if dimension % 2 == 1:
            raise ValueError(
                f"the dimension of the manifold must be even but it is {dimension}"
            )
        dim_half = dimension // 2

        if names is not None and symbols is None:
            symbols = " ".join(names)

        if symbols is None:
            if dim_half == 1:
                symbols = r"q:q p:p"
            else:
                symbols_list = [
                    f"q{i}:q^{i} p{i}:p_{i}" for i in range(1, dim_half + 1)
                ]
                symbols = " ".join(symbols_list)

        if name is None:
            name = f"R{dimension}"

        category = Manifolds(RR).Smooth()

        EuclideanSpace.__init__(
            self,
            dimension,
            name,
            latex_name=latex_name,
            coordinates=coordinates,
            symbols=symbols,
            start_index=start_index,
            base_manifold=base_manifold,
            category=category,
            init_coord_methods=None,
        )

        self._symplectic_form = SymplecticFormParal(
            self, symplectic_name, symplectic_latex_name
        )
        for i in range(0, dim_half):
            q_index = 2 * i + 1
            self._symplectic_form.set_comp()[q_index, q_index + 1] = -1

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: V.<q, p> = manifolds.StandardSymplecticSpace(2, symplectic_name='omega'); V
            Standard symplectic space R2
        """
        return f"Standard symplectic space {self._name}"

    def symplectic_form(self) -> SymplecticForm:
        r"""
        Return the symplectic form.

        EXAMPLES:

        Standard symplectic form on `\RR^2`::

            sage: M.<q, p> = manifolds.StandardSymplecticSpace(2, symplectic_name='omega')
            sage: omega = M.symplectic_form()
            sage: omega.display()
            omega = -dq∧dp
        """
        return self._symplectic_form
