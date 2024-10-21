# sage_setup: distribution = sagemath-symbolics
r"""
Manifold Structures

These classes encode the structure of a manifold.

AUTHORS:

- Travis Scrimshaw (2015-11-25): Initial version
- Eric Gourgoulhon (2015): add :class:`DifferentialStructure` and
  :class:`RealDifferentialStructure`
- Eric Gourgoulhon (2018): add :class:`PseudoRiemannianStructure`,
  :class:`RiemannianStructure` and :class:`LorentzianStructure`
"""

#*****************************************************************************
#       Copyright (C) 2015, 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.fast_methods import Singleton
from sage.manifolds.chart import Chart, RealChart
from sage.manifolds.scalarfield_algebra import ScalarFieldAlgebra
from sage.manifolds.manifold_homset import TopologicalManifoldHomset
from sage.manifolds.differentiable.chart import DiffChart, RealDiffChart
from sage.manifolds.differentiable.scalarfield_algebra import \
                                                         DiffScalarFieldAlgebra
from sage.manifolds.differentiable.manifold_homset import \
                                                   DifferentiableManifoldHomset

# This is a slight abuse by making this a Singleton, but there is no
#    need to have different copies of this object.
class TopologicalStructure(Singleton):
    """
    The structure of a topological manifold over a general topological field.
    """
    chart = Chart
    name = "topological"
    scalar_field_algebra = ScalarFieldAlgebra
    homset = TopologicalManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import TopologicalStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: TopologicalStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat


class RealTopologicalStructure(Singleton):
    r"""
    The structure of a topological manifold over `\RR`.
    """
    chart = RealChart
    name = "topological"
    scalar_field_algebra = ScalarFieldAlgebra
    homset = TopologicalManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import RealTopologicalStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: RealTopologicalStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat

class DifferentialStructure(Singleton):
    """
    The structure of a differentiable manifold over a general topological
    field.
    """
    chart = DiffChart
    name = "differentiable"
    scalar_field_algebra = DiffScalarFieldAlgebra
    homset = DifferentiableManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import DifferentialStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: DifferentialStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat


class RealDifferentialStructure(Singleton):
    r"""
    The structure of a differentiable manifold over `\RR`.
    """
    chart = RealDiffChart
    name = "differentiable"
    scalar_field_algebra = DiffScalarFieldAlgebra
    homset = DifferentiableManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import RealDifferentialStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: RealDifferentialStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat

class PseudoRiemannianStructure(Singleton):
    """
    The structure of a pseudo-Riemannian manifold.
    """
    chart = RealDiffChart
    name = "pseudo-Riemannian"
    scalar_field_algebra = DiffScalarFieldAlgebra
    homset = DifferentiableManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import PseudoRiemannianStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: PseudoRiemannianStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat

class RiemannianStructure(Singleton):
    """
    The structure of a Riemannian manifold.
    """
    chart = RealDiffChart
    name = "Riemannian"
    scalar_field_algebra = DiffScalarFieldAlgebra
    homset = DifferentiableManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import RiemannianStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: RiemannianStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat

class LorentzianStructure(Singleton):
    """
    The structure of a Lorentzian manifold.
    """
    chart = RealDiffChart
    name = "Lorentzian"
    scalar_field_algebra = DiffScalarFieldAlgebra
    homset = DifferentiableManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import LorentzianStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: LorentzianStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat

class DegenerateStructure(Singleton):
    """
    The structure of a degenerate manifold.
    """
    chart = RealDiffChart
    name = "degenerate_metric"
    scalar_field_algebra = DiffScalarFieldAlgebra
    homset = DifferentiableManifoldHomset

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import DegenerateStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: DegenerateStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision
        """
        return cat
