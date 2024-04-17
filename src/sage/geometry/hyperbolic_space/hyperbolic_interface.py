# sage_setup: distribution = sagemath-symbolics
r"""
Interface to Hyperbolic Models

This module provides a convenient interface for interacting with models
of hyperbolic space as well as their points, geodesics, and isometries.

The primary point of this module is to allow the code that implements
hyperbolic space to be sufficiently decoupled while still providing a
convenient user experience.

The interfaces are by default given abbreviated names.  For example,
UHP (upper half plane model), PD (Poincaré disk model), KM (Klein disk
model), and HM (hyperboloid model).

.. NOTE::

    All of the current models of 2 dimensional hyperbolic space
    use the upper half plane model for their computations.  This can
    lead to some problems, such as long coordinate strings for symbolic
    points.  For example, the vector ``(1, 0, sqrt(2))`` defines a point
    in the hyperboloid model.  Performing mapping this point to the upper
    half plane and performing computations there may return with vector
    whose components are unsimplified strings have several ``sqrt(2)``'s.
    Presently, this drawback is outweighted by the rapidity with which new
    models can be implemented.

AUTHORS:

- Greg Laun (2013): Initial version.
- Rania Amer, Jean-Philippe Burelle, Bill Goldman, Zach Groton,
  Jeremy Lent, Leila Vaden, Derrick Wigglesworth (2011): many of the
  methods spread across the files.

EXAMPLES::

    sage: HyperbolicPlane().UHP().get_point(2 + I)
    Point in UHP I + 2

    sage: HyperbolicPlane().PD().get_point(1/2 + I/2)
    Point in PD 1/2*I + 1/2
"""
# **********************************************************************
#
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# **********************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.geometry.hyperbolic_space.hyperbolic_model import (
        HyperbolicModelUHP, HyperbolicModelPD,
        HyperbolicModelHM, HyperbolicModelKM)


def HyperbolicSpace(n):
    """
    Return ``n`` dimensional hyperbolic space.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicSpace
        sage: HyperbolicSpace(2)
        Hyperbolic plane
    """
    if n == 2:
        return HyperbolicPlane()
    raise NotImplementedError("currently only implemented in dimension 2")


class HyperbolicPlane(Parent, UniqueRepresentation):
    r"""
    The hyperbolic plane `\mathbb{H}^2`.

    Here are the models currently implemented:

    - ``UHP`` -- upper half plane
    - ``PD`` -- Poincaré disk
    - ``KM`` -- Klein disk
    - ``HM`` -- hyperboloid model
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: TestSuite(H).run()
        """
        Parent.__init__(self, category=Sets().Metric().WithRealizations())
        self.a_realization() # We create a realization so at least one is known

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: HyperbolicPlane()
            Hyperbolic plane
        """
        return "Hyperbolic plane"

    def a_realization(self):
        """
        Return a realization of ``self``.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: H.a_realization()
            Hyperbolic plane in the Upper Half Plane Model
        """
        return self.UHP()

    UHP = HyperbolicModelUHP
    UpperHalfPlane = UHP

    PD = HyperbolicModelPD
    PoincareDisk = PD

    KM = HyperbolicModelKM
    KleinDisk = KM

    HM = HyperbolicModelHM
    Hyperboloid = HM


class HyperbolicModels(Category_realization_of_parent):
    r"""
    The category of hyperbolic models of hyperbolic space.
    """
    def __init__(self, base):
        r"""
        Initialize the hyperbolic models of hyperbolic space.

        INPUT:

        - ``base`` -- a hyperbolic space

        TESTS::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
            sage: H = HyperbolicPlane()
            sage: models = HyperbolicModels(H)
            sage: H.UHP() in models
            True
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
            sage: H = HyperbolicPlane()
            sage: HyperbolicModels(H)
            Category of hyperbolic models of Hyperbolic plane
        """
        return "Category of hyperbolic models of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
            sage: H = HyperbolicPlane()
            sage: models = HyperbolicModels(H)
            sage: models.super_categories()
            [Category of metric spaces,
             Category of realizations of Hyperbolic plane]
        """
        return [Sets().Metric(), Realizations(self.base())]

    class ParentMethods:
        def _an_element_(self):
            """
            Return an element of ``self``.

            EXAMPLES::

                sage: H = HyperbolicPlane()
                sage: H.UHP().an_element()
                Point in UHP I
                sage: H.PD().an_element()
                Point in PD 0
                sage: H.KM().an_element()
                Point in KM (0, 0)
                sage: H.HM().an_element()
                Point in HM (0, 0, 1)
            """
            from sage.rings.integer_ring import ZZ
            return self(self.realization_of().PD().get_point(ZZ.zero()))
