r"""
Advanced matroid functionality

This module collects a number of advanced functions which are not directly
available to the end user by default. To import them into the main namespace,
type::

    sage: from sage.matroids.advanced import *

This adds the following to the main namespace:

    - Matroid classes:

        - :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`
        - :class:`CircuitsMatroid <sage.matroids.circuits_matroid.CircuitsMatroid>`
        - :class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`
        - :class:`DualMatroid <sage.matroids.dual_matroid.DualMatroid>`
        - :class:`FlatsMatroid <sage.matroids.flats_matroid.FlatsMatroid>`
        - :class:`GraphicMatroid <sage.matroids.graphic_matroid.GraphicMatroid>`
        - :class:`LinearMatroid <sage.matroids.linear_matroid.LinearMatroid>`
        - :class:`RegularMatroid <sage.matroids.linear_matroid.RegularMatroid>`
        - :class:`BinaryMatroid <sage.matroids.linear_matroid.BinaryMatroid>`
        - :class:`TernaryMatroid <sage.matroids.linear_matroid.TernaryMatroid>`
        - :class:`QuaternaryMatroid <sage.matroids.linear_matroid.QuaternaryMatroid>`
        - :class:`MinorMatroid <sage.matroids.minor_matroid.MinorMatroid>`
        - :class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>`

    Note that you can construct all of these through the
    :class:`Matroid() <sage.matroids.matroid.Matroid>` class, which is
    available on startup. Using the classes directly can sometimes be useful
    for faster code (e.g. if your code calls ``Matroid()`` frequently).

    - Other classes:

        - :class:`LinearSubclasses <sage.matroids.extension.LinearSubclasses>`
        - :class:`MatroidExtensions <sage.matroids.extension.MatroidExtensions>`

    Instances of these classes are returned by the methods
    :meth:`Matroid.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`
    and
    :meth:`Matroid.extensions() <sage.matroids.matroid.Matroid.extensions>`.

    - Useful functions:

        - :func:`setprint() <sage.matroids.utilities.setprint>`
        - :func:`newlabel() <sage.matroids.utilities.newlabel>`
        - :func:`get_nonisomorphic_matroids() <sage.matroids.utilities.get_nonisomorphic_matroids>`
        - :func:`lift_cross_ratios() <sage.matroids.utilities.lift_cross_ratios>`
        - :func:`lift_map() <sage.matroids.utilities.lift_map>`
        - :func:`cmp_elements_key() <sage.matroids.utilities.cmp_elements_key>`

AUTHORS:

- Stefan van Zwam (2013-04-01): initial version
"""

from sage.matroids import matroid, basis_exchange_matroid, lean_matrix

from sage.matroids.basis_matroid import BasisMatroid
from sage.matroids.circuits_matroid import CircuitsMatroid
from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
from sage.matroids.dual_matroid import DualMatroid
from sage.matroids.flats_matroid import FlatsMatroid
from sage.matroids.linear_matroid import LinearMatroid, RegularMatroid, BinaryMatroid, TernaryMatroid, QuaternaryMatroid
from sage.matroids.minor_matroid import MinorMatroid
from sage.matroids.rank_matroid import RankMatroid
from sage.matroids.union_matroid import MatroidUnion, MatroidSum, PartitionMatroid

# lazy import of GraphicMatroid for modularization purposes
from sage.misc.lazy_import import lazy_import
lazy_import('sage.matroids.graphic_matroid', 'GraphicMatroid')

from sage.matroids.extension import LinearSubclasses, MatroidExtensions
from sage.matroids.utilities import setprint, newlabel, get_nonisomorphic_matroids, lift_cross_ratios, lift_map, cmp_elements_key
