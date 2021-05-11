r"""
Advanced matroid functionality.

This module collects a number of advanced functions which are not directly
available to the end user by default. To import them into the main namespace,
type::

    sage: from sage.matroids.advanced import *

This adds the following to the main namespace:

    - Matroid classes:
        - :class:`MinorMatroid <sage.matroids.minor_matroid.MinorMatroid>`
        - :class:`DualMatroid <sage.matroids.dual_matroid.DualMatroid>`
        - :class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>`
        - :class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`
        - :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`
        - :class:`LinearMatroid <sage.matroids.linear_matroid.LinearMatroid>`
        - :class:`RegularMatroid <sage.matroids.linear_matroid.RegularMatroid>`
        - :class:`BinaryMatroid <sage.matroids.linear_matroid.BinaryMatroid>`
        - :class:`TernaryMatroid <sage.matroids.linear_matroid.TernaryMatroid>`
        - :class:`QuaternaryMatroid <sage.matroids.linear_matroid.QuaternaryMatroid>`
        - :class:`GraphicMatroid <sage.matroids.linear_matroid.GraphicMatroid>`

    Note that you can construct all of these through the
    :func:`Matroid() <sage.matroids.constructor.Matroid>` function, which is
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
        - :func:`lift_cross_ratios() <sage.matroids.linear_matroid.lift_cross_ratios>`
        - :func:`lift_map() <sage.matroids.linear_matroid.lift_map>`

AUTHORS:

- Stefan van Zwam (2013-04-01): initial version
"""
import sage.matroids.matroid
import sage.matroids.basis_exchange_matroid
from .minor_matroid import MinorMatroid
from .dual_matroid import DualMatroid
from .rank_matroid import RankMatroid
from .circuit_closures_matroid import CircuitClosuresMatroid
from .basis_matroid import BasisMatroid
from .linear_matroid import LinearMatroid, RegularMatroid, BinaryMatroid, TernaryMatroid, QuaternaryMatroid
from .utilities import setprint, newlabel, get_nonisomorphic_matroids, lift_cross_ratios, lift_map
from . import lean_matrix
from .extension import LinearSubclasses, MatroidExtensions
from .union_matroid import MatroidUnion, MatroidSum, PartitionMatroid
from .graphic_matroid import GraphicMatroid
