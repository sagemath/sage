r"""
Catalog of designs

This module gathers all designs that can be reached through
``designs.<tab>``. Example with the the Witt design on 24 points::

    sage: designs.WittDesign(24) # optional - gap_packages
    Incidence structure with 24 points and 759 blocks

Or a Steiner Triple System on 19 points::

    sage: designs.steiner_triple_system(19)
    (19,3,1)-Balanced Incomplete Block Design

**La Jolla Covering Repository**

The La Jolla Covering Repository (LJCR, see [1]_) is an online database of
covering designs. As it is frequently updated, it is not included in Sage, but
one can query it through :meth:`designs.best_known_covering_design_from_LJCR
<sage.combinat.designs.covering_design.best_known_covering_design_www>`::

    sage: C = designs.best_known_covering_design_from_LJCR(7, 3, 2)   # optional - internet
    sage: C                            # optional - internet
    (7,3,2)-covering design of size 7
    Lower bound: 7
    Method: lex covering
    Submitted on: 1996-12-01 00:00:00
    sage: C.incidence_structure()      # optional - internet
    Incidence structure with 7 points and 7 blocks

**Design constructors**

This module gathers the following designs :

.. csv-table::
    :class: contentstable
    :widths: 100
    :delim: |

    :meth:`~sage.combinat.designs.block_design.ProjectiveGeometryDesign`
    :meth:`~sage.combinat.designs.block_design.DesarguesianProjectivePlaneDesign`
    :meth:`~sage.combinat.designs.block_design.HughesPlane`
    :meth:`~sage.combinat.designs.database.HigmanSimsDesign`
    :meth:`~sage.combinat.designs.bibd.balanced_incomplete_block_design`
    :meth:`~sage.combinat.designs.resolvable_bibd.resolvable_balanced_incomplete_block_design`
    :meth:`~sage.combinat.designs.resolvable_bibd.kirkman_triple_system`
    :meth:`~sage.combinat.designs.block_design.AffineGeometryDesign`
    :meth:`~sage.combinat.designs.block_design.CremonaRichmondConfiguration`
    :meth:`~sage.combinat.designs.block_design.WittDesign`
    :meth:`~sage.combinat.designs.block_design.HadamardDesign`
    :meth:`~sage.combinat.designs.block_design.Hadamard3Design`
    :meth:`~sage.combinat.designs.latin_squares.mutually_orthogonal_latin_squares`
    :meth:`~sage.combinat.designs.orthogonal_arrays.transversal_design`
    :meth:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
    :meth:`~sage.combinat.designs.orthogonal_arrays.incomplete_orthogonal_array`
    :meth:`~sage.combinat.designs.difference_family.difference_family`
    :meth:`~sage.combinat.designs.difference_matrices.difference_matrix`
    :meth:`~sage.combinat.designs.bibd.steiner_triple_system`
    :meth:`~sage.combinat.designs.steiner_quadruple_systems.steiner_quadruple_system`
    :meth:`~sage.combinat.designs.block_design.projective_plane`

And the :meth:`designs.best_known_covering_design_from_LJCR
<sage.combinat.designs.covering_design.best_known_covering_design_www>` function
which queries the LJCR.

.. TODO::

    Implement DerivedDesign and ComplementaryDesign.

REFERENCES:

.. [1] La Jolla Covering Repository,
  http://www.ccrwest.org/cover.html
"""
from __future__ import absolute_import
from sage.combinat.designs.block_design import (BlockDesign,
                                                ProjectiveGeometryDesign,
                                                DesarguesianProjectivePlaneDesign,
                                                projective_plane,
                                                AffineGeometryDesign,
                                                WittDesign,
                                                HadamardDesign,
                                                Hadamard3Design,
                                                HughesPlane,
                                                CremonaRichmondConfiguration)

from .database import HigmanSimsDesign

from sage.combinat.designs.steiner_quadruple_systems import steiner_quadruple_system

from sage.combinat.designs.covering_design import best_known_covering_design_www as best_known_covering_design_from_LJCR

from sage.combinat.designs.latin_squares import mutually_orthogonal_latin_squares

from sage.combinat.designs.orthogonal_arrays import transversal_design, incomplete_orthogonal_array


from sage.combinat.designs.difference_family import difference_family
from .difference_matrices import difference_matrix

from sage.misc.superseded import deprecated_callable_import
deprecated_callable_import(19096,
                           'sage.combinat.designs.incidence_structures',
                           globals(),
                           locals(),
                           ["IncidenceStructure"],
                           ("This alias will soon be removed. You can call the same object by removing 'designs.' in your command"))

Hypergraph = BlockDesign = IncidenceStructure    # just an alias
from sage.combinat.designs.bibd import balanced_incomplete_block_design, steiner_triple_system
from sage.combinat.designs.resolvable_bibd import resolvable_balanced_incomplete_block_design, kirkman_triple_system
from sage.combinat.designs.group_divisible_designs import group_divisible_design

from .orthogonal_arrays import OAMainFunctions as orthogonal_arrays

# When this deprecated function is removed, remove the handling of k=None in the
# function orthogonal_arrays.orthogonal_array()
deprecated_callable_import(17034,
                           'sage.combinat.designs.orthogonal_arrays',
                           globals(),
                           locals(),
                           ["orthogonal_array"],
                           ("This function will soon be removed. Use the designs.orthogonal_arrays.* functions instead"))

# We don't want this to appear in designs.<tab>
del deprecated_callable_import
