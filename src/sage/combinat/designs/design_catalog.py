r"""
Catalog of designs

This module gathers all designs that can be reached through
``designs.<tab>``. Example with the Witt design on 24 points::

    sage: designs.WittDesign(24)                                # optional - gap_package_design
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
    (7, 3, 2)-covering design of size 7
    Lower bound: 7
    Method: lex covering
    Submitted on: 1996-12-01 00:00:00
    sage: C.incidence_structure()      # optional - internet
    Incidence structure with 7 points and 7 blocks

**Design constructors**

This module gathers the following designs:

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
    :meth:`~sage.combinat.designs.biplane`
    :meth:`~sage.combinat.designs.gen_quadrangles_with_spread`

And the :meth:`designs.best_known_covering_design_from_LJCR
<sage.combinat.designs.covering_design.best_known_covering_design_www>` function
which queries the LJCR.

.. TODO::

    Implement DerivedDesign and ComplementaryDesign.

REFERENCES:

.. [1] La Jolla Covering Repository,
  https://dmgordon.org/cover
"""
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.designs.block_design',
            ('BlockDesign',
             'ProjectiveGeometryDesign',
             'DesarguesianProjectivePlaneDesign',
             'projective_plane',
             'AffineGeometryDesign',
             'WittDesign',
             'HadamardDesign',
             'Hadamard3Design',
             'HughesPlane',
             'CremonaRichmondConfiguration'))

lazy_import('sage.combinat.designs.database', 'HigmanSimsDesign')

lazy_import('sage.combinat.designs.steiner_quadruple_systems',
            'steiner_quadruple_system')

lazy_import('sage.combinat.designs.covering_design',
            'best_known_covering_design_www',
            as_='best_known_covering_design_from_LJCR')

lazy_import('sage.combinat.designs.latin_squares',
            'mutually_orthogonal_latin_squares')

lazy_import('sage.combinat.designs.orthogonal_arrays',
            ('transversal_design', 'incomplete_orthogonal_array'))

lazy_import('sage.combinat.designs.difference_family', 'difference_family')
lazy_import('sage.combinat.designs.difference_matrices', 'difference_matrix')

lazy_import('sage.combinat.designs.bibd',
            ('balanced_incomplete_block_design', 'steiner_triple_system', 'biplane'))
lazy_import('sage.combinat.designs.resolvable_bibd',
            ('resolvable_balanced_incomplete_block_design',
             'kirkman_triple_system'))
lazy_import('sage.combinat.designs.group_divisible_designs',
            'group_divisible_design')

lazy_import('sage.combinat.designs.orthogonal_arrays',
            'OAMainFunctions', as_='orthogonal_arrays')

lazy_import('sage.combinat.designs.gen_quadrangles_with_spread',
            ('generalised_quadrangle_with_spread', 'generalised_quadrangle_hermitian_with_ovoid'))
