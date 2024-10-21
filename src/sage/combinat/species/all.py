# sage_setup: distribution = sagemath-combinat
r"""
Combinatorial species

.. TODO:: Short blurb about species

.. TODO:: Proofread / point to the main classes rather than the modules?

Introductory material
---------------------

- :ref:`section-examples-catalan`
- :ref:`section-generic-species`

Basic Species
-------------

- :ref:`sage.combinat.species.species`
- :ref:`sage.combinat.species.empty_species`
- :ref:`sage.combinat.species.recursive_species`
- :ref:`sage.combinat.species.characteristic_species`
- :ref:`sage.combinat.species.cycle_species`
- :ref:`sage.combinat.species.partition_species`
- :ref:`sage.combinat.species.permutation_species`
- :ref:`sage.combinat.species.linear_order_species`
- :ref:`sage.combinat.species.set_species`
- :ref:`sage.combinat.species.subset_species`
- :ref:`sage.combinat.species.library`

Operations on Species
---------------------

- :ref:`sage.combinat.species.sum_species`
- :ref:`sage.combinat.species.product_species`
- :ref:`sage.combinat.species.composition_species`
- :ref:`sage.combinat.species.functorial_composition_species`

Miscellaneous
-------------

- :ref:`sage.combinat.species.structure`
- :ref:`sage.combinat.species.misc`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import
lazy_import("sage.combinat.species.recursive_species", "CombinatorialSpecies")
lazy_import("sage.combinat.species", "library", as_='species')
del lazy_import
del install_doc
