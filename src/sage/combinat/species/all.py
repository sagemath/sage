r"""
Combinatorial species

.. TODO:: Short blurb about species

.. TODO:: Proofread / point to the main classes rather than the modules?

Introductory material
---------------------

- :ref:`section-examples-catalan`
- :ref:`section-generic-species`

Lazy Power Series
-----------------

- :ref:`sage.combinat.species.stream`
- :ref:`sage.combinat.species.series_order`
- :ref:`sage.combinat.species.series`
- :ref:`sage.combinat.species.generating_series`

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

from .series import LazyPowerSeriesRing
from .recursive_species import CombinatorialSpecies
from . import library as species
