r"""
Crystals

Introductory material
---------------------

- :ref:`sage.combinat.crystals.crystals`
- The `Lie Methods and Related Combinatorics <../../../../../thematic_tutorials/lie.html>`_ thematic tutorial

Catalogs of crystals
--------------------

- :ref:`sage.combinat.crystals.catalog`

See also
--------

- The categories for crystals: :class:`Crystals`, :class:`HighestWeightCrystals`,
  :class:`FiniteCrystals`, :class:`ClassicalCrystals`, :class:`RegularCrystals`,
  :class:`~sage.categories.regular_supercrystals.RegularSuperCrystals`
  -- The categories for crystals
- :ref:`sage.combinat.root_system.all`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.crystals', 'catalog', 'crystals')

