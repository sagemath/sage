r"""
Crystals

- :ref:`sage.combinat.crystals.crystals`
- The `Lie Methods and Related Combinatorics <../../../../../thematic_tutorials/lie.html>`_ thematic tutorial
- :ref:`sage.combinat.crystals.catalog`

See also the categories for crystals: :class:`Crystals`,
:class:`HighestWeightCrystals`, :class:`FiniteCrystals`,
:class:`ClassicalCrystals`, :class:`RegularCrystals`,
:class:`~sage.categories.regular_supercrystals.RegularSuperCrystals`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.crystals', 'catalog', 'crystals')
del lazy_import
del install_doc
