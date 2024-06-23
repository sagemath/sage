# sage_setup: distribution = sagemath-combinat
r"""
Combinatorics

Introductory material
---------------------

- :ref:`sage.combinat.quickref`
- :ref:`sage.combinat.tutorial`

Thematic indexes
----------------

- :ref:`sage.combinat.algebraic_combinatorics`

  - :ref:`sage.combinat.chas.all`
  - :ref:`sage.combinat.cluster_algebra_quiver.all`
  - :ref:`sage.combinat.crystals.all`
  - :ref:`sage.combinat.root_system.all`
  - :ref:`sage.combinat.sf.all`
  - :class:`~sage.combinat.fully_commutative_elements.FullyCommutativeElements`

- :ref:`sage.combinat.counting`
- :ref:`sage.combinat.enumerated_sets`
- :ref:`sage.combinat.catalog_partitions`
- :ref:`sage.combinat.finite_state_machine`
- :ref:`sage.combinat.species.all`
- :ref:`sage.combinat.designs.all`
- :ref:`sage.combinat.posets.all`
- :ref:`sage.combinat.words.all`
- :ref:`sage.combinat.bijectionist`

Utilities
---------

- :ref:`sage.combinat.output`
- :ref:`sage.combinat.ranker`
- :func:`Combinatorial maps <sage.combinat.combinatorial_map.combinatorial_map>`
- :ref:`sage.combinat.misc`

Related topics
--------------

- :ref:`sage.coding`
- :ref:`sage.dynamics`
- :ref:`sage.graphs`
"""
from sage.misc.namespace_package import install_doc, install_dict
# install the docstring of this module to the containing package
install_doc(__package__, __doc__)

# install modules quickref and tutorial to the containing package
from sage.combinat import quickref, tutorial
install_dict(__package__, {'quickref': quickref, 'tutorial': tutorial})
del quickref, tutorial

from sage.combinat.all__sagemath_categories import *

from sage.combinat.all__sagemath_combinat import *
from sage.combinat.all__sagemath_modules import *
from sage.combinat.all__sagemath_graphs import *
del lazy_import
del install_dict
del install_doc
