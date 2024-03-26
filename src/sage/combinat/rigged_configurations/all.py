r"""
Rigged configurations

.. TODO:: Proofread / point to the main classes rather than the modules?

- :ref:`sage.combinat.rigged_configurations.rc_crystal`
- :ref:`sage.combinat.rigged_configurations.rc_infinity`

- :ref:`sage.combinat.rigged_configurations.rigged_configurations`
- :ref:`sage.combinat.rigged_configurations.rigged_configuration_element`

- :ref:`sage.combinat.rigged_configurations.tensor_product_kr_tableaux`
- :ref:`sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element`
- :ref:`sage.combinat.rigged_configurations.kr_tableaux`

- :ref:`sage.combinat.rigged_configurations.kleber_tree`

- :ref:`sage.combinat.rigged_configurations.rigged_partition`

Bijections
----------

- :ref:`sage.combinat.rigged_configurations.bijection`
- :ref:`sage.combinat.rigged_configurations.bij_abstract_class`
- :ref:`sage.combinat.rigged_configurations.bij_type_A`
- :ref:`sage.combinat.rigged_configurations.bij_type_B`
- :ref:`sage.combinat.rigged_configurations.bij_type_C`
- :ref:`sage.combinat.rigged_configurations.bij_type_D`
- :ref:`sage.combinat.rigged_configurations.bij_type_A2_odd`
- :ref:`sage.combinat.rigged_configurations.bij_type_A2_even`
- :ref:`sage.combinat.rigged_configurations.bij_type_A2_dual`
- :ref:`sage.combinat.rigged_configurations.bij_type_D_twisted`
- :ref:`sage.combinat.rigged_configurations.bij_type_D_tri`
- :ref:`sage.combinat.rigged_configurations.bij_infinity`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.rigged_configurations.rigged_configurations',
            'RiggedConfigurations')
del lazy_import
del install_doc
