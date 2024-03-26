# sage_setup: distribution = sagemath-combinat
"""
Combinatorial Hopf algebras

- :ref:`sage.combinat.sf.all`
- :ref:`sage.combinat.ncsf_qsym.all`
- :ref:`sage.combinat.ncsym.all`
- :ref:`sage.combinat.schubert_polynomial`
- :ref:`sage.combinat.chas.fsym`
- :ref:`sage.combinat.fqsym`
- :ref:`sage.combinat.grossman_larson_algebras`
- :ref:`sage.combinat.chas.wqsym`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.chas.fsym', ['FreeSymmetricFunctions'])
lazy_import('sage.combinat.chas.wqsym', ['WordQuasiSymmetricFunctions'])
del lazy_import
del install_doc
