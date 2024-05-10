r"""
Posets

Common posets can be accessed through ``posets.<tab>`` and are listed in the
posets catalog:

- :ref:`Catalog of posets and lattices <sage.combinat.posets.poset_examples>`

Poset-related classes:

- :ref:`sage.combinat.posets.posets`
- :ref:`sage.combinat.posets.lattices`

- :ref:`sage.combinat.posets.linear_extensions`
- :ref:`sage.combinat.posets.d_complete`
- :ref:`sage.combinat.posets.forest`
- :ref:`sage.combinat.posets.mobile`
- :ref:`sage.combinat.posets.incidence_algebras`

- :ref:`sage.combinat.posets.cartesian_product`

- :ref:`sage.combinat.posets.moebius_algebra`

- :ref:`sage.combinat.tamari_lattices`
- :ref:`sage.combinat.interval_posets`
- :ref:`sage.combinat.shard_order`

If you are looking for Poset-related :mod:`categories
<sage.categories.category>`, see
:class:`~sage.categories.posets.Posets`,
:class:`~sage.categories.finite_posets.FinitePosets`,
:class:`~sage.categories.lattice_posets.LatticePosets` and
:class:`~sage.categories.finite_lattice_posets.FiniteLatticePosets`.
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from .posets import Poset

from .lattices import LatticePoset, MeetSemilattice, JoinSemilattice

from .poset_examples import posets, Posets
