r"""
Combinatorial species

A combinatorial species, as introduced by Joyal, is a functor from
the category of finite sets with bijections to the category of finite
sets with bijections.  Alternatively, we can regard a combinatorial
species as a formal sum of group actions of the symmetric groups
`\mathfrak S_n`, for `n\in\NN`.  For example, the trivial action of
`\mathfrak S_n` corresponds to the species `\mathcal E_n` of sets of
cardinality `n`.

More generally, a weighted multisort species in `k` sorts is a
functor from the category of `k`-tuples of finite sets with
bijections to the category of weighted finite sets with
weight-preserving bijections.  We may think of the sorts as
variables, traditionally denoted by `X, Y, Z`.

Such a species is equivalent to a formal sum of group actions of
groups `\mathfrak S_{n_1} \times \cdots \times \mathfrak S_{n_k}` with
`(n_1,\ldots,n_k)\in\NN^k`, together with a weight on the orbits of
each group action.  Yet more generally, a virtual weighted multisort
species is a formal difference of weighted multisort species.

We regard a combinatorial species as a sequence of group actions of
the symmetric groups `\mathfrak S_n`, for `n\in\NN`.

Coefficients of lazy species are computed on demand.  They have
infinite precision, although equality can only be decided in special
cases.

There are currently two implementations of combinatorial species, one
of which is deprecated and will be removed once all of the
functionality has been ported to the new framework.

The recommended implementation is :ref:`sage.rings.lazy_species` and
used as the following examples illustrate.

EXAMPLES:

We define rooted ordered trees with leaves having weight `q`::

    sage: R.<q> = QQ[]
    sage: L.<X> = LazyCombinatorialSpecies(R)
    sage: leaf = X
    sage: node = q * X
    sage: linorder = X/(1 - X)
    sage: T = L.undefined(valuation=1)
    sage: T.define(leaf + node * linorder(T))
    sage: T.isotype_generating_series()
    X + q*X^2 + ((q^2+q)*X^3) + ((q^3+3*q^2+q)*X^4) + ((q^4+6*q^3+6*q^2+q)*X^5)
    + ((q^5+10*q^4+20*q^3+10*q^2+q)*X^6) + O(X^7)

We define rooted unordered trees with leaves of sort `Y`.  The
standard representation of a species is its molecular decomposition::

    sage: L = LazyCombinatorialSpecies(R, "X")
    sage: E = L.Sets()
    sage: Ep = E.restrict(1)
    sage: M.<X, Y> = LazyCombinatorialSpecies(R)
    sage: A = M.undefined(valuation=1)
    sage: A.define(Y + X * Ep(A))
    sage: A.truncate(5)
    Y + X*Y + (X^2*Y+X*E_2(Y)) + (X^3*Y+X^2*E_2(Y)+X^2*Y^2+X*E_3(Y))

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
lazy_import("sage.combinat.species.recursive_species", "CombinatorialSpecies",
            deprecation=(38544, "combinat.species is superseded by LazyCombinatorialSpecies"))

lazy_import("sage.combinat.species", "library", as_='species',
            deprecation=(38544, "combinat.species is superseded by LazyCombinatorialSpecies"))
del lazy_import
del install_doc
