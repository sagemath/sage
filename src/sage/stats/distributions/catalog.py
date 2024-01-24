r"""
Index of distributions

This catalogue includes the samplers for statistical distributions listed below.

Let ``<tab>`` indicate pressing the :kbd:`Tab` key.  So begin by typing
``algebras.<tab>`` to the see the currently implemented named algebras.

- :class:`distributions.discrete_gaussian_integer.DiscreteGaussianDistributionIntegerSampler
  <sage.stats.distributions.discrete_gaussian_integer.DiscreteGaussianDistributionIntegerSampler>`
- :class:`distributions.discrete_gaussian_lattice.DiscreteGaussianDistributionLatticeSampler
  <sage.stats.distributions.discrete_gaussian_lattice.DiscreteGaussianDistributionLatticeSampler>`
- :class:`distributions.discrete_gaussian_polynomial.DiscreteGaussianDistributionPolynomialSampler
  <sage.stats.distributions.discrete_gaussian_polynomial.DiscreteGaussianDistributionPolynomialSampler>`

To import these names into the global namespace, use::

    sage: from sage.stats.channels_catalog import *

"""
#*****************************************************************************
#       Copyright (C) 2024 Gareth Ma <grhkm21@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import("sage.stats.distributions.discrete_gaussian_integer", ["DiscreteGaussianDistributionIntegerSampler"])
_lazy_import("sage.stats.distributions.discrete_gaussian_lattice", ["DiscreteGaussianDistributionLatticeSampler"])
_lazy_import("sage.stats.distributions.discrete_gaussian_polynomial", ["DiscreteGaussianDistributionPolynomialSampler"])
