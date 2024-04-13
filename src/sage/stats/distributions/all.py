# sage_setup: distribution = sagemath-modules
# We lazy_import the following modules since they import numpy which
# slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.stats.distributions.discrete_gaussian_integer", ["DiscreteGaussianDistributionIntegerSampler"])
lazy_import("sage.stats.distributions.discrete_gaussian_lattice", ["DiscreteGaussianDistributionLatticeSampler"])
lazy_import("sage.stats.distributions.discrete_gaussian_polynomial", ["DiscreteGaussianDistributionPolynomialSampler"])
