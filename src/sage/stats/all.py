import sage.stats.distributions.catalog as distributions

from sage.stats.r import ttest
from sage.stats.basic_stats import (mean, mode, std, variance, median, moving_average)
from sage.stats.hmm import all as hmm
from sage.stats.statistical_tests import sigmas_from_uniform

# We lazy_import the following modules since they import numpy which
# slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.stats.time_series", ["TimeSeries", "autoregressive_fit"])
lazy_import("sage.stats.intlist", ["IntList"])
