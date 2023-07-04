
from .r import ttest
from .basic_stats import (mean, mode, std, variance, median, moving_average)
from .hmm import all as hmm

# We lazy_import the following modules since they import numpy which
# slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.stats.time_series", ["TimeSeries", "autoregressive_fit"])
lazy_import("sage.stats.intlist", ["IntList"])
