# sage_setup: distribution = sagemath-modules

from sage.misc.lazy_import import lazy_import

lazy_import("sage.stats.r", "ttest")

# We lazy_import the following modules since they import numpy which
# slows down sage startup

lazy_import("sage.stats.time_series", ["TimeSeries", "autoregressive_fit"])
lazy_import("sage.stats.intlist", ["IntList"])
del lazy_import
