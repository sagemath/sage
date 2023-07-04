r"""
Time Series (deprecated module)

This module consists only of deprecated lazy imports from
:mod:`sage.stats.time_series`.
"""


from sage.misc.lazy_import import lazy_import
lazy_import('sage.stats.time_series', ('TimeSeries', 'autoregressive_fit'), deprecation=32427)
