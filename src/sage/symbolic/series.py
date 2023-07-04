r"""
Symbolic Series

This module consists only of deprecated lazy imports from
:mod:`sage.symbolic.expression`.
"""


from sage.misc.lazy_import import lazy_import
lazy_import('sage.symbolic.expression', 'SymbolicSeries', deprecation=32386)
