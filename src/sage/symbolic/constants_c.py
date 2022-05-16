r"""
The constant `e` (deprecated module)

This module consists only of deprecated lazy imports from
:mod:`sage.symbolic.expression`.
"""


from sage.misc.lazy_import import lazy_import
lazy_import('sage.symbolic.expression', 'E', deprecation=32386)
