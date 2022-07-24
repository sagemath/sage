r"""
Symbolic constants (deprecated module)

This module consists only of deprecated lazy imports from
:mod:`sage.symbolic.expression`.
"""


from sage.misc.lazy_import import lazy_import
lazy_import('sage.symbolic.expression', 'PynacConstant', deprecation=32386)
