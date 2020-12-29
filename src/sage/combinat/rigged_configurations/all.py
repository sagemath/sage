r"""
Features that are imported by default in the interpreter namespace
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.rigged_configurations.rigged_configurations',
            'RiggedConfigurations')
