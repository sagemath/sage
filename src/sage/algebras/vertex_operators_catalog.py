r"""
Catalog of Vertex Operators
"""
from sage.misc.lazy_import import lazy_import

lazy_import('sage.algebras.vertex_operators', 'BosonicFockSpace')
lazy_import('sage.algebras.vertex_operators', 'CreationOperator')
lazy_import('sage.algebras.vertex_operators', 'AnnihilationOperator')

del lazy_import  # We remove the object from here so it doesn't appear under tab completion
