"""
Quantum Groups
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.quantum_groups.fock_space', 'FockSpace')
lazy_import('sage.algebras.quantum_groups.quantum_group_gap', 'QuantumGroup')
del lazy_import
