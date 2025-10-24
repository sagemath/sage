r"""
Algebraic properties of hypergeometric functions.

AUTHORS:

- Xavier Caruso, Florian Fürnsinn (2025-10): initial version
"""

from sage.structure.sage_object import SageObject

class HypergeometricAlgebraic(SageObject):
    def __init__(self, a, b, x):
        self._a = tuple(a)
        self._b = tuple(b)
        self._x = x

    def _repr_(self):
        return "hypergeometric(%s, %s, %s)" % (self._a, self._b, self._x)
