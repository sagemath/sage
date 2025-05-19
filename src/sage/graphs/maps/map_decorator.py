"""Define the internal CheckValid function."""

from typing import Type


def CheckValid(cls: Type) -> Type:
    """
    INPUT:

        -``cls`` ; a class

    OUTPUT:

        The same class but with a check before each method

    EXAMPLES::

        sage: from sage.graphs.maps.map_decorator import CheckValid
        sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
        sage: CheckValid(TopologicalDemiEdge)
        <class 'sage.graphs.maps.topological_demi_edge.TopologicalDemiEdge'>

    NOTE:

        Used internally
    """
    originalMethods = {name: method for name,
                       method in cls.__dict__.items() if callable(method)}

    def wrapper(name):
        def functionWithValidChecking(self, *args, **kwargs):

            if name != "__init__" and not self._isValid:
                raise ValueError(
                    "self isn't anymore a valid TopologicalDemiEdge")
            return originalMethods[name](self, *args, **kwargs)
        return functionWithValidChecking

    for name in originalMethods.keys():
        setattr(cls, name, wrapper(name))

    return cls
