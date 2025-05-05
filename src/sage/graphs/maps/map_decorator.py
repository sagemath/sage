

def CheckValid(cls):
    """
    INPUT:
        cls a class

    OUTPUT:
        the same class but with a check before each method 

    EXAMPLES::

        sage: CheckValid(TopologicalDemiEdge)
        <class 'TopologicalDemiEdge.TopologicalDemiEdge'>

    .. NOTE::
        Used internaly
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
