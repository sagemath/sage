"""
Base Class to Support Method Decorators

AUTHOR:

- Martin Albrecht (2009-05): inspired by a conversation with and code by Mike Hansen
"""
from sage.structure.sage_object import SageObject


class MethodDecorator(SageObject):
    def __init__(self, f):
        """
        EXAMPLES::

            sage: from sage.misc.method_decorator import MethodDecorator
            sage: class Foo:
            ....:     @MethodDecorator
            ....:     def bar(self, x):
            ....:         return x**2
            sage: J = Foo()
            sage: J.bar
            <sage.misc.method_decorator.MethodDecorator object at ...>
        """
        self.f = f
        if hasattr(f, "__doc__"):
            self.__doc__ = f.__doc__
        else:
            self.__doc__ = f.__doc__
        if hasattr(f, "__name__"):
            self.__name__ = f.__name__
        self.__module__ = f.__module__

    def _sage_src_(self):
        """
        Return the source code for the wrapped function.

        EXAMPLES:

        This class is rather abstract so we showcase its features
        using one of its subclasses::

            sage: P.<x,y,z> = PolynomialRing(ZZ)
            sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
            sage: "primary" in I.primary_decomposition._sage_src_() # indirect doctest
            True
        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.f)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES:

        This class is rather abstract so we showcase its features
        using one of its subclasses::

            sage: P.<x,y,z> = PolynomialRing(Zmod(126))
            sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
            sage: I.primary_decomposition() # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring must be a field for function 'primary_decomposition'.
        """
        return self.f(self._instance, *args, **kwds)

    def __get__(self, inst, cls=None):
        """
        EXAMPLES:

        This class is rather abstract so we showcase its features
        using one of its subclasses::

            sage: P.<x,y,z> = PolynomialRing(Zmod(126))
            sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
            sage: I.primary_decomposition() # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring must be a field for function 'primary_decomposition'.
        """
        self._instance = inst
        return self
