from collections.abc import Iterable
from typing import Any

from sage.features import FeatureNotPresentError

class MockConstructor(type):
    _name: str

    def __new__(cls, name: str, bases: tuple[type, ...], attrs: dict[str, Any], module: str | None = None):
        print(f"Constructor: __new__")
        attrs["_name"] = name
        attrs["_module"] = module
        return super().__new__(cls, name, bases, attrs)

        # if len(args) == 3 and isinstance(args[1], tuple):
        #     super().__new__(cls, *args)
        # return cls
    
    # def __init__(self, *args, name: str, **kwargs):
    #     self._name = name

    # def __getattr__(self, attr: str) -> Any:
    #     # For debugging:
    #     print(f"__getattr__ {attr}")
    #     return mock(f"{self._name}.{attr}")

class MockInstance:#(metaclass=Mock):
    _name: str
    _module: str
    #__pyx_vtable__ = 1

    def __getattr__(self, attr: str) -> Any:
        # For debugging:
        print(f"__getattr__ instance {attr}")
        if attr in ("_name", "_module"):
            raise NotImplementedError(attr)
        # if attr == '__pyx_vtable__':
        #     return 1

        return mock(self._module, f"{self._name}.{attr}")

# class Wrapper(type):
#     _key: str

#     def __new__(cls, *args, **kwargs):
#         # For debugging:
#         print(f"__new__ {len(args)}")
#     #    return super().__new__(cls)
#         if len(args) == 3 and isinstance(args[1], tuple):
#             super().__new__(cls, *args)
#             #superclass = args[1][-1].__class__
#             #print(f"__new__ {args[0], superclass.__display_name__}")
#             #if superclass is cls:
#             #    return super().__new__(
#             #        args[0],
#             #        superclass.__display_name__,
#             #        superclass=superclass,
#             #        attributes=args[2])
#         #cls._key = args[0]
#         return cls
#         #raise FeatureNotPresentError(
#         #    reason=f"A dependency of the package {__name__} is not available, trying to access {args[0]}"
#         #)

#     def __init__(self, *args, **kwargs):
#         print(f"__init__")
#         #self._key = name

#     def raise_error(self):
#         raise FeatureNotPresentError(
#             reason=f"A dependency of the package {__name__} is not available, trying to access {self._key}"
#         )

#     def __basicsize__(self):
#         print(f"__basicsize__")
#         self.raise_error()

    def __call__(self, *args, **kwargs):
        # For debugging:
        print(f"__call__{self._name}")
        #for arg in args:
        #   for subarg in arg:
        #       print(type(subarg))
        self.raise_error()

    # def __mro_entries__(self, bases):
    #     """
    #     In case someone tries to inherit from this class, we substitute the bases with 
    #     a Wrapper.
    #     See https://docs.python.org/3/reference/datamodel.html#object.__mro_entries__
    #     """
    #     return (self.__class__, )

#     def __getattribute__(self, attr: str) -> Any:
#         if attr in ('_key', 'raise_error'):
#             return super().__getattribute__(attr)
#         # For debugging:
#         print(f"__getattr__ {attr}")
#         return self

    def __repr__(self):
        print(f"__repr__")
        return f"<Mock name='{self._name}' module='{self._module}'>"
        #return super().__repr__()
        #self.raise_error()

#     def __eq__(self, value: object) -> bool:
#         print(f"__eq__ {value}")
#         return super().__eq__(value)
    
#     def __subclasscheck__(self, subclass: type) -> bool:
#         print(f"__subclasscheck__ {subclass}")
#         return True
    
#     def __instancecheck__(self, instance: object) -> bool:
#         print(f"__instancecheck__ {instance}")
#         return True

#     @classmethod
#     def __subclasshook__(cls, subclass: type) -> bool:
#         print(f"__subclasshook__ {subclass}")
#         return True

#     def __init_subclass__(cls, **kwargs):
#         print(f"__init_subclass__ {kwargs}")
#         return super().__init_subclass__(**kwargs)


# class Mock:#(metaclass=Wrapper):
#     _key: str

#     def __init__(self, *args, **kwargs):
#         print(f"__init__")
#         self._key = args[0]

    def raise_error(self):
        raise FeatureNotPresentError(
            reason=f"A dependency of the package '{self._module}' is not available, trying to access '{self._name}'"
        )

#     def __basicsize__(self):
#         print(f"__basicsize__")
#         self.raise_error()

#     def __call__(self, *args, **kwargs):
#         # For debugging:
#         print(f"__call__{self._key}")
#         for arg in args:
#            for subarg in arg:
#                print(type(subarg))
#         self.raise_error()

    def __mro_entries__(self, bases):
        """
        In case someone tries to inherit from this class, we substitute the bases with 
        a Wrapper.
        See https://docs.python.org/3/reference/datamodel.html#object.__mro_entries__
        """
        print(f"__mro_entries__ {bases}")
        #MockCls = MockConstructor(self._name, (MockInstance, ), {}, module=self._module)
        #return (MockCls, )
        MockCls = self.__class__
        MockCls._module = self._module
        MockCls._name = self._name
        return (MockCls, )

#     def __getattribute__(self, attr: str) -> Any:
#         if attr in ('_key', 'raise_error'):
#             return super().__getattribute__(attr)
#         # For debugging:
#         print(f"__getattr__ {attr}")
#         return self

#     def __repr__(self):
#         print(f"__repr__")
#         return self._key
#         #return super().__repr__()
#         #self.raise_error()

#     def __eq__(self, value: object) -> bool:
#         print(f"__eq__ {value}")
#         return super().__eq__(value)
    
#     def __subclasscheck__(self, subclass: type) -> bool:
#         print(f"__subclasscheck__ {subclass}")
#         return True
    
#     def __instancecheck__(self, instance: object) -> bool:
#         print(f"__instancecheck__ {instance}")
#         return True

#     @classmethod
#     def __subclasshook__(cls, subclass: type) -> bool:
#         print(f"__subclasshook__ {subclass}")
#         return True

#     def __init_subclass__(cls, **kwargs):
#         print(f"__init_subclass__ {kwargs}")
#         return super().__init_subclass__(**kwargs)


def mock(module: str, name: str) -> MockInstance:
    r"""
    Return a mock object.

    EXAMPLES::

        sage: from sage.features.mock import mock
        sage: mock('foo')
        <Mock name='foo'>
    """
    print(f"mock {module}.{name}")
    mock = MockInstance()
    mock._module = module
    mock._name = name
    return mock
    #MockCls = MockConstructor(name, (MockInstance, ), {}, module=module)
    #return MockCls()
    #return #MockInstance(module, name)
    #return Mock.__new__(MockInstance, (), {})
