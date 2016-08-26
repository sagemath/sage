"""
The purpose of this file is to allow instancemethods to be pickable.
This is needed in setup.py only and not installed anywhere, thus is not
a library file. It solves the issue from :trac:`11874`.
"""
import types
import copy
from six.moves import copyreg
from six import get_method_function, get_method_self

# The following four methods are taken from twisted.persisted.styles
# into sage.misc.fpickle, and then here. It was needed due to
# chicken vs egg issue, as we cannot use sage.misc.fpickle in
# setup.py of sage-***.spkg


def pickleMethod(method):
    'support function for copyreg to pickle method refs'
    return unpickleMethod, (get_method_function(method).__name__,
                            get_method_self(method),
                            get_method_self(method).__class__)


def unpickleMethod(im_name,
                   im_self,
                   im_class):
    'support function for copyreg to unpickle method refs'
    try:
        unbound = getattr(im_class, im_name)
        if im_self is None:
            return unbound
        bound = types.MethodType(get_method_function(unbound),
                                 im_self)
        return bound
    except AttributeError:
        # assert im_self is not None,"No recourse: no instance to guess from."
        # Attempt a common fix before bailing -- if classes have
        # changed around since we pickled this method, we may still be
        # able to get it by looking on the instance's current class.
        unbound = getattr(im_self.__class__, im_name)
        if im_self is None:
            return unbound
        bound = types.MethodType(get_method_function(unbound),
                                 im_self)
        return bound

copyreg.pickle(types.MethodType,
                pickleMethod,
                unpickleMethod)

oldModules = {}


def pickleModule(module):
    'support function for copyreg to pickle module refs'
    return unpickleModule, (module.__name__,)


def unpickleModule(name):
    'support function for copyreg to unpickle module refs'
    if name in oldModules:
        name = oldModules[name]
    return __import__(name, {}, {}, 'x')


copyreg.pickle(types.ModuleType,
               pickleModule,
               unpickleModule)
