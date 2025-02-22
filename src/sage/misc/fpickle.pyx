# cython: old_style_globals=True
# cython: binding=True
"""
Function pickling

REFERENCE: The python cookbook.
"""
import copyreg
import pickle
import sys
import types


def code_ctor(*args):
    """
    EXAMPLES:

    This indirectly tests this function. ::

        sage: def foo(a, b, c=10): return a+b+c
        sage: sage.misc.fpickle.reduce_code(foo.__code__)
        (<cyfunction code_ctor at ...>, ...)
        sage: unpickle_function(pickle_function(foo))
        <function foo at ...>
    """
    return types.CodeType(*args)


def reduce_code(co):
    """
    EXAMPLES::

        sage: def foo(N): return N+1
        sage: sage.misc.fpickle.reduce_code(foo.__code__)
        (<cyfunction code_ctor at ...>, ...)

    Test that the constructed code matches the original code::

        sage: ctor, args = sage.misc.fpickle.reduce_code(foo.__code__)
        sage: ctor(*args) == foo.__code__
        True
    """
    if co.co_freevars or co.co_cellvars:
        raise ValueError("Cannot pickle code objects from closures")

    co_args = (co.co_argcount,)
    co_args += (co.co_posonlyargcount,)
    co_args += (co.co_kwonlyargcount, co.co_nlocals,
                co.co_stacksize, co.co_flags, co.co_code,
                co.co_consts, co.co_names, co.co_varnames, co.co_filename,
                co.co_name)
    co_args += (co.co_qualname, co.co_firstlineno,
                co.co_linetable, co.co_exceptiontable)
    return (code_ctor, co_args)


copyreg.pickle(types.CodeType, reduce_code)


def pickle_function(func):
    """
    Pickle the Python function func.  This is not a normal pickle; you
    must use the unpickle_function method to unpickle the pickled
    function.

    NOTE: This does not work on all functions, but does work on
    'surprisingly' many functions.  In particular, it does not
    work on functions that includes nested functions.

    INPUT:

    - ``func`` -- a Python function

    OUTPUT: string

    EXAMPLES::

        sage: def f(N): return N+1
        ...
        sage: g = pickle_function(f)
        sage: h = unpickle_function(g)
        sage: h(10)
        11
    """
    return pickle.dumps(func.__code__)


def unpickle_function(pickled):
    """
    Unpickle a pickled function.

    EXAMPLES::

        sage: def f(N, M): return N*M
        ...
        sage: unpickle_function(pickle_function(f))(3,5)
        15
    """
    recovered = pickle.loads(pickled)
    return types.FunctionType(recovered, globals())


def call_pickled_function(fpargs):
    import sage.all as toplevel
    (fp, (args, kwds)) = fpargs
    f = eval("unpickle_function(fp)", toplevel.__dict__, {'fp': fp})
    res = eval("f(*args, **kwds)", toplevel.__dict__,
               {'args': args, 'kwds': kwds, 'f': f})
    return ((args, kwds), res)


# The following four methods are taken from twisted.persisted.styles - the
# import of twisted.persisted.styles takes a long time and we do not use
# most functionality it provides
def pickleMethod(method):
    'support function for copyreg to pickle method refs'

    if isinstance(method.__self__, type):
        # This is a class method, so get it from the type directly
        return (getattr, (method.__self__, method.__func__.__name__))
    else:
        cls = method.__self__.__class__
        return (unpickleMethod, (method.__func__.__name__, method.__self__, cls))


def unpickleMethod(im_name,
                   __self__,
                   im_class):
    'support function for copyreg to unpickle method refs'
    try:
        unbound = getattr(im_class, im_name)
        if __self__ is None:
            return unbound

        bound = types.MethodType(unbound, __self__)
        return bound
    except AttributeError:
        assert __self__ is not None, "No recourse: no instance to guess from."
        # Attempt a common fix before bailing -- if classes have
        # changed around since we pickled this method, we may still be
        # able to get it by looking on the instance's current class.
        unbound = getattr(__self__.__class__, im_name)

        bound = types.MethodType(unbound, __self__)
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
