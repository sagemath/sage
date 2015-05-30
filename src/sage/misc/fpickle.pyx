"""
Function pickling

REFERENCE: The python cookbook.
"""

import types, copy_reg, cPickle

def code_ctor(*args):
    """
    EXAMPLES:

    This indirectly tests this function. ::

        sage: def foo(a,b,c=10): return a+b+c
        sage: sage.misc.fpickle.reduce_code(foo.__code__)
        (<built-in function code_ctor>, ...)
        sage: unpickle_function(pickle_function(foo))
        <function foo at ...>
    """
    return types.CodeType(*args)

def reduce_code(co):
    """
    EXAMPLES::

        sage: def foo(N): return N+1
        sage: sage.misc.fpickle.reduce_code(foo.__code__)
        (<built-in function code_ctor>, ...)
    """
    if co.co_freevars or co.co_cellvars:
        raise ValueError("Cannot pickle code objects from closures")
    return code_ctor, (co.co_argcount, co.co_nlocals, co.co_stacksize,
                       co.co_flags, co.co_code, co.co_consts, co.co_names,
                       co.co_varnames, co.co_filename, co.co_name,
                       co.co_firstlineno, co.co_lnotab)

copy_reg.pickle(types.CodeType, reduce_code)

def pickle_function(func):
    """
    Pickle the Python function func.  This is not a normal pickle; you
    must use the unpickle_function method to unpickle the pickled
    function.

    NOTE: This does not work on all functions, but does work on
    'surprisingly' many functions.  In particular, it does not
    work on functions that includes nested functions.

    INPUT:

        func -- a Python function

    OUTPUT:

        a string

    EXAMPLES::

        sage: def f(N): return N+1
        ...
        sage: g = pickle_function(f)
        sage: h = unpickle_function(g)
        sage: h(10)
        11
    """
    return cPickle.dumps(func.__code__)

def unpickle_function(pickled):
    """
    Unpickle a pickled function.

    EXAMPLES:

        sage: def f(N,M): return N*M
        ...
        sage: unpickle_function(pickle_function(f))(3,5)
        15
    """
    recovered = cPickle.loads(pickled)
    return types.FunctionType(recovered, globals())



def call_pickled_function(fpargs):
    import sage.all
    from sage.misc.fpickle import unpickle_function
    (fp, (args, kwds)) = fpargs
    f = eval("unpickle_function(fp)", sage.all.__dict__, {'fp':fp})
    res = eval("f(*args, **kwds)",sage.all.__dict__, {'args':args, 'kwds':kwds, 'f':f})
    return ((args, kwds), res)

# The following four methods are taken from twisted.persisted.styles - the
# import of twisted.persisted.styles takes a long time and we do not use
# most functionality it provides
def pickleMethod(method):
    'support function for copy_reg to pickle method refs'
    return unpickleMethod, (method.__func__.__name__,
                             method.__self__,
                             method.im_class)

def unpickleMethod(im_name,
                    __self__,
                    im_class):
    'support function for copy_reg to unpickle method refs'
    try:
        unbound = getattr(im_class,im_name)
        if __self__ is None:
            return unbound
        bound=types.MethodType(unbound.__func__,
                                 __self__)
        return bound
    except AttributeError:
        assert __self__ is not None,"No recourse: no instance to guess from."
        # Attempt a common fix before bailing -- if classes have
        # changed around since we pickled this method, we may still be
        # able to get it by looking on the instance's current class.
        unbound = getattr(__self__.__class__,im_name)
        if __self__ is None:
            return unbound
        bound=types.MethodType(unbound.__func__,
                                 __self__)
        return bound

copy_reg.pickle(types.MethodType,
                pickleMethod,
                unpickleMethod)

oldModules = {}

def pickleModule(module):
    'support function for copy_reg to pickle module refs'
    return unpickleModule, (module.__name__,)

def unpickleModule(name):
    'support function for copy_reg to unpickle module refs'
    if name in oldModules:
        name = oldModules[name]
    return __import__(name,{},{},'x')

copy_reg.pickle(types.ModuleType,
                pickleModule,
                unpickleModule)
