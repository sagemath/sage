r"""
Decorate interface for parallel computation
"""

import types

from sage.rings.integer import Integer

from .reference import parallel_iter as p_iter_reference
from .use_fork import p_iter_fork
from . import multiprocessing_sage
from sage.misc.instancedoc import instancedoc


def normalize_input(a):
    r"""
    Convert ``a`` to a pair ``(args, kwds)`` using some rules:

    - if already of that form, leave that way.
    - if ``a`` is a tuple make ``(a,{})``
    - if ``a`` is a dict make ``(tuple(),a)``
    - otherwise make ``((a,),{})``

    INPUT:

    - ``a`` -- object

    OUTPUT:

    - ``args`` -- tuple
    - ``kwds`` -- dictionary

    EXAMPLES::

        sage: sage.parallel.decorate.normalize_input( (2, {3:4}) )
        ((2, {3: 4}), {})
        sage: sage.parallel.decorate.normalize_input( (2,3) )
        ((2, 3), {})
        sage: sage.parallel.decorate.normalize_input( {3:4} )
        ((), {3: 4})
        sage: sage.parallel.decorate.normalize_input( 5 )
        ((5,), {})
    """
    if isinstance(a, tuple) and len(a) == 2 and isinstance(a[0], tuple) and isinstance(a[1], dict):
        return a
    elif isinstance(a, tuple):
        return (a, {})
    elif isinstance(a, dict):
        return (tuple(), a)
    else:
        return ((a,), {})


class Parallel:
    r"""
    Create a ``parallel``-decorated function.
    This is the object created by :func:`parallel`.
    """
    def __init__(self, p_iter='fork', ncpus=None, **kwds):
        """
        EXAMPLES::

            sage: P = sage.parallel.decorate.Parallel(); P
            <sage.parallel.decorate.Parallel object at 0x...>
        """
        # The default p_iter is currently the 'fork' implementation.
        # This has changed.

        self.p_iter = None

        if isinstance(p_iter, (int, Integer)):
            p_iter, ncpus = 'fork', p_iter

        if ncpus is None:
            from .ncpus import ncpus as compute_ncpus
            ncpus = compute_ncpus()

        if p_iter == 'fork':
            self.p_iter = p_iter_fork(ncpus, **kwds)
        elif p_iter == 'multiprocessing':
            self.p_iter = multiprocessing_sage.pyprocessing(ncpus)
        elif p_iter == 'reference':
            self.p_iter = p_iter_reference
        elif isinstance(p_iter, str):
            raise ValueError("unknown iterator '%s'" % p_iter)
        else:
            if self.p_iter is None:
                self.p_iter = p_iter

    def __call__(self, f):
        r"""
        Create a callable object that wraps ``f`` and that when called
        with a list of inputs returns an iterator over pairs ``(x,
        f(x))`` in possibly random order. Here ``x`` is replaced by
        its normalized form ``(args, kwds)`` using
        :func:`normalize_inputs`.

        INPUT:

        - ``f`` -- Python callable object or function

        OUTPUT: decorated version of ``f``

        EXAMPLES::

            sage: from sage.parallel.decorate import Parallel
            sage: p = Parallel()
            sage: f = x^2 - 1                                                           # needs sage.symbolic
            sage: p(f)                                                                  # needs sage.symbolic
            <sage.parallel.decorate.ParallelFunction object at ...>

            sage: P = sage.parallel.decorate.Parallel()
            sage: def g(n,m): return n+m
            sage: h = P(g)          # indirect doctest
            sage: list(h([(2,3)]))
            [(((2, 3), {}), 5)]
        """
        return ParallelFunction(self, f)


@instancedoc
class ParallelFunction:
    """
    Class which parallelizes a function or class method.
    This is typically accessed indirectly through
    :meth:`Parallel.__call__`.
    """
    def __init__(self, parallel, func):
        """
        .. NOTE::

            This is typically accessed indirectly through
            :meth:`Parallel.__call__`.

        INPUT:

        - ``parallel`` -- a :class:`Parallel` object which controls
          how the parallel execution will be done

        - ``func`` -- Python callable object or function
        """
        self.parallel = parallel
        self.func = func

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: from sage.parallel.decorate import Parallel
            sage: p = Parallel()
            sage: def f(x):
            ....:     return x*x
            sage: pf = p(f); pf
            <sage.parallel.decorate.ParallelFunction object at ...>
            sage: pf(2)
            4
            sage: sorted(pf([2,3]))
            [(((2,), {}), 4), (((3,), {}), 9)]
        """
        if len(args) > 0 and isinstance(args[0], (list, types.GeneratorType)):
            return self.parallel.p_iter(self.func, (normalize_input(a)
                                                    for a in args[0]))
        return self.func(*args, **kwds)

    def __get__(self, instance, owner):
        """
        Implement part of the descriptor protocol for
        :class:`ParallelFunction` objects.

        .. NOTE::

            This is the key to fixing :issue:`11461`.

        EXAMPLES:

        We verify that the decorated functions work correctly on
        methods, classmethods, and staticmethods, for both the
        parallel and non-parallel versions::

            sage: class Foo():
            ....:     @parallel(2)
            ....:     def square(self, n):
            ....:         return n*n
            ....:     @parallel(2)
            ....:     @classmethod
            ....:     def square_classmethod(cls, n):
            ....:         return n*n
            ....:     @parallel(2)
            ....:     @staticmethod
            ....:     def square_staticmethod(n):
            ....:         return n*n
            sage: a = Foo()
            sage: a.square(3)
            9
            sage: sorted(a.square([2,3]))
            [(((2,), {}), 4), (((3,), {}), 9)]
            sage: a.square_classmethod(3)
            9
            sage: sorted(a.square_classmethod([2,3]))
            [(((2,), {}), 4), (((3,), {}), 9)]
            sage: Foo.square_classmethod(3)
            9
            sage: sorted(Foo.square_classmethod([2,3]))
            [(((2,), {}), 4), (((3,), {}), 9)]
            sage: a.square_staticmethod(3)
            9
            sage: sorted(a.square_staticmethod([2,3]))
            [(((2,), {}), 4), (((3,), {}), 9)]
            sage: Foo.square_staticmethod(3)
            9
            sage: sorted(Foo.square_staticmethod([2,3]))
            [(((2,), {}), 4), (((3,), {}), 9)]
        """
        try:
            # If this ParallelFunction object is accessed as an
            # attribute of a class or instance, the underlying function
            # should be "accessed" in the same way.
            new_func = self.func.__get__(instance, owner)
        except AttributeError:
            # This will happen if a non-function attribute is
            # decorated.  For example, an expression that's an
            # attribute of a class.
            new_func = self.func
        return ParallelFunction(self.parallel, new_func)

    def _sage_argspec_(self):
        """
        Return the argument specification for this object, which is
        just the argument specification for the underlying function.
        See :mod:`sage.misc.sageinspect` for more information on
        this convention.

        EXAMPLES::

            sage: from sage.parallel.decorate import Parallel
            sage: p = Parallel(2)
            sage: def f(x, y):
            ....:     return x + y
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(p(f))
            FullArgSpec(args=['x', 'y'], varargs=None, varkw=None, defaults=None,
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
        """
        from sage.misc.sageinspect import sage_getargspec
        return sage_getargspec(self.func)

    def _sage_src_(self):
        """
        Return the source code for this object, which is just the
        source code for the underlying function.  See
        :mod:`sage.misc.sageinspect` for more information on this
        convention.

        EXAMPLES::

            sage: from sage.parallel.decorate import Parallel
            sage: p = Parallel(2)
            sage: def f(x, y):
            ....:     return x + y
            sage: from sage.misc.sageinspect import sage_getsource
            sage: 'return x + y' in sage_getsource(p(f))
            True
        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.func)

    def _instancedoc_(self):
        r"""
        Return the docstring for this object, which is just the
        docstring for the underlying function.  See
        :mod:`sage.misc.sageinspect` for more information on this
        convention.

        EXAMPLES::

            sage: from sage.parallel.decorate import Parallel
            sage: p = Parallel(2)
            sage: def f(x, y):
            ....:     '''Test docstring'''
            ....:     return x + y
            sage: from sage.misc.sageinspect import sage_getdoc
            sage: sage_getdoc(p(f))
            'Test docstring\n'
        """
        return self.func.__doc__


def parallel(p_iter='fork', ncpus=None, **kwds):
    r"""
    This is a decorator that gives a function a parallel interface,
    allowing it to be called with a list of inputs, whose values will
    be computed in parallel.

    .. warning::

         The parallel subprocesses will not have access to data
         created in pexpect interfaces.  This behavior with respect to
         pexpect interfaces is very important to keep in mind when
         setting up certain computations.  It's the one big limitation
         of this decorator.

    INPUT:

    - ``p_iter`` -- parallel iterator function or string:
      - ``'fork'`` -- (default) use a new forked subprocess for each input
      - ``'multiprocessing'`` -- use multiprocessing library
      - ``'reference'`` -- use a fake serial reference implementation
    - ``ncpus`` -- integer; maximal number of subprocesses to use at the same time
    - ``timeout`` -- number of seconds until each subprocess is killed (only supported
      by ``'fork'``; zero means not at all)

    .. warning::

         If you use anything but ``'fork'`` above, then a whole new
         subprocess is spawned, so none of your local state (variables,
         certain functions, etc.) is available.

    EXAMPLES:

    We create a simple decoration for a simple function.  The number
    of cpus (or cores, or hardware threads) is automatically detected::

        sage: @parallel
        ....: def f(n): return n*n
        sage: f(10)
        100
        sage: sorted(list(f([1,2,3])))
        [(((1,), {}), 1), (((2,), {}), 4), (((3,), {}), 9)]

    We use exactly two cpus::

        sage: @parallel(2)
        ....: def f(n): return n*n


    We create a decorator that uses three subprocesses, and times out
    individual processes after 10 seconds::

        sage: @parallel(ncpus=3, timeout=10)
        ....: def fac(n): return factor(2^n-1)
        sage: for X, Y in sorted(list(fac([101,119,151,197,209]))): print((X,Y))        # needs sage.libs.pari
        (((101,), {}), 7432339208719 * 341117531003194129)
        (((119,), {}), 127 * 239 * 20231 * 131071 * 62983048367 * 131105292137)
        (((151,), {}), 18121 * 55871 * 165799 * 2332951 * 7289088383388253664437433)
        (((197,), {}), 7487 * 26828803997912886929710867041891989490486893845712448833)
        (((209,), {}), 23 * 89 * 524287 * 94803416684681 * 1512348937147247 * 5346950541323960232319657)

        sage: @parallel('multiprocessing')
        ....: def f(N): return N^2
        sage: v = list(f([1,2,4])); v.sort(); v
        [(((1,), {}), 1), (((2,), {}), 4), (((4,), {}), 16)]
        sage: @parallel('reference')
        ....: def f(N): return N^2
        sage: v = list(f([1,2,4])); v.sort(); v
        [(((1,), {}), 1), (((2,), {}), 4), (((4,), {}), 16)]

    For functions that take multiple arguments, enclose the arguments in tuples
    when calling the parallel function::

        sage: @parallel
        ....: def f(a,b): return a*b
        sage: for X, Y in sorted(list(f([(2,3),(3,5),(5,7)]))): print((X, Y))
        (((2, 3), {}), 6)
        (((3, 5), {}), 15)
        (((5, 7), {}), 35)

    For functions that take a single tuple as an argument, enclose it in an
    additional tuple at call time, to distinguish it as the first argument,
    as opposed to a tuple of arguments::

        sage: @parallel
        ....: def firstEntry(aTuple): return aTuple[0]
        sage: for X, Y in sorted(list(firstEntry([((1,2,3,4),),((5,6,7,8),)]))): print((X, Y))
        ((((1, 2, 3, 4),), {}), 1)
        ((((5, 6, 7, 8),), {}), 5)

    The parallel decorator also works with methods, classmethods, and
    staticmethods.  Be sure to apply the parallel decorator after ("above")
    either the ``classmethod`` or ``staticmethod`` decorators::

        sage: class Foo():
        ....:     @parallel(2)
        ....:     def square(self, n):
        ....:         return n*n
        ....:     @parallel(2)
        ....:     @classmethod
        ....:     def square_classmethod(cls, n):
        ....:         return n*n
        sage: a = Foo()
        sage: a.square(3)
        9
        sage: sorted(a.square([2,3]))
        [(((2,), {}), 4), (((3,), {}), 9)]
        sage: Foo.square_classmethod(3)
        9
        sage: sorted(Foo.square_classmethod([2,3]))
        [(((2,), {}), 4), (((3,), {}), 9)]
        sage: Foo.square_classmethod(3)
        9

    .. warning::

       Currently, parallel methods do not work with the
       multiprocessing implementation.
    """
    if isinstance(p_iter, types.FunctionType):
        return Parallel()(p_iter)
    return Parallel(p_iter, ncpus, **kwds)


###################################################################
# The @fork decorator -- evaluate a function with no side effects
# in memory, so the only side effects (if any) are on disk.
#
# We have both a function and a class below, so that the decorator
# can be used with or without options:
#
#   @fork
#   def f(...): ...
# and
#   @fork(...options...):
#   def f(...): ...
###################################################################

class Fork:
    """
    A ``fork`` decorator class.
    """
    def __init__(self, timeout=0, verbose=False):
        """
        INPUT:

        - ``timeout`` -- (default: 0) kill the subprocess after it has run this
          many seconds (wall time), or if ``timeout`` is zero, do not kill it.
        - ``verbose`` -- boolean (default: ``False``); whether to print anything about
          what the decorator does (e.g., killing the subprocess)

        EXAMPLES::

            sage: sage.parallel.decorate.Fork()
            <sage.parallel.decorate.Fork object at 0x...>
            sage: sage.parallel.decorate.Fork(timeout=3)
            <sage.parallel.decorate.Fork object at 0x...>
        """
        self.timeout = timeout
        self.verbose = verbose

    def __call__(self, f):
        """
        INPUT:

        - ``f`` -- a function

        OUTPUT: a decorated function

        EXAMPLES::

            sage: F = sage.parallel.decorate.Fork(timeout=3)
            sage: def g(n,m): return n+m
            sage: h = F(g)     # indirect doctest
            sage: h(2,3)
            5
        """
        P = Parallel(p_iter='fork', ncpus=1, timeout=self.timeout,
                     verbose=self.verbose)
        g = P(f)

        def h(*args, **kwds):
            return list(g([(args, kwds)]))[0][1]
        return h


def fork(f=None, timeout=0, verbose=False):
    """
    Decorate a function so that when called it runs in a forked subprocess.

    This means that it will not have any in-memory side effects on the
    parent Sage process. The pexpect interfaces are all reset.

    INPUT:

    - ``f`` -- a function
    - ``timeout`` -- (default: 0) if positive, kill the subprocess after
      this many seconds (wall time)
    - ``verbose`` -- boolean (default: ``False``); whether to print anything
      about what the decorator does (e.g., killing the subprocess)

    .. warning::

        The forked subprocess will not have access to data created
        in pexpect interfaces.  This behavior with respect to pexpect
        interfaces is very important to keep in mind when setting up
        certain computations.  It's the one big limitation of this
        decorator.

    EXAMPLES:

    We create a function and run it with the ``fork`` decorator.  Note
    that it does not have a side effect.  Despite trying to change
    the global variable ``a`` below in ``g``, the variable ``a`` does not
    get changed::

        sage: a = 5
        sage: @fork
        ....: def g(n, m):
        ....:     global a
        ....:     a = 10
        ....:     return factorial(n).ndigits() + m
        sage: g(5, m=5)
        8
        sage: a
        5

    We use ``fork`` to make sure that the function terminates after 100 ms,
    no matter what::

        sage: @fork(timeout=0.1, verbose=True)
        ....: def g(n, m): return factorial(n).ndigits() + m
        sage: g(10^7, m=5)
        Killing subprocess ... with input ((10000000,), {'m': 5}) which took too long
        'NO DATA (timed out)'

    We illustrate that the state of the pexpect interface is not altered by
    forked functions (they get their own new pexpect interfaces!)::

        sage: # needs sage.libs.pari
        sage: gp.eval('a = 5')
        '5'
        sage: @fork()
        ....: def g():
        ....:     gp.eval('a = 10')
        ....:     return gp.eval('a')
        sage: g()
        '10'
        sage: gp.eval('a')
        '5'

    We illustrate that the forked function has its own pexpect
    interface::

        sage: gp.eval('a = 15')                                                         # needs sage.libs.pari
        '15'
        sage: @fork()
        ....: def g(): return gp.eval('a')
        sage: g()                                                                       # needs sage.libs.pari
        'a'

    We illustrate that segfaulting subprocesses are no trouble at all::

        sage: cython('def f(): print(<char*>0)')                                        # needs sage.misc.cython
        sage: @fork
        ....: def g():
        ....:     os.environ["CYSIGNALS_CRASH_NDEBUG"]="yes" # skip enhanced backtrace (it is slow)
        ....:     f()
        sage: print("this works"); g()                                                  # needs sage.misc.cython
        this works...
        <BLANKLINE>
        ------------------------------------------------------------------------
        Unhandled SIG...
        ------------------------------------------------------------------------
        'NO DATA'
    """
    F = Fork(timeout=timeout, verbose=verbose)
    return F(f) if f else F
