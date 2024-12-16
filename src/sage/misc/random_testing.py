r"""
Random testing

Some Sage modules do random testing in their doctests; that is, they
construct test cases using a random number generator.  To get the
broadest possible test coverage, we want everybody who runs the
doctests to use a different random seed; but we also want to be able
to reproduce the problems when debugging.  This module provides a
decorator to help write random testers that meet these goals.
"""

from functools import wraps


def random_testing(fn):
    r"""
    This decorator helps create random testers.  These can be run as
    part of the standard Sage test suite; everybody who runs the test
    will use a different random number seed, so many different random
    tests will eventually be run.

    INPUT:

    - ``fn`` -- the function that we are wrapping for random testing

    The resulting function will take two additional arguments, *seed*
    (default: ``None``) and *print_seed* (default: ``False``).  The
    result will set the random number seed to the given seed value (or
    to a truly random value, if *seed* is not specified), then call
    the original function.  If *print_seed* is true, then the seed will
    be printed before calling the original function.  If the original
    function raises an exception, then the random seed that was used
    will be displayed, along with a message entreating the user to
    submit a bug report.  All other arguments will be passed through
    to the original function.

    Here is a set of recommendations for using this wrapper.

    The function to be tested should take arguments specifying the
    difficulty of the test (size of the test cases, number of
    iterations, etc.), as well as an argument *verbose* (defaulting to
    false).  With *verbose* true, it should print the values being
    tested.  Suppose ``test_foo()`` takes an argument for number of
    iterations.  Then the doctests could be::

        test_foo(2, verbose=True, seed=0)
        test_foo(10)
        test_foo(100) # long time

    The first doctest, with the specified seed and ``verbose=True``, simply
    verifies that the tests really are reproducible (that ``test_foo``
    is correctly using the :mod:`randstate` framework).  The next two tests
    use truly random seeds, and will print out the seed used if the test
    fails (raises an exception).

    If you want a very long-running test using this setup, you should do
    something like::

        for _ in range(10^10): test_foo(100)

    instead of::

        test_foo(10^12)

    If the test fails after several hours, the latter snippet would
    make you rerun the test for several hours while reproducing and
    debugging the problem.  With the former snippet, you only need to
    rerun ``test_foo(100)`` with a known-failing random seed.

    See :func:`sage.misc.random_testing.check_add_commutes` for a
    simple example using this decorator, and :mod:`sage.rings.tests`
    for realistic uses.

    Setting *print_seed* to true is useless in doctests, because the
    random seed printed will never match the expected doctest result
    (and using ``# random`` means the doctest framework will never
    report an error even if one happens).  However, it is useful if
    you have a random test that sometimes segfaults.  The normal
    print-the-random-seed-on-exceptions won't work then, so you can
    run::

        while True: test_foo(print_seed=True)

    and look at the last seed that was printed before it crashed.


    TESTS::

        sage: from sage.misc.random_testing import random_testing
        sage: def foo(verbose=False):
        ....:     'oh look, a docstring'
        ....:     n = ZZ.random_element(2^50)
        ....:     if verbose:
        ....:         print("Random value: %s" % n)
        ....:     assert(n == 49681376900427)
        sage: foo = random_testing(foo)
        sage: foo(seed=0, verbose=True)
        Random value: 49681376900427
        sage: foo(seed=15, verbose=True)
        Random value: 1049538412064764
        Random testing has revealed a problem in foo
        Please report this bug!  You may be the first
        person in the world to have seen this problem.
        Please include this random seed in your bug report:
        Random seed: 15
        AssertionError()
        sage: foo() # random
        Random testing has revealed a problem in foo
        Please report this bug!  You may be the first
        person in the world to have seen this problem.
        Please include this random seed in your bug report:
        Random seed: 272500700755151445506092479579811710040
        AssertionError()
        sage: foo.__doc__
        'oh look, a docstring'
        sage: foo.__name__
        'foo'
        sage: def bar(): pass
        sage: bar = random_testing(bar)
        sage: bar(print_seed=True) # random
        Random seed: 262841091890156346923539765543814146051
    """
    from sage.misc.randstate import seed, initial_seed
    from sys import stdout

    @wraps(fn)
    def wrapped_fun(*args, **kwargs):
        arg_seed = None
        if 'seed' in kwargs:
            arg_seed = kwargs['seed']
            del kwargs['seed']
        with seed(arg_seed):
            used_seed = initial_seed()
            if 'print_seed' in kwargs:
                if kwargs['print_seed']:
                    print("Random seed: {}".format(used_seed))
                    del kwargs['print_seed']
                # I don't know if this line is necessary, but it can't
                # hurt; and it would be a real pity to lose the
                # information you need to reproduce a segfault because
                # it was missing...
                stdout.flush()
            try:
                fn(*args, **kwargs)
            except Exception as e:
                # We treat any sort of Exception as a doctest
                # failure.  (We have to eat the exception, because if
                # doctesting sees an exception, it doesn't display
                # whatever was printed before the exception happened
                # -- so the text we print here would be lost.)  Note
                # that KeyboardInterrupt is not an Exception, so
                # pressing Control-C doesn't print this message.
                print("Random testing has revealed a problem in " + fn.__name__)
                print("Please report this bug!  You may be the first")
                print("person in the world to have seen this problem.")
                print("Please include this random seed in your bug report:")
                print("Random seed: {}".format(used_seed))
                print(repr(e))
    return wrapped_fun


@random_testing
def check_add_commutes(trials, verbose=False):
    r"""
    This is a simple demonstration of the :func:`random_testing` decorator and
    its recommended usage.

    We test that addition is commutative over rationals.

    EXAMPLES::

        sage: from sage.misc.random_testing import check_add_commutes
        sage: check_add_commutes(2, verbose=True, seed=0)
        a == -4, b == 0 ...
        Passes!
        a == -1/2, b == -1/95 ...
        Passes!
        sage: check_add_commutes(10)
        sage: check_add_commutes(1000) # long time
    """
    from sage.rings.rational_field import QQ
    for _ in range(trials):
        a = QQ.random_element()
        b = QQ.random_element()
        if verbose:
            print("a == {}, b == {} ...".format(a, b))
        assert a + b == b + a
        if verbose:
            print("Passes!")


@random_testing
def check_add_is_mul(trials, verbose=False):
    r"""
    This example demonstrates a failing :func:`random_testing` test,
    and shows how to reproduce the error.

    DO NOT USE THIS AS AN EXAMPLE OF HOW TO USE
    :func:`random_testing`!  Instead, look at
    :func:`sage.misc.random_testing.check_add_commutes`.

    We test that ``a+b == a*b``, for *a*, *b* rational.  This is of
    course false, so the test will almost always fail.

    EXAMPLES::

        sage: from sage.misc.random_testing import test_add_is_mul

    We start by testing that we get reproducible results when setting
    *seed* to 0.

    ::

        sage: test_add_is_mul(2, verbose=True, seed=0)
        a == -4, b == 0 ...
        Random testing has revealed a problem in test_add_is_mul
        Please report this bug!  You may be the first
        person in the world to have seen this problem.
        Please include this random seed in your bug report:
        Random seed: 0
        AssertionError()

    Normally in a ``@random_testing`` doctest, we would leave off the
    ``verbose=True`` and the ``# random``.  We put it in here so that we can
    verify that we are seeing the exact same error when we reproduce
    the error below.

    ::

        sage: test_add_is_mul(10, verbose=True) # random
        a == -2/7, b == 1 ...
        Random testing has revealed a problem in test_add_is_mul
        Please report this bug!  You may be the first
        person in the world to have seen this problem.
        Please include this random seed in your bug report:
        Random seed: 216390410596009428782506007128692114173
        AssertionError()

    OK, now assume that some user has reported a
    :func:`test_add_is_mul` failure.  We can specify the same
    *random_seed* that was found in the bug report, and we will get the
    exact same failure so that we can debug the "problem".

    ::

        sage: test_add_is_mul(10, verbose=True, seed=216390410596009428782506007128692114173)
        a == -2/7, b == 1 ...
        Random testing has revealed a problem in test_add_is_mul
        Please report this bug!  You may be the first
        person in the world to have seen this problem.
        Please include this random seed in your bug report:
        Random seed: 216390410596009428782506007128692114173
        AssertionError()
    """
    from sage.rings.rational_field import QQ
    for _ in range(trials):
        a = QQ.random_element()
        b = QQ.random_element()
        if verbose:
            print("a == {}, b == {} ...".format(a, b))
        assert a + b == a * b
        if verbose:
            print("Passes!")
