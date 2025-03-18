r"""
Logging Backend

It records, for debugging and unit testing purposes, all calls to
backend methods in one of three ways.

See :class:`LoggingBackendFactory` for more information.
"""

# ****************************************************************************
#       Copyright (C) 2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.numerical.backends.generic_backend import GenericBackend


def _format_function_call(fn_name, *v, **k):
    """
    Return a Python function call as a string.

    EXAMPLES::

        sage: from sage.numerical.backends.logging_backend import _format_function_call
        sage: _format_function_call('foo', 17, hellooooo='goodbyeeee')
        "foo(17, hellooooo='goodbyeeee')"
    """
    args = [repr(a) for a in v] + ["%s=%r" % (arg, val) for arg, val in k.items()]
    return "{}({})".format(fn_name, ", ".join(args))


def _make_wrapper(backend, attr):
    """
    Return a wrapper for the backend method named by ``attr`` that does the logging.

    Documentation and other metadata is copied from the `backend` method named ``attr``.

    EXAMPLES::

        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: from sage.numerical.backends.logging_backend import _make_wrapper, LoggingBackend
        sage: backend = get_solver(solver='GLPK')
        sage: w = _make_wrapper(backend, 'ncols')
        sage: logging_backend = LoggingBackend(backend)
        sage: w(logging_backend)
        # p.ncols()
        # result: 0
        0
    """
    def m(self, *args, **kwdargs):
        funcall = _format_function_call("p." + attr, *args, **kwdargs)
        a = getattr(self._backend, attr)
        if self._printing:
            print("# {}".format(funcall))
        if self._doctest:
            self._doctest.write("        sage: {}\n".format(funcall))
        try:
            result = a(*args, **kwdargs)
        except Exception as e:
            if self._printing:
                print("# exception: {}".format(e))
            if self._doctest:
                self._doctest.write("        Traceback (most recent call last):\n"
                                    "        ...\n"
                                    "        MIPSolverException: {}\n".format(e))
            if self._test_method:
                self._test_method.write(("        with tester.assertRaises({}) as cm:\n" +
                                         "            {}\n").format(type(e).__name__, funcall))
            raise
        else:
            if self._printing:
                print("# result: {}".format(result))
            if self._doctest:
                self._doctest.write("        {}\n".format(result))
            if self._test_method:
                if result is None:
                    self._test_method.write("        tester.assertIsNone({})\n".format(funcall))
                elif type(result) is float:
                    # TODO: by default assertAlmostEqual does 7 decimal places (not significant digits)
                    # better perhaps to compute an appropriate 'places' or 'delta' parameter from result.
                    self._test_method.write("        tester.assertAlmostEqual({}, {})\n".format(funcall, result))
                else:
                    self._test_method.write("        tester.assertEqual({}, {})\n".format(funcall, result))
        return result
    from functools import update_wrapper
    update_wrapper(m, getattr(backend, attr))
    return m


class LoggingBackend(GenericBackend):
    """
    See :class:`LoggingBackendFactory` for documentation.

    EXAMPLES::

        sage: import sage.numerical.backends.logging_backend
        sage: from sage.numerical.backends.logging_backend import LoggingBackend
        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: b = get_solver(solver = "GLPK")
        sage: lb = LoggingBackend(backend=b)
        sage: lb.add_variable(obj=42, name='Helloooooo')
        # p.add_variable(obj=42, name='Helloooooo')
        # result: 0
        0
        sage: lb.add_variable(obj=1789)
        # p.add_variable(obj=1789)
        # result: 1
        1

    .. :no-undoc-members:
    """

    def __init__(self, backend, printing=True, doctest=None, test_method=None,
                 base_ring=None):
        """
        See :class:`LoggingBackendFactory` for documentation.

        EXAMPLES::

            sage: import sage.numerical.backends.logging_backend
            sage: from sage.numerical.backends.logging_backend import LoggingBackend
            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: b = get_solver(solver = "GLPK")
            sage: lb = LoggingBackend(backend=b)
        """
        self._backend = backend
        self._printing = printing
        self._doctest = doctest
        self._test_method = test_method
        self._base_ring = base_ring

    def __getattr__(self, attr):
        """
        Look up an attribute in the instance.

        It is provided to dynamically create delegating methods for all methods of
        the backend that are not part of the :class:`GenericBackend` interface,
        from which :class:`LoggingBackend` inherits.

        EXAMPLES::

            sage: import sage.numerical.backends.logging_backend
            sage: from sage.numerical.backends.logging_backend import LoggingBackend
            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: b = get_solver(solver = "GLPK")
            sage: lb = LoggingBackend(backend=b)
            sage: lb.print_ranges
            <bound method ...>
        """
        _a = getattr(self._backend, attr)
        if callable(_a):
            # make a bound method
            import types
            _mm = types.MethodType(_make_wrapper(self._backend, attr), self)
            # cache it
            setattr(self, attr, _mm)
            return _mm
        else:
            return _a

    def base_ring(self):
        """
        Return the base ring.

        The backend's base ring can be overridden.  It is best to run
        the tests with GLPK and override the base ring to ``QQ``.  Then
        default input to backend methods, prepared by
        :class:`MixedIntegerLinearProgram`, depends on the base ring.
        This way input will be rational and so suitable for both exact
        and inexact methods; whereas output will be float and will thus
        trigger :func:`assertAlmostEqual` tests.

        EXAMPLES::

            sage: import sage.numerical.backends.logging_backend
            sage: from sage.numerical.backends.logging_backend import LoggingBackend
            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: b = get_solver(solver = "GLPK")
            sage: lb = LoggingBackend(backend=b)
            sage: lb.base_ring()
            Real Double Field
            sage: from sage.rings.rational_field import QQ
            sage: lb = LoggingBackend(backend=b, base_ring=QQ)
            sage: lb.base_ring()
            Rational Field
        """
        if self._base_ring is not None:
            return self._base_ring
        else:
            return self._backend.base_ring()


# Override all methods that we inherited from GenericBackend
# by delegating methods
def _override_attr(attr):
    """
    Override a method by a delegating method.
    """
    a = getattr(LoggingBackend, attr)
    if callable(a):
        m = _make_wrapper(GenericBackend(), attr)
        setattr(LoggingBackend, attr, m)


for attr in dir(LoggingBackend):
    if not attr.startswith("_") and attr not in ("zero", "base_ring"):
        _override_attr(attr)

test_method_template = \
r'''
    @classmethod
    def _test_{name}(cls, tester=None, **options):
        """
        Run tests on ...

        TESTS::

            SAGE: from sage.numerical.backends.generic_backend import GenericBackend
            SAGE: p = GenericBackend()
            SAGE: p._test_{name}()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        p = cls()                         # fresh instance of the backend
        if tester is None:
            tester = p._tester(**options)
'''.replace("SAGE:", "sage:") # so that the above test does not get picked up by the doctester

from sage.rings.rational_field import QQ


def LoggingBackendFactory(solver=None, printing=True, doctest_file=None, test_method_file=None,
                          test_method=None, base_ring=QQ):
    """
    Factory that constructs a :class:`LoggingBackend` for debugging and testing.

    An instance of it can be passed as the solver argument of
    :func:`sage.numerical.backends.generic_backend.get_solver` and
    :class:`MixedIntegerLinearProgram`.

    EXAMPLES:

    Assume that we have the following function that does some
    computation using :class:`MixedIntegerLinearProgram` (or MIP
    backend methods), and suppose we have observed that it works with
    the GLPK backend, but not with the COIN backend::

        sage: def compute_something(solver='GLPK'):
        ....:     from sage.numerical.mip import MIPSolverException
        ....:     mip = MixedIntegerLinearProgram(solver=solver)
        ....:     lb = mip.get_backend()
        ....:     lb.add_variable(obj=42, name='Helloooooo')
        ....:     lb.add_variable(obj=1789)
        ....:     try:
        ....:         lb.solve()
        ....:     except MIPSolverException:
        ....:         return 4711
        ....:     else:
        ....:         return 91

    We can investigate what the backend methods are doing by running a
    :class:`LoggingBackend` in its in-terminal logging mode::

        sage: import sage.numerical.backends.logging_backend
        sage: from sage.numerical.backends.logging_backend import LoggingBackendFactory
        sage: compute_something(solver = LoggingBackendFactory(solver='GLPK'))
        # p = get_solver(solver='GLPK')
        # p.add_variable(obj=42, name='Helloooooo')
        # result: 0
        # p.add_variable(obj=1789)
        # result: 1
        # p.solve()
        # exception: GLPK: The LP (relaxation) problem has no dual feasible solution
        4711

    By replacing 'GLPK' by 'COIN' above, we can then compare the two
    logs and see where they differ.

    Imagine that we have now fixed the bug in the COIN backend, and we
    want to add a doctest that documents this fact.  We do not want to
    call ``compute_something`` in the doctest, but rather just have a
    sequence of calls to backend methods.

    We can have the doctest autogenerated by running a
    :class:`LoggingBackend` in its doctest-writing mode::

        sage: fname = tmp_filename()
        sage: compute_something(solver = LoggingBackendFactory(solver='GLPK', printing=False,
        ....:                                                  doctest_file=fname))
        4711
        sage: with open(fname) as f:
        ....:     for line in f.readlines(): _ = sys.stdout.write('|{}'.format(line))
        |        sage: p = get_solver(solver='GLPK')
        |        sage: p.add_variable(obj=42, name='Helloooooo')
        |        0
        |        sage: p.add_variable(obj=1789)
        |        1
        |        sage: p.solve()
        |        Traceback (most recent call last):
        |        ...
        |        MIPSolverException: GLPK: The LP (relaxation) problem has no dual feasible solution

    We then copy from the generated file and paste into the source
    code of the COIN backend.

    If this test seems valuable enough that all backends should be
    tested against it, we should create a test method instead of a
    docstring.

    We can have the test method autogenerated by running a
    :class:`LoggingBackend` in its test-method-writing mode::

        sage: fname = tmp_filename()
        sage: compute_something(solver= LoggingBackendFactory(solver='GLPK', printing=False,
        ....:                                                 test_method_file=fname,
        ....:                                                 test_method='something'))
        4711
        sage: with open(fname) as f:
        ....:     for line in f.readlines(): _ = sys.stdout.write('|{}'.format(line))
        |
        |    @classmethod
        |    def _test_something(cls, tester=None, **options):
        |        ...
        |        Run tests on ...
        |
        |        TESTS::
        |
        |            sage: from sage.numerical.backends.generic_backend import GenericBackend
        |            sage: p = GenericBackend()
        |            sage: p._test_something()
        |            Traceback (most recent call last):
        |            ...
        |            NotImplementedError
        |        ...
        |        p = cls()                         # fresh instance of the backend
        |        if tester is None:
        |            tester = p._tester(**options)
        |        tester.assertEqual(p.add_variable(obj=42, name='Helloooooo'), 0)
        |        tester.assertEqual(p.add_variable(obj=1789), 1)
        |        with tester.assertRaises(MIPSolverException) as cm:
        |            p.solve()

    We then copy from the generated file and paste into the source
    code of the generic backend, where all test methods are defined.

    If ``test_method_file`` is not provided, a default output file name
    will be computed from ``test_method``.
    """

    if test_method is not None:
        if test_method_file is None:
            # Construct output file name from method name.
            test_method_file = "test_{}.py".format(test_method)
    else:
        test_method = 'CHANGE'
        # Will have to be edited by user in generated file.

    if doctest_file is not None:
        doctest = open(doctest_file, "w", 1)  # line-buffered
    else:
        doctest = None
    if test_method_file is not None:
        test_method_output = open(test_method_file, "w", 1)  # line-buffered
    else:
        test_method_output = None

    def logging_solver(**kwds):
        """
        Create an instance of :class:`LoggingBackend`.
        """
        construct = "p = {}".format(_format_function_call('get_solver', solver=solver, **kwds))
        if printing:
            print("# {}".format(construct))
        if doctest is not None:
            doctest.write("        sage: {}\n".format(construct))
        if test_method_output is not None:
            test_method_output.write(test_method_template.format(name=test_method))
        from sage.numerical.backends.generic_backend import get_solver
        return LoggingBackend(backend=get_solver(solver=solver, **kwds),
                              printing=printing, doctest=doctest, test_method=test_method_output,
                              base_ring=base_ring)

    return logging_solver
