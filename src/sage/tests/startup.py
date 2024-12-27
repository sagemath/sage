r"""
Ensure that certain modules are not loaded on startup.

Check that IPython is not imported at startup (:issue:`18726`). It is
imported by the doctest framework, so the simple test like above would
not work. Instead, we test this by starting a new Python process::

    sage: from sage.tests.cmdline import check_executable
    sage: environment = "sage.all"
    sage: cmd = f"from {environment} import *\nprint('IPython' in sys.modules)\n"
    sage: print(check_executable(["sage", "--python"], cmd)[0])  # long time
    False

Check that numpy (:issue:`11714`) and pyparsing are not imported on startup
as they increase the startup time. Since :issue:`23696` those are imported
by the doctest framework via a matplotlib import. Again the simple test
would not work (but we don't have to avoid loading IPython)::

    sage: from sage.tests.cmdline import check_executable
    sage: cmd = "print('numpy' in sys.modules)\n"
    sage: print(check_executable(["sage", "-c", cmd])[0])  # long time
    False
    sage: cmd = "print('pyparsing' in sys.modules)\n"
    sage: print(check_executable(["sage", "-c", cmd])[0])  # long time
    False
"""
