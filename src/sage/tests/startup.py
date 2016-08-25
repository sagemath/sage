r"""
Ensure that certain modules are not loaded on startup.

EXAMPLES::

    sage: 'numpy' in sys.modules
    False
    sage: 'sage.libs.gap.libgap' in sys.modules
    False

Check that IPython is not imported at startup (:trac:`18726`). It is
imported by the doctest framework, so the simple test like above would
not work. Instead, we test this by starting a new Python process::

    sage: from sage.tests.cmdline import test_executable
    sage: cmd = "from sage.all import *\nprint('IPython' in sys.modules)\n"
    sage: print(test_executable(["sage", "--python"], cmd)[0])  # long time
    False
"""
