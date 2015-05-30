"""
Deprecated attach module

TESTS::

    sage: import sage.misc.attached_files
    sage: sage.misc.attached_files.load_attach_mode()
    doctest:...: DeprecationWarning: 
    Importing load_attach_mode from here is deprecated. If you need to use it, please import it directly from sage.repl.attach
    See http://trac.sagemath.org/17396 for details.
    (False, True)
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.repl.attach', '*', deprecation=17396)
