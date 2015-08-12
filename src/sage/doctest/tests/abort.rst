This causes an "Unhandled SIGABRT..." which should be handled by
the doctester::

    sage: from sage.ext.interrupt.tests import *
    sage: unguarded_abort()
