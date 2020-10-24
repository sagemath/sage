This causes an "Unhandled SIGABRT..." which should be handled by
the doctester::

    sage: from cysignals.tests import unguarded_abort
    sage: unguarded_abort()
