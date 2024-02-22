"""
CPU Detection
"""

import os


def ncpus():
    """
    Detects the number of available CPUs in the system.

    ALGORITHM: :func:`os.cpu_count`

    EXAMPLES::

        sage: sage.parallel.ncpus.ncpus()  # random output -- depends on machine
        2
    """
    # Support Sage environment variable SAGE_NUM_THREADS
    # NOTE: while doctesting, this is forced to be 2 by the
    # sage-runtests script
    try:
        n = os.environ["SAGE_NUM_THREADS"]
    except KeyError:
        pass
    else:
        return int(n)

    return os.cpu_count() or 1
