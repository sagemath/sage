r"""
TESTS:

Sage integers can index NumPy matrices (see :trac:`10928`)::

    sage: import numpy as np
    sage: m = np.matrix(np.arange(4).reshape(2, 2))
    sage: a = m[:, int(0)]
    sage: b = m[:, Integer(0)]
    sage: a.shape, b.shape
    ((2, 1), (2, 1))
    sage: print('{}\n{}\n{}'.format(m, a, b))
    [[0 1]
     [2 3]]
    [[0]
     [2]]
    [[0]
     [2]]
"""
