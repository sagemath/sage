This should raise a ``MemoryError``::

    sage: 256 ^ (2000 << 20)  # optional - memlimit
