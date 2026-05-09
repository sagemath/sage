Test combining various modifiers::

    sage: sys.maxsize  # long time, abs tol 0.001, needs 32_bit
    2147483646.999
    sage: sys.maxsize  # long time, abs tol 0.001, needs !32_bit
    9223372036854775806.999
