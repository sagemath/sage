These doctests are intentionally failing. They test that, even though
a tolerance is specified, not everything is accepted.

Tolerance but no numbers::

    sage: print(":-(")    # abs tol 0.1
    :-)

Incorrect amount of numbers::

    sage: print("1.0 2.0 3.0")  # abs tol 0.1
    4.0 5.0
    sage: print("Hello")  # abs tol 0.1
    1.0
    sage: print("1.0")  # abs tol 0.1
    Hello

Good number but ordinary doctest failure::

    sage: print("Hello 1.1")  # abs tol 0.1
    Goodbye 1.0

Wrong number and ordinary doctest failure (both errors are reported)::

    sage: print("Hello 1.0")  # rel tol 1e-6
    Goodbye 0.999999

Hiding numbers under ellipsis (...) is not supported::

    sage: print("Hello 1.0")  # rel tol 1e-6
    Hello ...
