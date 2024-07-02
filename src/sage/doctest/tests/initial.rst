One initial typo causes a :exc:`NameError` in the first test and many
following failures::

    sage: a = binomiak(10,5)  # random to test that we still get the exception
    sage: a
    252
    sage: a == factorial(10)/factorial(5)^2
    True
    sage: a//12
    21

But this is a new failure::

    sage: binomial(10,5)
    255
