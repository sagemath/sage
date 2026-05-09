Integer Factorization
=====================


Quadratic Sieve
---------------

Bill Hart's quadratic sieve is included with Sage. The quadratic sieve
is one of the best algorithms for factoring numbers of the form
:math:`pq` up to around 100 digits. It involves searching for
relations, solving a linear algebra problem modulo :math:`2`, then
factoring :math:`n` using a relation :math:`x^2 \equiv y^2 \mod n`.
Using the qsieve algorithm can be faster than the default, which
uses PARI.

::

    sage: n = next_prime(2^90)*next_prime(2^91)
    sage: n.factor(algorithm="qsieve")
    doctest:... RuntimeWarning: the factorization returned
    by qsieve may be incomplete (the factors may not be prime)
    or even wrong; see qsieve? for details
    1237940039285380274899124357 * 2475880078570760549798248507
    sage: n.factor()  # uses PARI at the time of writing
    1237940039285380274899124357 * 2475880078570760549798248507


GMP-ECM
-------

Paul Zimmerman's GMP-ECM is included in Sage. The elliptic curve
factorization (ECM) algorithm is the best algorithm for factoring
numbers of the form :math:`n=pm`, where :math:`p` is not "too
big". ECM is an algorithm due to Hendrik Lenstra, which works by
"pretending" that :math:`n` is prime, choosing a random elliptic curve
over :math:`\ZZ/n\ZZ`, and doing arithmetic on that
curve--if something goes wrong when doing arithmetic, we factor
:math:`n`.

In the following example, GMP-ECM is much faster than Sage's generic
factor function. Again, this emphasizes that the best factoring
algorithm may depend on your specific problem.

::

    sage: n = next_prime(2^40) * next_prime(2^300)
    sage: n.factor(algorithm="ecm")
    1099511627791 * 2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533
    sage: n.factor()  # uses PARI at the time of writing
    1099511627791 * 2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533
