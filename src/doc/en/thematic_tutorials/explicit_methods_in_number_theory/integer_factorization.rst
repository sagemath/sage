Integer Factorization
=====================


Quadratic Sieve
---------------

Bill Hart's quadratic sieve is included with Sage. The quadratic sieve
is the best algorithm for factoring numbers of the form :math:`pq` up
to around 100 digits. It involves searching for relations, solving a
linear algebra problem modulo :math:`2`, then factoring :math:`n`
using a relation :math:`x^2 \equiv y^2 \mod n`.

::

    sage: qsieve(next_prime(2^90)*next_prime(2^91), time=True)   # not tested
    ([1237940039285380274899124357, 2475880078570760549798248507],
     '14.94user 0.53system 0:15.72elapsed 98%CPU (0avgtext+0avgdata 0maxresident)k')

Using qsieve is twice as fast as Sage's general factor command in
this example. Note that Sage's general factor command does nothing
but call Pari's factor C library function.

::

    sage: time factor(next_prime(2^90)*next_prime(2^91))     # not tested
    CPU times: user 28.71 s, sys: 0.28 s, total: 28.98 s
    Wall time: 29.38 s
    1237940039285380274899124357 * 2475880078570760549798248507

Obviously, Sage's factor command should not just call Pari, but
nobody has gotten around to rewriting it yet.

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

In the following example, GMP-ECM is over 10 times faster than
Sage's generic factor function. Again, this emphasizes that Sage's
generic factor command would benefit from a rewrite that uses
GMP-ECM and qsieve.

::

    sage: time ecm.factor(next_prime(2^40) * next_prime(2^300))    # not tested
    CPU times: user 0.85 s, sys: 0.01 s, total: 0.86 s
    Wall time: 1.73 s
    [1099511627791,
     2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533]
    sage: time factor(next_prime(2^40) * next_prime(2^300))        # not tested
    CPU times: user 23.82 s, sys: 0.04 s, total: 23.86 s
    Wall time: 24.35 s
    1099511627791 * 2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533
