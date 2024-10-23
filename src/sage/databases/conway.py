r"""
Frank Lübeck's tables of Conway polynomials over finite fields
"""
# ****************************************************************************
#
#       Copyright (C) 2005-2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2010      Alexandru Ghitza
#       Copyright (C) 2013      R. Andrew Ohana <andrew.ohana@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from collections.abc import Mapping


class DictInMapping(Mapping):
    def __init__(self, dict):
        """
        Places dict into a non-mutable mapping.

        TESTS::

            sage: from sage.databases.conway import DictInMapping
            sage: d = {}
            sage: m = DictInMapping(d); m
            {}
            sage: d[0] = 1; m
            {0: 1}
            sage: m[2] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'DictInMapping' object does not support item assignment
        """
        self._store = dict

    def __getitem__(self, key):
        """
        TESTS::

            sage: from sage.databases.conway import DictInMapping
            sage: DictInMapping({'foo': 'bar'})['foo']
            'bar'
        """
        return self._store[key]

    def __len__(self):
        """
        TESTS::

            sage: from sage.databases.conway import DictInMapping
            sage: d = {}
            sage: m = DictInMapping(d); len(m)
            0
            sage: d['foo'] = 'bar'; len(m)
            1
        """
        return len(self._store)

    def __iter__(self):
        """
        TESTS::

            sage: from sage.databases.conway import DictInMapping
            sage: next(iter(DictInMapping({'foo': 'bar'})))
            'foo'
        """
        return iter(self._store)

    def __repr__(self):
        """
        TESTS::

            sage: from sage.databases.conway import DictInMapping
            sage: DictInMapping({'foo': 'bar'})
            {'foo': 'bar'}
        """
        return repr(self._store)


class ConwayPolynomials(Mapping):
    def __init__(self):
        """
        Initialize the database.

        TESTS::

            sage: c = ConwayPolynomials()
            sage: c
            Frank Lübeck's database of Conway polynomials
        """
        import conway_polynomials
        self._store = conway_polynomials.database()

    def __repr__(self):
        """
        Return a description of this database.

        TESTS::

            sage: c = ConwayPolynomials()
            sage: c.__repr__()
            "Frank Lübeck's database of Conway polynomials"
        """
        return "Frank Lübeck's database of Conway polynomials"

    def __getitem__(self, key):
        """
        If key is a pair of integers ``p,n``, return the Conway
        polynomial of degree ``n`` over ``GF(p)``.

        If key is an integer ``p``, return a non-mutable mapping
        whose keys are the degrees of the polynomial values.

        TESTS::

            sage: c = ConwayPolynomials()
            sage: c[60859]
            {1: (60856, 1), 2: (3, 60854, 1),
                    3: (60856, 8, 0, 1), 4: (3, 32881, 3, 0, 1)}
            sage: c[60869, 3]
            (60867, 2, 0, 1)
        """
        try:
            return DictInMapping(self._store[key])
        except KeyError as err:
            try:
                if isinstance(key, (tuple, list)):
                    if len(key) == 2:
                        return self._store[key[0]][key[1]]
            except KeyError:
                pass
            raise err

    def __len__(self):
        """
        Return the number of polynomials in this database.

        TESTS:

        The database currently contains `35357` polynomials, but due to
        :issue:`35357` it will be extended by Conway polynomials of
        degrees `1`, `2` and `3` for primes between `65537` and `110000`,
        thus leading to a new total of `47090` entries::

            sage: c = ConwayPolynomials()
            sage: len(c) in [35357, 47090]
            True
        """
        try:
            return self._len
        except AttributeError:
            pass
        self._len = sum(len(a) for a in self._store.values())
        return self._len

    def __iter__(self):
        """
        Return an iterator over the keys of this database.

        TESTS::

            sage: c = ConwayPolynomials()
            sage: itr = iter(c)
            sage: next(itr)  # random
            (65537, 4)
        """
        for a, b in self._store.items():
            for c in b:
                yield a, c

    def polynomial(self, p, n):
        """
        Return the Conway polynomial of degree ``n`` over ``GF(p)``,
        or raise a :exc:`RuntimeError` if this polynomial is not in the
        database.

        .. NOTE::

            See also the global function ``conway_polynomial`` for
            a more user-friendly way of accessing the polynomial.

        INPUT:

        - ``p`` -- prime number

        - ``n`` -- positive integer

        OUTPUT:

        List of Python int's giving the coefficients of the corresponding
        Conway polynomial in ascending order of degree.

        EXAMPLES::

            sage: c = ConwayPolynomials()
            sage: c.polynomial(3, 21)
            (1, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
            sage: c.polynomial(97, 128)
            Traceback (most recent call last):
            ...
            RuntimeError: Conway polynomial over F_97 of degree 128 not in database.
        """
        try:
            return self[p, n]
        except KeyError:
            raise RuntimeError("Conway polynomial over F_%s of degree %s not in database." % (p, n))

    def has_polynomial(self, p, n):
        """
        Return ``True`` if the database of Conway polynomials contains the
        polynomial of degree ``n`` over ``GF(p)``.

        INPUT:

        - ``p`` -- prime number

        - ``n`` -- positive integer

        EXAMPLES::

            sage: c = ConwayPolynomials()
            sage: c.has_polynomial(97, 12)
            True
            sage: c.has_polynomial(60821, 5)
            False
        """
        return (p,n) in self

    def primes(self):
        """
        Return the list of prime numbers ``p`` for which the database of
        Conway polynomials contains polynomials over ``GF(p)``.

        EXAMPLES::

            sage: c = ConwayPolynomials()
            sage: P = c.primes()
            sage: 2 in P
            True
            sage: next_prime(10^7) in P
            False
        """
        return self._store.keys()

    def degrees(self, p):
        """
        Return the list of integers ``n`` for which the database of Conway
        polynomials contains the polynomial of degree ``n`` over ``GF(p)``.

        EXAMPLES::

            sage: c = ConwayPolynomials()
            sage: c.degrees(60821)
            [1, 2, 3, 4]
            sage: c.degrees(next_prime(10^7))
            []
        """
        if p not in self._store:
            return []
        return list(self._store[p])

    def __reduce__(self):
        """
        TESTS::

            sage: c = ConwayPolynomials()
            sage: loads(dumps(c)) == c
            True
        """
        return (ConwayPolynomials, ())
