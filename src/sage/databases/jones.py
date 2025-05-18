r"""
John Jones's tables of number fields

In order to use the Jones database, the optional :ref:`database_jones_numfield
<spkg_database_jones_numfield>` package must be installed using the Sage command ::

    sage -i database_jones_numfield

This is a table of number fields with bounded ramification and
degree `\leq 6`. You can query the database for all number
fields in Jones's tables with bounded ramification and degree.

EXAMPLES: First load the database::

    sage: J = JonesDatabase()
    sage: J
    John Jones's table of number fields with bounded ramification and degree <= 6

List the degree and discriminant of all fields in the database that
have ramification at most at 2::

    sage: [(k.degree(), k.disc()) for k in J.unramified_outside([2])]    # optional - database_jones_numfield
    [(1, 1), (2, -4), (2, -8), (2, 8), (4, 256), (4, 512), (4, -1024), (4, -2048), (4, 2048), (4, 2048), (4, 2048)]

List the discriminants of the fields of degree exactly 2 unramified
outside 2::

    sage: [k.disc() for k in J.unramified_outside([2],2)]                # optional - database_jones_numfield
    [-4, -8, 8]

List the discriminants of cubic field in the database ramified
exactly at 3 and 5::

    sage: [k.disc() for k in J.ramified_at([3,5],3)]                     # optional - database_jones_numfield
    [-135, -675, -6075, -6075]
    sage: factor(6075)
    3^5 * 5^2
    sage: factor(675)
    3^3 * 5^2
    sage: factor(135)
    3^3 * 5

List all fields in the database ramified at 101::

    sage: J.ramified_at(101)                                             # optional - database_jones_numfield
    [Number Field in a with defining polynomial x^2 - 101,
     Number Field in a with defining polynomial x^4 - x^3 + 13*x^2 - 19*x + 361,
     Number Field in a with defining polynomial x^5 + x^4 - 6*x^3 - x^2 + 18*x + 4,
     Number Field in a with defining polynomial x^5 + 2*x^4 + 7*x^3 + 4*x^2 + 11*x - 6,
     Number Field in a with defining polynomial x^5 - x^4 - 40*x^3 - 93*x^2 - 21*x + 17]
"""

# ****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os

from sage.rings.number_field.number_field import NumberField
from sage.rings.rational_field import RationalField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.subset import powerset
from sage.env import SAGE_SHARE

from sage.misc.persist import load, save

from sage.features.databases import DatabaseJones


def sortkey(K):
    """
    A completely deterministic sorting key for number fields.

    EXAMPLES::

        sage: from sage.databases.jones import sortkey
        sage: sortkey(QuadraticField(-3))
        (2, 3, False, x^2 + 3)
    """
    return K.degree(), abs(K.discriminant()), K.discriminant() > 0, K.polynomial()


class JonesDatabase:
    def __init__(self):
        self.root = None

    def __repr__(self):
        return "John Jones's table of number fields with bounded ramification and degree <= 6"

    def _load(self, path, filename):
        print(filename)
        i = 0
        while filename[i].isalpha():
            i += 1
        j = len(filename) - 1
        while filename[j].isalpha() or filename[j] in [".", "_"]:
            j -= 1
        S = sorted([eval(z) for z in filename[i:j + 1].split("-")])
        with open(path + "/" + filename) as f:
            data = f.read()
        data = data.replace("^", "**")
        x = PolynomialRing(RationalField(), 'x').gen()  # used next line
        v = eval(data)
        s = tuple(S)
        if s in self.root:
            self.root[s] += v
            self.root[s].sort()
        else:
            self.root[s] = v

    def _init(self, path):
        """
        Create the database from scratch from the PARI files on John Jones's
        web page, downloaded (e.g., via wget) to a local directory, which
        is specified as path above.

        INPUT:

        - ``path`` -- (default: works on William Stein install)
          path must be the path to Jones's Number_Fields directory
          (http://hobbes.la.asu.edu/Number_Fields); these files should have
          been downloaded using wget

        EXAMPLES: This is how to create the database from scratch, assuming
        that the number fields are in the default directory above: From a
        cold start of Sage::

                sage: J = JonesDatabase()
                sage: J._init()   # not tested
                ...

        This takes about 5 seconds.
        """
        x = PolynomialRing(RationalField(), 'x').gen()
        self.root = {}
        self.root[tuple()] = [x - 1]
        if not os.path.exists(path):
            raise OSError("Path %s does not exist." % path)
        for X in os.listdir(path):
            if X[-4:] == "solo":
                Z = path + "/" + X
                print(X)
                for Y in os.listdir(Z):
                    if Y[-3:] == ".gp":
                        self._load(Z, Y)

        data_dir = os.path.dirname(DatabaseJones().absolute_filename())
        os.makedirs(data_dir, exist_ok=True)
        save(self.root, os.path.join(data_dir, "jones.sobj"))

    def unramified_outside(self, S, d=None, var='a'):
        """
        Return all fields in the database of degree d unramified
        outside S. If d is omitted, return fields of any degree up to 6.
        The fields are ordered by degree and discriminant.

        INPUT:

        - ``S`` -- list or set of primes, or a single prime

        - ``d`` -- ``None`` (default, in which case all fields of degree <= 6 are returned)
          or a positive integer giving the degree of the number fields returned

        - ``var`` -- the name used for the generator of the number fields (default: ``'a'``)

        EXAMPLES::

            sage: J = JonesDatabase()             # optional - database_jones_numfield
            sage: J.unramified_outside([101,109]) # optional - database_jones_numfield
            [Number Field in a with defining polynomial x - 1,
             Number Field in a with defining polynomial x^2 - 101,
             Number Field in a with defining polynomial x^2 - 109,
             Number Field in a with defining polynomial x^3 - x^2 - 36*x + 4,
             Number Field in a with defining polynomial x^4 - x^3 + 13*x^2 - 19*x + 361,
             Number Field in a with defining polynomial x^4 - x^3 + 14*x^2 + 34*x + 393,
             Number Field in a with defining polynomial x^5 + x^4 - 6*x^3 - x^2 + 18*x + 4,
             Number Field in a with defining polynomial x^5 + 2*x^4 + 7*x^3 + 4*x^2 + 11*x - 6,
             Number Field in a with defining polynomial x^5 - x^4 - 40*x^3 - 93*x^2 - 21*x + 17]
        """
        try:
            S = list(S)
        except TypeError:
            S = [S]
        Z = []
        for X in powerset(S):
            Z += self.ramified_at(X, d=d, var=var)
        return sorted(Z, key=sortkey)

    def __getitem__(self, S):
        return self.get(S)

    def get(self, S, var='a'):
        """
        Return all fields in the database ramified exactly at
        the primes in S.

        INPUT:

        - ``S`` -- list or set of primes, or a single prime

        - ``var`` -- the name used for the generator of the number fields (default: ``'a'``)

        EXAMPLES::

            sage: J = JonesDatabase()              # optional - database_jones_numfield
            sage: J.get(163, var='z')              # optional - database_jones_numfield
            [Number Field in z with defining polynomial x^2 + 163,
             Number Field in z with defining polynomial x^3 - x^2 - 54*x + 169,
             Number Field in z with defining polynomial x^4 - x^3 - 7*x^2 + 2*x + 9]
            sage: J.get([3, 4])                    # optional - database_jones_numfield
            Traceback (most recent call last):
            ...
            ValueError: S must be a list of primes
        """
        if self.root is None:
            self.root = load(DatabaseJones().absolute_filename())
        try:
            S = list(S)
        except TypeError:
            S = [S]
        if not all(p.is_prime() for p in S):
            raise ValueError("S must be a list of primes")
        S.sort()
        s = tuple(S)
        if s not in self.root:
            return []
        return [NumberField(f, var, check=False) for f in self.root[s]]

    def ramified_at(self, S, d=None, var='a'):
        """
        Return all fields in the database of degree d ramified exactly at
        the primes in S.  The fields are ordered by degree and discriminant.

        INPUT:

        - ``S`` -- list or set of primes

        - ``d`` -- ``None`` (default, in which case all fields of degree <= 6 are returned)
          or a positive integer giving the degree of the number fields returned

        - ``var`` -- the name used for the generator of the number fields (default: ``'a'``)

        EXAMPLES::

            sage: # optional - database_jones_numfield
            sage: J = JonesDatabase()
            sage: J.ramified_at([101,109])
            []
            sage: J.ramified_at([109])
            [Number Field in a with defining polynomial x^2 - 109,
             Number Field in a with defining polynomial x^3 - x^2 - 36*x + 4,
             Number Field in a with defining polynomial x^4 - x^3 + 14*x^2 + 34*x + 393]
            sage: J.ramified_at(101)
            [Number Field in a with defining polynomial x^2 - 101,
             Number Field in a with defining polynomial x^4 - x^3 + 13*x^2 - 19*x + 361,
             Number Field in a with defining polynomial x^5 + x^4 - 6*x^3 - x^2 + 18*x + 4,
             Number Field in a with defining polynomial x^5 + 2*x^4 + 7*x^3 + 4*x^2 + 11*x - 6,
             Number Field in a with defining polynomial x^5 - x^4 - 40*x^3 - 93*x^2 - 21*x + 17]
            sage: J.ramified_at((2, 5, 29), 3, 'c')
            [Number Field in c with defining polynomial x^3 - x^2 - 8*x - 28,
             Number Field in c with defining polynomial x^3 - x^2 + 10*x + 102,
             Number Field in c with defining polynomial x^3 - x^2 - 48*x - 188,
             Number Field in c with defining polynomial x^3 - x^2 + 97*x - 333]
        """
        Z = self.get(S, var=var)
        if d is not None:
            Z = [k for k in Z if k.degree() == d]
        return sorted(Z, key=sortkey)
