"""
Database of Modular Polynomials
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2016 Vincent Delecroix <vincent.delecroix@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import bz2
import os
from sage.cpython.string import bytes_to_str


def _dbz_to_string(name):
    r"""
    TESTS::

        sage: from sage.databases.db_modular_polynomials import _dbz_to_string
        sage: _dbz_to_string('PolMod/Atk/pol.002.dbz')                      # optional - database_kohel
        '3 0 1 \n2 1 -1 \n2 0 744 \n1 1 -1 \n1 0 184512 \n0 2 1 \n0 1 7256 \n0 0 15252992 \n'
        sage: _dbz_to_string('PolMod/Cls/pol.001.dbz')                      # optional - database_kohel
        '1 0 1 \n'
        sage: _dbz_to_string('PolMod/Eta/pol.002.dbz')                      # optional - database_kohel
        '3 0 1 \n2 0 48 \n1 1 -1 \n1 0 768 \n0 0 4096 \n'
        sage: _dbz_to_string('PolMod/EtaCrr/crr.02.002.dbz')                # optional - database_kohel
        '2 1 1 \n2 0 -48 \n1 1 2304 \n0 2 -4096 \n0 1 196608 \n'
        sage: _dbz_to_string('PolHeeg/Cls/0000001-0005000/pol.0000003.dbz') # optional - database_kohel
        '0\n1\n'
    """
    from sage.env import SAGE_SHARE
    dblocation = os.path.join(SAGE_SHARE, 'kohel')
    filename = os.path.join(dblocation, name)
    try:
        with open(filename, 'rb') as f:
            data = bz2.decompress(f.read())
    except IOError:
        raise ValueError('file not found in the Kohel database')
    return bytes_to_str(data)


def _dbz_to_integer_list(name):
    r"""
    TESTS::

        sage: from sage.databases.db_modular_polynomials import _dbz_to_integer_list
        sage: _dbz_to_integer_list('PolMod/Atk/pol.002.dbz') # optional - database_kohel
        [[3, 0, 1],
         [2, 1, -1],
         [2, 0, 744],
         [1, 1, -1],
         [1, 0, 184512],
         [0, 2, 1],
         [0, 1, 7256],
         [0, 0, 15252992]]
        sage: _dbz_to_integer_list('PolMod/Cls/pol.001.dbz') # optional - database_kohel
        [[1, 0, 1]]
        sage: _dbz_to_integer_list('PolMod/Eta/pol.002.dbz') # optional - database_kohel
        [[3, 0, 1], [2, 0, 48], [1, 1, -1], [1, 0, 768], [0, 0, 4096]]
    """
    from sage.rings.integer import Integer
    data = _dbz_to_string(name)
    return [[Integer(v) for v in row.strip().split(" ")]
            for row in data.split("\n")[:-1]]


def _dbz_to_integers(name):
    r"""
    TESTS::

        sage: from sage.databases.db_modular_polynomials import _dbz_to_integers
        sage: _dbz_to_integers('PolHeeg/Cls/0000001-0005000/pol.0000003.dbz') # optional - database_kohel
        [0, 1]
    """
    from sage.rings.integer import Integer
    return [Integer(i) for i in _dbz_to_string(name).split()]


class ModularPolynomialDatabase:
    def _dbpath(self, level):
        r"""
        TESTS::

            sage: C = ClassicalModularPolynomialDatabase()
            sage: C._dbpath(3)
            'PolMod/Cls/pol.003.dbz'
            sage: C._dbpath(8)
            'PolMod/Cls/pol.008.dbz'
        """
        return "PolMod/%s/pol.%03d.dbz" % (self.model, level)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: ClassicalModularPolynomialDatabase()
            Classical modular polynomial database

            sage: DedekindEtaModularPolynomialDatabase()
            Dedekind eta modular polynomial database
            sage: DedekindEtaModularPolynomialDatabase()
            Dedekind eta modular polynomial database

            sage: AtkinModularPolynomialDatabase()
            Atkin modular polynomial database
        """
        if self.model.startswith("Cls"):
            head = "Classical"
        elif self.model.startswith("Atk"):
            head = "Atkin"
        elif self.model.startswith("Eta"):
            head = "Dedekind eta"

        if self.model.endswith("Crr"):
            poly = "correspondence"
        else:
            poly = "polynomial"

        return "%s modular %s database" % (head, poly)

    def __getitem__(self, level):
        """
        Return the modular polynomial of given level, or an error if
        there is no such polynomial in the database.

        EXAMPLES::

            sage: DBMP = ClassicalModularPolynomialDatabase()
            sage: f = DBMP[29]                                 # optional - database_kohel
            sage: f.degree()                                   # optional - database_kohel
            58
            sage: f.coefficient([28,28])                       # optional - database_kohel
            400152899204646997840260839128

            sage: DBMP[50]                                     # optional - database_kohel
            Traceback (most recent call last):
            ...
            ValueError: file not found in the Kohel database
        """
        from sage.rings.integer import Integer
        from sage.rings.integer_ring import IntegerRing
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if self.model in ("Atk", "Eta"):
            level = Integer(level)
            if not level.is_prime():
                raise TypeError("Argument level (= %s) must be prime." % level)
        elif self.model in ("AtkCrr", "EtaCrr"):
            N = Integer(level[0])
            if N not in (2, 3, 5, 7, 13):
                raise TypeError("Argument level (= %s) must be prime." % N)
        modpol = self._dbpath(level)
        coeff_list = _dbz_to_integer_list(modpol)
        if self.model == "Cls":
            P = PolynomialRing(IntegerRing(), 2, "j")
        else:
            P = PolynomialRing(IntegerRing(), 2, "x,j")
        poly = {}
        if self.model == "Cls":
            if level == 1:
                return P({(1, 0): 1, (0, 1): -1})
            for cff in coeff_list:
                i = cff[0]
                j = cff[1]
                poly[(i, j)] = Integer(cff[2])
                if i != j:
                    poly[(j, i)] = Integer(cff[2])
        else:
            for cff in coeff_list:
                poly[(cff[0], cff[1])] = Integer(cff[2])
        return P(poly)


class ModularCorrespondenceDatabase(ModularPolynomialDatabase):
    def _dbpath(self, level):
        r"""
        TESTS::

            sage: DB = DedekindEtaModularCorrespondenceDatabase()
            sage: DB._dbpath((2,4))
            'PolMod/EtaCrr/crr.02.004.dbz'
        """
        (Nlevel, crrlevel) = level
        return "PolMod/%s/crr.%02d.%03d.dbz" % (self.model, Nlevel, crrlevel)


class ClassicalModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of classical modular polynomials, i.e. the polynomials
    Phi_N(X,Y) relating the j-functions j(q) and j(q^N).
    """
    model = "Cls"


class DedekindEtaModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of modular polynomials Phi_N(X,Y) relating a quotient
    of Dedekind eta functions, well-defined on X_0(N), relating x(q) and
    the j-function j(q).
    """
    model = "Eta"


class DedekindEtaModularCorrespondenceDatabase(ModularCorrespondenceDatabase):
    r"""
    The database of modular correspondences in $X_0(p) \times X_0(p)$, where
    the model of the curves $X_0(p) = \Bold{P}^1$ are specified by quotients of
    Dedekind's eta function.
    """
    model = "EtaCrr"


class AtkinModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of modular polynomials Phi(x,j) for $X_0(p)$, where
    x is a function on invariant under the Atkin-Lehner invariant,
    with pole of minimal order at infinity.
    """
    model = "Atk"


class AtkinModularCorrespondenceDatabase(ModularCorrespondenceDatabase):
    model = "AtkCrr"
