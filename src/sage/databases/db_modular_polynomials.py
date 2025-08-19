"""
Database of modular polynomials

This module gives access to the database of modular polynomials. To use the
database, you need to install the optional :ref:`database_kohel
<spkg_database_kohel>` package by the Sage command ::

    sage -i database_kohel

EXAMPLES::

    sage: # optional - database_kohel
    sage: DBMP = ClassicalModularPolynomialDatabase()
    sage: f = DBMP[29]
    sage: f.degree()
    58
    sage: f.coefficient([28,28])
    400152899204646997840260839128

AUTHORS:

- David Kohel (2006-08-04): initial version
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
from pathlib import Path


def _dbz_to_string(name) -> str:
    r"""
    TESTS::

        sage: # optional - database_kohel
        sage: from sage.databases.db_modular_polynomials import _dbz_to_string
        sage: _dbz_to_string('PolMod/Atk/pol.002.dbz')
        '3 0 1 \n2 1 -1 \n2 0 744 \n1 1 -1 \n1 0 184512 \n0 2 1 \n0 1 7256 \n0 0 15252992 \n'
        sage: _dbz_to_string('PolMod/Cls/pol.001.dbz')
        '1 0 1 \n'
        sage: _dbz_to_string('PolMod/Eta/pol.002.dbz')
        '3 0 1 \n2 0 48 \n1 1 -1 \n1 0 768 \n0 0 4096 \n'
        sage: _dbz_to_string('PolMod/EtaCrr/crr.02.002.dbz')
        '2 1 1 \n2 0 -48 \n1 1 2304 \n0 2 -4096 \n0 1 196608 \n'
        sage: _dbz_to_string('PolHeeg/Cls/0000001-0005000/pol.0000003.dbz')
        '0\n1\n'
    """
    from sage.env import SAGE_SHARE
    filename = Path(SAGE_SHARE) / 'kohel' / name
    try:
        with open(filename, 'rb') as f:
            data = bz2.decompress(f.read())
    except OSError:
        raise FileNotFoundError('file not found in the Kohel database')
    return data.decode()


def _dbz_to_integer_list(name) -> list[list]:
    r"""
    TESTS::

        sage: # optional - database_kohel
        sage: from sage.databases.db_modular_polynomials import _dbz_to_integer_list
        sage: _dbz_to_integer_list('PolMod/Atk/pol.002.dbz')
        [[3, 0, 1],
         [2, 1, -1],
         [2, 0, 744],
         [1, 1, -1],
         [1, 0, 184512],
         [0, 2, 1],
         [0, 1, 7256],
         [0, 0, 15252992]]
        sage: _dbz_to_integer_list('PolMod/Cls/pol.001.dbz')
        [[1, 0, 1]]
        sage: _dbz_to_integer_list('PolMod/Eta/pol.002.dbz')
        [[3, 0, 1], [2, 0, 48], [1, 1, -1], [1, 0, 768], [0, 0, 4096]]
    """
    from sage.rings.integer import Integer
    data = _dbz_to_string(name)
    return [[Integer(v) for v in row.strip().split(" ")]
            for row in data.split("\n")[:-1]]


def _dbz_to_integers(name) -> list:
    r"""
    TESTS::

        sage: from sage.databases.db_modular_polynomials import _dbz_to_integers
        sage: _dbz_to_integers('PolHeeg/Cls/0000001-0005000/pol.0000003.dbz') # optional - database_kohel
        [0, 1]
    """
    from sage.rings.integer import Integer
    return [Integer(i) for i in _dbz_to_string(name).split()]


class ModularPolynomialDatabase:
    def _dbpath(self, level) -> Path:
        r"""
        TESTS::

            sage: C = ClassicalModularPolynomialDatabase()
            sage: C._dbpath(3)
            PosixPath('PolMod/Cls/pol.003.dbz')
            sage: C._dbpath(8)
            PosixPath('PolMod/Cls/pol.008.dbz')
        """
        path = Path("PolMod")
        return path / self.model / ("pol.%03d.dbz" % level)

    def __repr__(self) -> str:
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

            sage: # optional - database_kohel
            sage: DBMP = ClassicalModularPolynomialDatabase()
            sage: f = DBMP[29]
            sage: f.degree()
            58
            sage: f.coefficient([28,28])
            400152899204646997840260839128
            sage: DBMP[50]
            Traceback (most recent call last):
            ...
            FileNotFoundError: file not found in the Kohel database
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
    def _dbpath(self, level) -> Path:
        r"""
        TESTS::

            sage: DB = DedekindEtaModularCorrespondenceDatabase()
            sage: DB._dbpath((2,4))
            PosixPath('PolMod/EtaCrr/crr.02.004.dbz')
        """
        path = Path("PolMod")
        return path / self.model / ("crr.%02d.%03d.dbz" % level)


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
    The database of modular correspondences in `X_0(p) \times X_0(p)`, where
    the model of the curves `X_0(p) = \Bold{P}^1` are specified by quotients of
    Dedekind's eta function.
    """
    model = "EtaCrr"


class AtkinModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of modular polynomials Phi(x,j) for `X_0(p)`, where
    x is a function on invariant under the Atkin-Lehner invariant,
    with pole of minimal order at infinity.
    """
    model = "Atk"


class AtkinModularCorrespondenceDatabase(ModularCorrespondenceDatabase):
    model = "AtkCrr"
