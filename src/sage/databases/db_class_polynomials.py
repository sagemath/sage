"""
Database of Hilbert Polynomials
"""
# ****************************************************************************
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2016 Vincent Delecroix <vincent.delecroix@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .db_modular_polynomials import _dbz_to_integers

disc_format = "%07d"  #  disc_length = 7
level_format = "%03d" #  level_length = 3


class ClassPolynomialDatabase:
    def _dbpath(self, disc, level=1):
        """
        TESTS::

            sage: db = HilbertClassPolynomialDatabase()
            sage: db._dbpath(5, 1)
            'PolHeeg/Cls/0000001-0005000/pol.0000005.dbz'
            sage: db._dbpath(15000, 1)
            'PolHeeg/Cls/0010001-0015000/pol.0015000.dbz'
            sage: db._dbpath(3, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Level (= 2) > 1 not yet implemented
        """
        if level != 1:
            raise NotImplementedError("Level (= %s) > 1 not yet implemented"%level)
        n1 = 5000*((abs(disc)-1)//5000)
        s1 = disc_format % (n1+1) #_pad_int(n1+1, disc_length)
        s2 = disc_format % (n1+5000)
        subdir = "%s-%s" % (s1, s2)
        discstr = disc_format % abs(disc)
        return "PolHeeg/%s/%s/pol.%s.dbz" % (self.model, subdir, discstr)

    def __getitem__(self, disc):
        r"""
        TESTS::

            sage: db = HilbertClassPolynomialDatabase()
            sage: db[32]  # optional - database_kohel
            x^2 - 52250000*x + 12167000000
            sage: db[123913912]
            Traceback (most recent call last):
            ...
            ValueError: file not found in the Kohel database
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        classpol = self._dbpath(disc)
        coeff_list = _dbz_to_integers(classpol)
        return PolynomialRing(ZZ, 'x')(coeff_list)


class HilbertClassPolynomialDatabase(ClassPolynomialDatabase):
    """
    The database of Hilbert class polynomials.

    EXAMPLES::

        sage: db = HilbertClassPolynomialDatabase()
        sage: db[-4]                     # optional - database_kohel
        x - 1728
        sage: db[-7]                     # optional - database_kohel
        x + 3375
        sage: f = db[-23]; f             # optional - database_kohel
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
        sage: f.discriminant().factor()  # optional - database_kohel
        -1 * 5^18 * 7^12 * 11^4 * 17^2 * 19^2 * 23
        sage: db[-23]                    # optional - database_kohel
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
    """
    model = "Cls"

    def __repr__(self):
        return "Hilbert class polynomial database"

######################################################
# None of the following are implemented yet.
######################################################

class AtkinClassPolynomialDatabase(ClassPolynomialDatabase):
    """
    The database of Atkin class polynomials.
    """
    model = "Atk"

    def __repr__(self):
        return "Atkin class polynomial database"


class WeberClassPolynomialDatabase(ClassPolynomialDatabase):
    """
    The database of Weber class polynomials.
    """
    def __repr__(self):
        return "Weber class polynomial database"


class DedekindEtaClassPolynomialDatabase(ClassPolynomialDatabase):
    """
    The database of Dedekind eta class polynomials.
    """
    model = "Eta"

    def __repr__(self):
        return "Dedekind eta class polynomial database"
