"""
Database of Hilbert Polynomials
"""

#######################################################################

#  Sage: System for Algebra and Geometry Experimentation
#
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#######################################################################


import bz2, os
from sage.env import SAGE_SHARE
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

dblocation = os.path.join(SAGE_SHARE, 'kohel')

disc_length = 7
level_length = 3

def _dbz_to_integers(name):
    file = '%s/%s'%(dblocation, name)
    if not os.path.exists(file):
        raise RuntimeError("Class polynomial database file %s not available"%file)
    data = bz2.decompress(open(file).read())
    data = "[" + data.replace("\n",",") + "]"
    return eval(data)

def _pad_int_str(s,n):
    return "0"*(n-len(str(s))) + str(s)

class ClassPolynomialDatabase:
    def _dbpath(self,disc,level=1):
        """
        TESTS:
            sage: db = HilbertClassPolynomialDatabase()
            sage: db.__getitem__(-23,level=2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Level (= 2) > 1 not yet implemented.
        """
        if level != 1:
            raise NotImplementedError("Level (= %s) > 1 not yet implemented."%level)
        n1 = 5000*((abs(disc)-1)//5000)
        s1 = _pad_int_str(n1+1,disc_length)
        s2 = _pad_int_str(n1+5000,disc_length)
        subdir = "%s-%s"%(s1,s2)
        discstr = _pad_int_str(abs(disc),disc_length)
        return "PolHeeg/%s/%s/pol.%s.dbz"%(self.model,subdir,discstr)

    def __getitem__(self,disc,level=1,var='x'):
        classpol = self._dbpath(disc,level)
        try:
            coeff_list = _dbz_to_integers(classpol)
        except RuntimeError as msg:
            print(msg)
            raise RuntimeError("No database entry for class polynomial of discriminant %s"%disc)
        P = PolynomialRing(IntegerRing(),names=var)
        return P(list(coeff_list))

class HilbertClassPolynomialDatabase(ClassPolynomialDatabase):
    """
    The database of Hilbert class polynomials.

    EXAMPLES::

        sage: db = HilbertClassPolynomialDatabase()
        sage: db[-4]                 # optional - database_kohel
        x - 1728
        sage: db[-7]                 # optional
        x + 3375
        sage: f = db[-23]; f         # optional
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
        sage: f.discriminant().factor()    # optional
        -1 * 5^18 * 7^12 * 11^4 * 17^2 * 19^2 * 23
        sage: db[-23]                      # optional
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
    """
    def __init__(self):
        """
        Initialize the database.
        """
        self.model = "Cls"

    def __repr__(self):
        return "Hilbert class polynomial database"

######################################################
# None of the following are implemented yet.
######################################################

class AtkinClassPolynomialDatabase(ClassPolynomialDatabase):
    """
    The database of Atkin class polynomials.
    """
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
    def __repr__(self):
        return "Dedekind eta class polynomial database"

