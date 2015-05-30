"""
Database of Modular Polynomials
"""

#######################################################################
#
#  Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#######################################################################


import bz2, os
from sage.env import SAGE_SHARE
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

dblocation = os.path.join(SAGE_SHARE, 'kohel')

def _dbz_to_integer_list(name):
    file = os.path.join(dblocation, name)
    if not os.path.exists(file):
        raise RuntimeError("Modular polynomial database file %s not available"%file)
    data = bz2.decompress(open(file).read())
    data = "[[" + data.replace("\n","],[").replace(" ",",")[:-2] + "]"
    return eval(data)

def _pad_int(s,n):
    return "0"*(n-len(str(s))) + str(s)

class ModularPolynomialDatabase:
    def _dbpath(self,level):
        return "PolMod/%s/pol.%s.dbz"%(self.model, _pad_int(level,3))

    def __repr__(self):
        if self.model == "Cls":
            head = "Classical"
            poly = "polynomial"
        elif self.model == "Atk":
            head = "Atkin"
            poly = "polynomial"
        elif self.model == "AtkCrr":
            head = "Atkin"
            poly = "correspondence"
        elif self.model == "Eta":
            head = "Dedekind eta"
            poly = "polynomial"
        elif self.model == "EtaCrr":
            head = "Dedekind eta"
            poly = "correspondence"
        return "%s modular %s database"%(head,poly)

    def __getitem__(self,level):
        """
        Return the modular polynomial of given level, or an error if
        there is no such polynomial in the database.

        EXAMPLES:
            sage: DBMP = ClassicalModularPolynomialDatabase()  # optional - database_kohel
            sage: f = DBMP[29]                                 #optional
            sage: f.degree()                                   #optional
            58
            sage: f.coefficient([28,28])                       #optional
            400152899204646997840260839128

            sage: DBMP[50]                                     #optional
            Traceback (most recent call last):
            ...
            RuntimeError: No database entry for modular polynomial of level 50
        """
        if self.model in ("Atk","Eta"):
            level = Integer(level)
            if not level.is_prime():
                raise TypeError("Argument level (= %s) must be prime."%level)
        elif self.model in ("AtkCrr","EtaCrr"):
            N = Integer(level[0])
            if not N in (2,3,5,7,13):
                raise TypeError("Argument level (= %s) must be prime."%N)
        modpol = self._dbpath(level)
        try:
            coeff_list = _dbz_to_integer_list(modpol)
        except RuntimeError as msg:
            print(msg)
            raise RuntimeError("No database entry for modular polynomial of level %s"%level)
        if self.model == "Cls":
            P = PolynomialRing(IntegerRing(),2,"j")
        else:
            P = PolynomialRing(IntegerRing(),2,"x,j")
        poly = {}
        if self.model == "Cls":
            if level == 1:
                return P({(1,0):1,(0,1):-1})
            for cff in coeff_list:
                i = cff[0]
                j = cff[1]
                poly[(i,j)] = Integer(cff[2])
                if i != j:
                    poly[(j,i)] = Integer(cff[2])
        else:
            for cff in coeff_list:
                poly[(cff[0],cff[1])] = Integer(cff[2])
        return P(poly)

class ModularCorrespondenceDatabase(ModularPolynomialDatabase):
    def _dbpath(self,level):
        (Nlevel,crrlevel) = level
        return "PolMod/%s/crr.%s.%s.dbz"%(
            self.model, _pad_int(Nlevel,2), _pad_int(crrlevel,3))

class ClassicalModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of classical modular polynomials, i.e. the polynomials
    Phi_N(X,Y) relating the j-functions j(q) and j(q^N).
    """
    def __init__(self):
        """
        Initialize the database.
        """
        self.model = "Cls"

class DedekindEtaModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of modular polynomials Phi_N(X,Y) relating a quotient
    of Dedekind eta functions, well-defined on X_0(N), relating x(q) and
    the j-function j(q).
    """
    def __init__(self):
        """
        Initialize the database.
        """
        self.model = "Eta"

class DedekindEtaModularCorrespondenceDatabase(ModularCorrespondenceDatabase):
    """
    The database of modular correspondences in $X_0(p) \times X_0(p)$, where
    the model of the curves $X_0(p) = \Bold{P}^1$ are specified by quotients of
    Dedekind's eta function.
    """
    def __init__(self):
        """
        Returns the
        """
        self.model = "EtaCrr"

class AtkinModularPolynomialDatabase(ModularPolynomialDatabase):
    """
    The database of modular polynomials Phi(x,j) for $X_0(p)$, where
    x is a function on invariant under the Atkin-Lehner invariant,
    with pole of minimal order at infinity.
    """
    def __init__(self):
        """
        Initialize the database.
        """
        self.model = "Atk"

class AtkinModularCorrespondenceDatabase(ModularCorrespondenceDatabase):
    def __init__(self):
        """
        Initialize the database.
        """
        self.model = "AtkCrr"
