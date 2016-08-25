"""
This file gathers together all the tables in Sage.

    * ConwayPolynomials() -- database of selected Conway polynomials.

    * CremonaDatabase() - Cremona's tables of elliptic curves and related data.

    * findstat -- The FindStat database (http://www.findstat.org/).

    * JonesDatabase() -- returns the John Jones table of number fields
      with bounded ramification and degree <= 6.

    * oeis -- The On-Line Encyclopedia of Integer Sequences (http://oeis.org/).

    * SloaneEncyclopedia -- Local copy of Sloane On-Line Encyclopedia of
      Integer Sequences.

    * SteinWatkinsAllData() and SteinWatkinsPrimeData() - The
      Stein-Watkins tables of elliptic curves and related data.

    * SymbolicData() -- many benchmark and testing ideals

EXAMPLES::

    sage: ConwayPolynomials()
    Frank Luebeck's database of Conway polynomials

    sage: CremonaDatabase()
    Cremona's database of elliptic curves with conductor...

    sage: JonesDatabase()
    John Jones's table of number fields with bounded ramification and degree <= 6

    sage: oeis
    The On-Line Encyclopedia of Integer Sequences (http://oeis.org/)

    sage: SymbolicData()
    SymbolicData with ... ideals
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from .sql_db import SQLQuery, SQLDatabase

from .conway import ConwayPolynomials

from .cremona import CremonaDatabase

from .jones import JonesDatabase

from .stein_watkins import SteinWatkinsAllData, SteinWatkinsPrimeData

from .sloane import sloane_sequence, sloane_find, SloaneEncyclopedia

from sage.misc.lazy_import import lazy_import
lazy_import('sage.databases.oeis', 'oeis')

from .symbolic_data import SymbolicData

lazy_import('sage.databases.odlyzko', 'zeta_zeros')

from .db_modular_polynomials import \
     ClassicalModularPolynomialDatabase, \
     DedekindEtaModularPolynomialDatabase, \
     DedekindEtaModularCorrespondenceDatabase, \
     AtkinModularPolynomialDatabase, \
     AtkinModularCorrespondenceDatabase

from .db_class_polynomials import \
     HilbertClassPolynomialDatabase

from .cunningham_tables import cunningham_prime_factors

lazy_import('sage.databases.findstat', 'findstat')
