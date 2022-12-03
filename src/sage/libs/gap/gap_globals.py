# sage_setup: distribution = sagemath-gap
"""Common globals defined by GAP."""

###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


from .gap_functions import common_gap_functions


# selected gap globals to use in tab completion
common_gap_globals = set([
  'Assert',
  'Cyclotomics',
  'GaussianIntegers',
  'GaussianRationals',
  'GlobalMersenneTwister',
  'GlobalRandomSource',
  'InfoAlgebra',
  'InfoAttributes',
  'InfoBckt',
  'InfoCharacterTable',
  'InfoCoh',
  'InfoComplement',
  'InfoCoset',
  'InfoFpGroup',
  'InfoGroebner',
  'InfoGroup',
  'InfoLattice',
  'InfoMatrix',
  'InfoMonomial',
  'InfoNumtheor',
  'InfoOptions',
  'InfoPackageLoading',
  'InfoPcSubgroup',
  'InfoWarning',
  'Integers',
  'NiceBasisFiltersInfo',
  'Primes',
  'Rationals',
  'TableOfMarksComponents'
]) | common_gap_functions
