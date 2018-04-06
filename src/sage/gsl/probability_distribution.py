"""
Probability Distributions

This module provides three types of probability distributions:

- ``RealDistribution``: various real-valued probability distributions.

- ``SphericalDistribution``: uniformly distributed points on the
  surface of an `n-1` sphere in `n` dimensional euclidean space.

- ``GeneralDiscreteDistribution``: user-defined discrete distributions.

AUTHORS:

- Josh Kantor (2007-02): first version

- William Stein (2007-02): rewrite of docs, conventions, etc.

- Carlo Hamalainen (2008-08): full doctest coverage, more documentation,
  GeneralDiscreteDistribution, misc fixes.

- Kwankyu Lee (2010-05-29): F-distribution support.

REFERENCES:

    GNU gsl library, General discrete distributions
    http://www.gnu.org/software/gsl/manual/html_node/General-Discrete-Distributions.html

    GNU gsl library, Random number distributions
    http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html
"""

#*****************************************************************************
#       Copyright (C) 2004, 2005, 2006 Joshua Kantor <kantor.jm@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(9084, "the module sage.gsl.probability_distribution has moved"
            " to sage.probability.probability_distribution")

from sage.probability.probability_distribution import *
