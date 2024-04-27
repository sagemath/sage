#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL), v2+.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.stats.hmm.hmm", ["DiscreteHiddenMarkovModel"])
lazy_import("sage.stats.hmm.chmm", ["GaussianHiddenMarkovModel","GaussianMixtureHiddenMarkovModel"])
lazy_import("sage.stats.hmm.distributions", ["GaussianMixtureDistribution"])
