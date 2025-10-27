"Exceptions raised by the coercion model"

###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################


class CoercionException(TypeError):
    """
    This is the baseclass of exceptions that the coercion model raises
    when trying to discover coercions. We do not use standard Python
    exceptions to avoid inadvertently catching and suppressing real errors.

    Usually one raises this to indicate the attempted action is not
    implemented/appropriate, but if there are other things to try not
    to immediately abort to the user.
    """
    pass
