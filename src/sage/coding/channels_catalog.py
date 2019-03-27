r"""
Index of Channels: the information theoretic notion of transmission

The ``channels`` object may be used to access the codes that Sage can build.

- :class:`channel_constructions.ErrorErasureChannel <sage.coding.channel_constructions.ErrorErasureChannel>`
- :class:`channel_constructions.QarySymmetricChannel <sage.coding.channel_constructions.QarySymmetricChannel>`
- :class:`channel_constructions.StaticErrorRateChannel <sage.coding.channel_constructions.StaticErrorRateChannel>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.channels_catalog import *

"""
#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import('sage.coding.channel_constructions', ['ErrorErasureChannel',
                                                   'QarySymmetricChannel',
                                                   'StaticErrorRateChannel'])
# We don't want this to appear in tab completion
del absolute_import
