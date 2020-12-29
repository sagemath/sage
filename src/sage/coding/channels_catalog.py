r"""
Index of channels

Channels in Sage implement the information theoretic notion of transmission of messages.

The ``channels`` object may be used to access the codes that Sage can build.

- :class:`channel.ErrorErasureChannel <sage.coding.channel.ErrorErasureChannel>`
- :class:`channel.QarySymmetricChannel <sage.coding.channel.QarySymmetricChannel>`
- :class:`channel.StaticErrorRateChannel <sage.coding.channel.StaticErrorRateChannel>`

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

from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import('sage.coding.channel', ['ErrorErasureChannel',
                                     'QarySymmetricChannel',
                                     'StaticErrorRateChannel'])
