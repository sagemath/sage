r"""
Index of channels

The ``channels`` object may be used to access the codes that Sage can build.

- :func:`channel_constructions.ErrorErasureChannel <sage.coding.channel_constructions.ErrorErasureChannel>`
- :func:`channel_constructions.StaticErrorRateChannel <sage.coding.channel_constructions.StaticErrorRateChannel>`

- :func:`channel_constructions.QarySymmetricChannel <sage.coding.channel_constructions.QarySymmetricChannel>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.channels_catalog import *

"""
from __future__ import absolute_import

from .channel_constructions import (ErrorErasureChannel, StaticErrorRateChannel, QarySymmetricChannel)
