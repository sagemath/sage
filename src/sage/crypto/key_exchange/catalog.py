"""
Index of key exchange schemes

This catalogue includes implementations of key exchange schemes.

Let ``<tab>`` indicate pressing the :kbd:`Tab` key. So begin by typing
``key_exchange.<tab>`` to the see the currently implemented key exchange
schemes.

This catalogue includes the following key exchange schemes:

- :class:`sage.crypto.key_exchange.diffie_hellman.DiffieHellman`

To import these names into the global namespace, use::

    sage: from sage.crypto.key_exchange.catalog import *
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.key_exchange.diffie_hellman', 'DiffieHellman')

del lazy_import
