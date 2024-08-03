from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.key_exchange.diffie_hellman', 'DiffieHellman')
lazy_import('sage.crypto.key_exchange.key_exchange_scheme', 'KeyExchangeScheme')

del lazy_import
