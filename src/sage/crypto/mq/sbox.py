from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.sbox', ['SBox',
                                 'feistel_construction',
                                 'misty_construction'],
            deprecation=22986)
