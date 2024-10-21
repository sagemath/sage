# sage_setup: distribution = sagemath-modules
import sage.crypto.sbox
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.crypto.mq.sbox', 'SBox', sage.crypto.sbox.SBox)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.classical', ['AffineCryptosystem',
                                      'HillCryptosystem',
                                      'SubstitutionCryptosystem',
                                      'ShiftCryptosystem',
                                      'TranspositionCryptosystem',
                                      'VigenereCryptosystem',
                                      ])

lazy_import('sage.crypto.stream', ['LFSRCryptosystem',
                                   'ShrinkingGeneratorCryptosystem',
                                   ])

lazy_import('sage.crypto.lfsr', ['lfsr_sequence',
                                 'lfsr_autocorrelation',
                                 'lfsr_connection_polynomial',
                                 ])


import sage.crypto.key_exchange.catalog as key_exchange

del lazy_import
