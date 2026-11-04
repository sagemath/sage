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

lazy_import('sage.crypto.ed25519', ['Ed25519',
                                    'Ed25519BasePoint',
                                    'ED25519_A',
                                    'ED25519_BASE_X',
                                    'ED25519_BASE_Y',
                                    'ED25519_COFACTOR',
                                    'ED25519_D',
                                    'ED25519_FIELD_SIZE',
                                    'ED25519_SQRT_M1',
                                    'ED25519_SUBGROUP_ORDER',
                                    'ed25519_decode',
                                    'ed25519_encode',
                                    'ed25519_public_key',
                                    'ed25519_sign',
                                    'ed25519_verify',
                                    ])

lazy_import('sage.crypto.lfsr', ['lfsr_sequence',
                                 'lfsr_autocorrelation',
                                 'lfsr_connection_polynomial',
                                 ])

lazy_import('sage.crypto.public_key.key_exchange', 'all', 'key_exchange')

del lazy_import
