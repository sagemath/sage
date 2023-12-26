from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.block_cipher.miniaes', 'MiniAES')
lazy_import('sage.crypto.block_cipher.sdes', 'SimplifiedDES')
lazy_import('sage.crypto.block_cipher.present', 'PRESENT')
del lazy_import
