# sage_setup: distribution = sagemath-modules
from sage.misc.lazy_import import lazy_import
lazy_import('sage.crypto.mq.rijndael_gf', 'RijndaelGF')
lazy_import('sage.crypto.mq.sr', 'SR')
