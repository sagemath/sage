from __future__ import absolute_import
from .sr import SR
from .sbox import SBox

from sage.misc.lazy_import import lazy_import
lazy_import('sage.crypto.mq.rijndael_gf', 'RijndaelGF')
