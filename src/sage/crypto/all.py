from __future__ import absolute_import

from .classical import (
    AffineCryptosystem,
    HillCryptosystem,
    SubstitutionCryptosystem,
    ShiftCryptosystem,
    TranspositionCryptosystem,
    VigenereCryptosystem)

from .stream import (
    LFSRCryptosystem,
    ShrinkingGeneratorCryptosystem)

from .lfsr import (
    lfsr_sequence,
    lfsr_autocorrelation,
    lfsr_connection_polynomial)
