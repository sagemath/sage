
from .free_monoid import FreeMonoid
from .string_monoid import (BinaryStrings, OctalStrings, HexadecimalStrings,
                            Radix64Strings, AlphabeticStrings)

from .free_abelian_monoid import FreeAbelianMonoid

from .string_ops import (
    strip_encoding,
    frequency_distribution,
    coincidence_index,
    coincidence_discriminant)
