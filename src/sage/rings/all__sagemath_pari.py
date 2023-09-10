# sage_setup: distribution = sagemath-pari

# Pseudo-ring of PARI objects.
from .pari_ring import PariRing, Pari

# p-adic field
from .padics.all import *
from .padics.padic_printing import _printer_defaults as padic_printing
