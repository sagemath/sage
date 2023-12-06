# sage_setup: distribution = sagemath-pari

from sage.rings.padics.all__sagemath_categories import *

from sage.rings.padics.factory import Zp, Zq, Zp as pAdicRing, ZpCR, ZpCA, ZpFM, ZpFP, ZpLC, ZpLF, ZqCR, ZqCA, ZqFM, ZqFP, ZpER
from sage.rings.padics.factory import Qp, Qq, Qp as pAdicField, QpCR, QpFP, QpLC, QpLF, QqCR, QqFP, QpER
from sage.rings.padics.factory import pAdicExtension

from sage.rings.padics.padic_printing import _printer_defaults as padic_printing

from sage.rings.padics.pow_computer import PowComputer
