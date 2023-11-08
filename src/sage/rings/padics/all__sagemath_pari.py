
from .all__sagemath_categories import *

from .factory import Zp, Zq, Zp as pAdicRing, ZpCR, ZpCA, ZpFM, ZpFP, ZpLC, ZpLF, ZqCR, ZqCA, ZqFM, ZqFP, ZpER
from .factory import Qp, Qq, Qp as pAdicField, QpCR, QpFP, QpLC, QpLF, QqCR, QqFP, QpER
from .factory import pAdicExtension

from .padic_printing import _printer_defaults as padic_printing

from .pow_computer import PowComputer
