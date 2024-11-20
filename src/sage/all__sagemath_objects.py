# sage_setup: distribution = sagemath-objects
import os
import sys
import operator
import math

import warnings

# TODO: More to be moved from all.py

# This import also sets up the interrupt handler
from cysignals.signals import (AlarmInterrupt, SignalError,
                               sig_on_reset as sig_on_count)

from time import sleep

from sage.misc.all__sagemath_objects import *
from sage.structure.all import *
from sage.arith.power import generic_power as power
from sage.categories.all__sagemath_objects import *

from sage.cpython.all import *

from copy import copy, deepcopy

true = True
false = False


# For doctesting. These are overwritten later

Integer = int
RealNumber = float
