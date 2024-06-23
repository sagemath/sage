# sage_setup: distribution = sagemath-objects
import os
import sys
import operator
import math
import warnings

from copy import copy, deepcopy
from functools import reduce  # in order to keep reduce in python3
from time import sleep

# TODO: More to be moved from all.py

# This import also sets up the interrupt handler
from cysignals.signals import (AlarmInterrupt, SignalError,
                               sig_on_reset as sig_on_count)

from cysignals.alarm import alarm, cancel_alarm

# Full and final
from sage.cpython.all import *
from sage.structure.all import *

# Partial
from sage.arith.all__sagemath_objects import *
from sage.categories.all__sagemath_objects import *
from sage.misc.all__sagemath_objects import *

true = True
false = False

# For doctesting. These are overwritten later

Integer = int
RealNumber = float
