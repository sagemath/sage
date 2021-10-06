# TODO: More to be moved from all.py

import os
import sys
import operator
import math

# This import also sets up the interrupt handler
from cysignals.signals import (AlarmInterrupt, SignalError,
        sig_on_reset as sig_on_count)

from time                import sleep
from functools import reduce  # in order to keep reduce in python3

from sage.misc.all__sagemath_objects       import *
from sage.structure.all  import *
from sage.categories.all__sagemath_objects import *


from cysignals.alarm import alarm, cancel_alarm

from copy import copy, deepcopy
