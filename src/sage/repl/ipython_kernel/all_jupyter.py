# sage_setup: distribution = sagemath-repl
"""
All imports for Jupyter
"""

from sage.all_cmdline import *

from sage.repl.ipython_kernel.widgets_sagenb import (input_box, text_control, slider,
                                                     range_slider, checkbox, selector, input_grid, color_selector)
from sage.repl.ipython_kernel.interact import interact
