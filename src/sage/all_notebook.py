"""
All imports for the Sage notebook
"""

import sys
from sage.all import *
from sagenb.notebook.all import *
from sage.calculus.predefined import x
from sage.misc.python import python
from sage.misc.html import html, pretty_print_default, MathJax
from sagenb.misc.support import help, automatic_names

sage_mode = 'notebook'

from sage.misc.latex import Latex
latex = Latex(density=130)
latex_debug = Latex(debug=True, density=130)
slide = Latex(slide=True, density=256)
slide_debug = Latex(slide=True, debug=True, density=256)
pdflatex = Latex(density=130, pdflatex=True)
pdflatex_debug = Latex(density=130, pdflatex=True, debug=True)

# Set user globals to notebook worksheet globals
from sage.repl.user_globals import set_globals
set_globals(sys._getframe(1).f_globals)

sage.misc.session.init()
