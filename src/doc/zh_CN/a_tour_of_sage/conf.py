# nodoctest
# Sage documentation build configuration file, created by
# sphinx-quickstart on Thu Aug 21 20:15:55 2008.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default; values that are commented out
# serve to show the default.

from sage_docbuild.conf import release
from sage_docbuild.conf import *  # NOQA

# Add any paths that contain custom static files (such as style sheets),
# relative to this directory to html_static_path. They are copied after the
# builtin static files, so a file named "default.css" will overwrite the
# builtin "default.css". html_common_static_path imported from sage_docbuild.conf
# contains common paths.
html_static_path = [] + html_common_static_path

# General information about the project.
project = "Sage Tutorial in Chinese (Simplified)"
name = 'tutorial-zh_CN'
language = "zh_CN"

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = project
html_short_title = project

# Output file base name for HTML help builder.
htmlhelp_basename = 'a_tour_of_sage'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [
  ('index', 'SageTutorial_zh_CN.tex', project,
   'The Sage Development Team', 'manual'),
]

latex_elements['babel'] = ''
latex_elements['inputenc'] = ''
latex_elements['fontenc'] = ''
latex_elements['utf8extra'] = ''

# https://tex.stackexchange.com/questions/223893/how-do-i-find-out-what-chinese-fonts-are-installed-with-my-mactex-installation
latex_elements['preamble'] = r"""
\usepackage[UTF8]{ctex}
\setCJKmainfont[BoldFont=FandolSong-Bold.otf]{FandolSong-Regular.otf}
\setCJKsansfont[BoldFont=FandolHei-Bold.otf]{FandolHei-Regular.otf}
\setCJKmonofont{FandolFang-Regular.otf}
\newCJKfontfamily\kaiti{FandolKai-Regular.otf}
\usepackage{amsmath}
\usepackage{amssymb}
"""
