# -*- coding: utf-8 -*-
#
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

import sys, os
from sage.env import SAGE_DOC_SRC, SAGE_DOC
sys.path.append(SAGE_DOC_SRC)
from common.conf import *

from six.moves import range

ref_src = os.path.join(SAGE_DOC_SRC, 'en', 'reference')
ref_out = os.path.join(SAGE_DOC, 'html', 'en', 'reference')

# We use the main document's title, if we can find it.
rst_file = open('index.rst', 'r')
rst_lines = rst_file.read().splitlines()
rst_file.close()

title = u''
for i in range(len(rst_lines)):
    if rst_lines[i].startswith('==') and i > 0:
        title = rst_lines[i-1].strip()
        break

# Otherwise, we use this directory's name.
name = os.path.basename(os.path.abspath('.'))
if not title:
    title = name.capitalize()
title = title.replace(u'`', u'$')

# General information about the project.
project = u'Sage Reference Manual: ' + title

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = u'Sage Reference Manual v' + release + ': ' + title

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = title

# HTML theme (e.g., 'default', 'sphinxdoc').  The pages for the
# reference manual use a custom theme, a slight variant on the 'sage'
# theme, to set the links in the top line.
html_theme = 'sageref'

# Output file base name for HTML help builder.
htmlhelp_basename = name

# Grouping the document tree into LaTeX files. List of tuples (source
# start file, target name, title, author, document class
# [howto/manual]).
latex_documents = [
('index', name + '.tex', project, u'The Sage Development Team', 'manual')
]

#Ignore all .rst in the _sage subdirectory
exclude_patterns = exclude_patterns + ['_sage']

multidocs_is_master = False
