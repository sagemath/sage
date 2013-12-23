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
sys.path.append(os.environ['SAGE_DOC'])
from common.conf import *

# settings for the intersphinx extension:

ref_src = os.path.join(SAGE_DOC, 'en', 'reference')
ref_out = os.path.join(SAGE_DOC, 'output', 'html', 'en', 'reference')
intersphinx_mapping[ref_out] = None

for doc in os.listdir(ref_src):
    if os.path.exists(os.path.join(ref_src, doc, 'index.rst')):
        intersphinx_mapping[os.path.join(ref_out, doc)] = None

# We use the main document's title, if we can find it.
rst_file = open('index.rst', 'r')
rst_lines = rst_file.read().splitlines()
rst_file.close()

title = u''
for i in xrange(len(rst_lines)):
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
exclude_trees = exclude_trees + ['_sage']

multidocs_is_master = False
