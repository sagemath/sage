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

import os
from sage.env import SAGE_DOC_SRC, SAGE_DOC
from sage.docs.conf import release, latex_elements, exclude_patterns
from sage.docs.conf import *

ref_src = os.path.join(SAGE_DOC_SRC, 'en', 'reference')
ref_out = os.path.join(SAGE_DOC, 'html', 'en', 'reference')

# General information about the project.
project = u"Sage Reference Manual"
name = "reference"

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = project + " v"+release

html_short_title = u'Sage Reference v' + release

# Output file base name for HTML help builder.
htmlhelp_basename = name

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [
  ('index', name + '.tex', u'Sage Reference Manual',
   u'The Sage Development Team', 'manual'),
]

latex_elements['preamble'] += r'''
% One-column index
\makeatletter
\renewenvironment{theindex}{
  \chapter*{\indexname}
  \markboth{\MakeUppercase\indexname}{\MakeUppercase\indexname}
  \setlength{\parskip}{0.1em}
  \relax
  \let\item\@idxitem
}{}
\makeatother
'''

#Ignore all .rst in the _sage subdirectory
exclude_patterns = exclude_patterns + ['_sage']

multidocs_is_master = True

# Sorted list of subdocs. Include all subdirectories of ref_src except
# for 'static' and 'templates', and to deal with upgrades: 'sage',
# 'media', and 'other'.
bad_directories = ['static', 'templates', 'sage', 'media', 'other']
multidocs_subdoc_list = sorted([x for x in os.listdir(ref_src)
                                if os.path.isdir(os.path.join(ref_src, x))
                                and x not in bad_directories])

# List of directories, relative to source directory, that shouldn't be
# searched for source files.
exclude_patterns += multidocs_subdoc_list + [
    'sage', 'options'
    ]
