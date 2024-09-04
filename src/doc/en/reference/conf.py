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

import os
from sage.env import SAGE_DOC_SRC, SAGE_DOC
from sage_docbuild.conf import release, latex_elements, exclude_patterns
from sage_docbuild.conf import *


for tag in feature_tags():
    tags.add(tag)


# Add any paths that contain custom static files (such as style sheets),
# relative to this directory to html_static_path. They are copied after the
# builtin static files, so a file named "default.css" will overwrite the
# builtin "default.css". html_common_static_path imported from sage_docbuild.conf
# contains common paths.
html_static_path = [] + html_common_static_path

ref_src = os.path.join(SAGE_DOC_SRC, 'en', 'reference')
ref_out = os.path.join(SAGE_DOC, 'html', 'en', 'reference')

# Add small view/edit buttons.
html_theme_options.update({
  'source_view_link': os.path.join(source_repository, 'blob/develop/src/doc/en/reference/index.rst'),
  'source_edit_link': os.path.join(source_repository, 'edit/develop/src/doc/en/reference/index.rst'),
})

# General information about the project.
project = "Reference Manual"

# The name for this set of Sphinx documents.
html_title = project
html_short_title = project

# Output file base name for HTML help builder.
htmlhelp_basename = 'reference'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [
  ('index', 'reference.tex', 'Reference Manual',
   'The Sage Development Team', 'manual'),
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
