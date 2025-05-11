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
from sage_docbuild.conf import release, exclude_patterns
from sage_docbuild.conf import *

# Add any paths that contain custom static files (such as style sheets),
# relative to this directory to html_static_path. They are copied after the
# builtin static files, so a file named "default.css" will overwrite the
# builtin "default.css". html_common_static_path imported from sage_docbuild.conf
# contains common paths.
html_static_path = [] + html_common_static_path

ref_src = os.path.join(SAGE_DOC_SRC, 'en', 'reference')
ref_out = os.path.join(SAGE_DOC, 'html', 'en', 'reference')

# We use the main document's title, if we can find it.
rst_file = open('index.rst', 'r')
rst_lines = rst_file.read().splitlines()
rst_file.close()

title = ''
for i in range(len(rst_lines)):
    if rst_lines[i].startswith('==') and i > 0:
        title = rst_lines[i-1].strip()
        break

# Otherwise, we use this directory's name.
name = os.path.basename(os.path.abspath('.'))
if not title:
    title = name.capitalize()
title = title.replace('`', '$')

# We use the directory's name to add small view/edit buttons.
html_theme_options.update({
  'source_view_link': os.path.join(source_repository, f'blob/develop/src/doc/en/reference/{name}', '{filename}'),
  'source_edit_link': os.path.join(source_repository, f'edit/develop/src/doc/en/reference/{name}', '{filename}'),
})

# General information about the project.
project = title

# The name for this set of Sphinx documents.
html_title = title
html_short_title = title

# Output file base name for HTML help builder.
htmlhelp_basename = name

# Grouping the document tree into LaTeX files. List of tuples (source
# start file, target name, title, author, document class
# [howto/manual]).
latex_documents = [
  ('index', name + '.tex', title,
   'The Sage Development Team', 'manual')
]

latex_elements['hyperref'] = r"""
\usepackage{xr}
\externaldocument[../references/]{../references/references}
% Include hyperref last.
\usepackage{hyperref}
% Fix anchor placement for figures with captions.
\usepackage{hypcap}% it must be loaded after hyperref.
% Set up styles of URL: it should be placed after hyperref.
\urlstyle{same}"""

#Ignore all .rst in the _sage subdirectory
exclude_patterns = exclude_patterns + ['_sage']

multidocs_is_master = False
