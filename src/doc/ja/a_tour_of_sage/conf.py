# nodoctest
# Sage documentation build configuration file, based on that created by
# sphinx-quickstart on Thu Aug 21 20:15:55 2008.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default; values that are commented out
# serve to show the default.

from sage_docbuild.conf import release, latex_elements
from sage_docbuild.conf import *  # NOQA


for tag in feature_tags():
    tags.add(tag)


# Add any paths that contain custom static files (such as style sheets),
# relative to this directory to html_static_path. They are copied after the
# builtin static files, so a file named "default.css" will overwrite the
# builtin "default.css". html_common_static_path imported from sage_docbuild.conf
# contains common paths.
html_static_path = [] + html_common_static_path

# Add small view/edit buttons.
html_theme_options.update({
  'source_view_link': os.path.join(source_repository, 'blob/develop/src/doc/ja/a_tour_of_sage', '{filename}'),
  'source_edit_link': os.path.join(source_repository, 'edit/develop/src/doc/ja/a_tour_of_sage', '{filename}'),
})

# General information about the project.
project = "Sage ガイドツアー"
name = 'a_tour_of_sage'
language = "ja"

# The LaTeX engine to build the docs in Japanese.
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-latex_engine
latex_engine = 'uplatex'

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = project + " v" + release

# Output file base name for HTML help builder.
htmlhelp_basename = name

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [
  ('index', name + '.tex', project,
   'The Sage Group', 'manual'),
]

# PDF output: let long decimal expansions in code-blocks wrap rather than
# overflow beyond page margin
latex_elements = {
    'sphinxsetup': 'verbatimforcewraps=true',
}
