sphinx: Python documentation generator
======================================

Description
-----------

Sphinx is a tool that makes it easy to create intelligent and beautiful
documentation for Python projects (or other documents consisting of
multiple reStructuredText sources), written by Georg Brandl. It was
originally created to translate the new Python documentation, but has
now been cleaned up in the hope that it will be useful to many other
projects.

License
-------

Modified BSD; see e.g. its egg-info file for other options


Upstream Contact
----------------

- Author: Georg Brandl
- Home Page: http://www.sphinx-doc.org
- see also http://pypi.python.org/pypi/Sphinx

Dependencies
------------

-  six >= 1.4
-  Jinja2 >= 2.3
-  Pygments >= 2.0
-  docutils >= 0.11
-  snowballstemmer >= 1.1
-  babel >= 1.3
-  setuptools / distribute
-  Python
-  GNU patch (shipped with Sage)


Special Update/Build Instructions
---------------------------------

-  The script create_grammar_pickle.py creates the file
   Grammar2.7.pickle in site-packages/Sphinx-.../sphinx/pycode/. This
   helps to avoid race conditions when building the documentation in
   parallel.
