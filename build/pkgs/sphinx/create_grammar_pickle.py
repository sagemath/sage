# Code taken from sphinx/pycode/__init__.py to generate the grammar pickle.

from os import path
from sphinx import package_dir
from sphinx.pycode.pgen2 import driver

_grammarfile = path.join(package_dir, 'pycode', 'Grammar-py2.txt')
pygrammar = driver.load_grammar(_grammarfile)
