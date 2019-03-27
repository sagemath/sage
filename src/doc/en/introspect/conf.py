# -*- coding: utf-8 -*-
# Sage introspection build configuration file.
# See sagenb.notebook.cell.Cell.set_introspect_html() for details.

import sys, os
from sage.env import SAGE_DOC_SRC
sys.path.append(SAGE_DOC_SRC)
from common.conf import *

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.mathjax', 'sphinx.ext.todo',
              'sphinx.ext.extlinks']

templates_path = ['templates']
html_static_path = ['static']

html_use_modindex = False
html_use_index = False
html_split_index = False
html_copy_source = False

todo_include_todos = True
