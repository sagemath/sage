"""
Selection of the members to document

This extension decides which members Sphinx includes in the reference
manual, and optionally reports the nested classes whose pickling is broken.
"""
# ****************************************************************************
#       Copyright (C) 2026 Chenxin Zhong <chenxin.zhong@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import sys

from sphinx.util import logging as sphinx_logging

logger = sphinx_logging.getLogger(__name__)


skip_picklability_check_modules = [
    #'sage.misc.test_nested_class', # for test only
    'sage.misc.latex',
    'sage.misc.explain_pickle',
    '__builtin__',
]


def check_nested_class_picklability(app, what, name, obj, skip, options):
    """
    Print a warning if pickling is broken for nested classes.
    """
    if hasattr(obj, '__dict__') and hasattr(obj, '__module__'):
        # Check picklability of nested classes.  Adapted from
        # sage.misc.nested_class.modify_for_nested_pickle.
        module = sys.modules[obj.__module__]
        for (nm, v) in obj.__dict__.items():
            if (isinstance(v, type) and
                v.__name__ == nm and
                v.__module__ == module.__name__ and
                getattr(module, nm, None) is not v and
                v.__module__ not in skip_picklability_check_modules):
                # OK, probably this is an *unpicklable* nested class.
                logger.warning('Pickling of nested class %r is probably broken. '
                               'Please set the metaclass of the parent class to '
                               'sage.misc.nested_class.NestedClassMetaclass.',
                               v.__module__ + '.' + name + '.' + nm)


def skip_member(app, what, name, obj, skip, options):
    """
    To suppress Sphinx warnings / errors, we

    - Don't include [aliases of] builtins.

    - Don't include the docstring for any nested class which has been
      inserted into its module by
      :class:`sage.misc.nested_class.NestedClassMetaclass` only for pickling.  The
      class will be properly documented inside its surrounding class.

    - Optionally, check whether pickling is broken for nested classes.

    - Optionally, include objects whose name begins with an underscore
      ('_'), i.e., "private" or "hidden" attributes, methods, etc.

    Otherwise, we abide by Sphinx's decision.  Note: The object
    ``obj`` is excluded (included) if this handler returns True
    (False).
    """
    if 'SAGE_CHECK_NESTED' in os.environ:
        check_nested_class_picklability(app, what, name, obj, skip, options)

    if getattr(obj, '__module__', None) == '__builtin__':
        return True

    objname = getattr(obj, "__name__", None)
    if objname is not None:
        # check if name was inserted to the module by NestedClassMetaclass
        if name.find('.') != -1 and objname.find('.') != -1:
            if objname.split('.')[-1] == name.split('.')[-1]:
                return True

    if 'SAGE_DOC_UNDERSCORE' in os.environ:
        if name.split('.')[-1].startswith('_'):
            return False

    return skip


def setup(app):
    """
    Register this extension with Sphinx.
    """
    app.connect('autodoc-skip-member', skip_member)
    return {'parallel_read_safe': True}
