"""
Cross-references between Sage and other projects

This extension configures intersphinx for the documentation of the projects
that Sage refers to, and resolves the references that Sphinx alone cannot,
either because the target is documented under a public name other than the
one that the reference names, or because it has no target to link to.
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

import builtins
import functools
import os
import sys
import typing

from sphinx.ext import intersphinx
from sphinx.util import logging as sphinx_logging
from sphinx.util.inspect import safe_getattr

from sage.env import PPLPY_DOCS, SAGE_DOC, SAGE_DOC_SRC

# Reading an attribute of a Sage module can run a lazy import, which raises if
# it names an optional package that is not installed.  An unresolved reference
# must not make the build fail over that, so every attribute of an object that
# is not ours is read with safe_getattr(), which answers with its default
# whatever goes wrong.

logger = sphinx_logging.getLogger(__name__)


SAGE_DOC_REMOTE_INVENTORIES = os.environ.get('SAGE_DOC_REMOTE_INVENTORIES', 'no') == 'yes'

_vendored_inventories_dir = os.path.join(SAGE_DOC_SRC, "common", "_vendor")


# Run "sage -python -m sage_docbuild.vendor" to update src/doc/common/_vendor/*.inv
_intersphinx_targets = {
    'cvxopt':     ['https://cvxopt.org/userguide/'],
    'cvxpy':      ['https://www.cvxpy.org/'],
    'cypari2':    ['https://cypari2.readthedocs.io/en/latest/'],
    'cysignals':  ['https://cysignals.readthedocs.io/en/latest/'],
    'flint':      ['https://flintlib.org/doc/'],
    'fpylll':     ['https://fpylll.readthedocs.io/en/latest/'],
    'gmpy2':      ['https://gmpy2.readthedocs.io/en/latest/'],
    'ipykernel':  ['https://ipykernel.readthedocs.io/en/stable/'],
    'ipython':    ['https://ipython.readthedocs.io/en/stable/'],
    'ipywidgets': ['https://ipywidgets.readthedocs.io/en/stable/'],
    'matplotlib': ['https://matplotlib.org/stable/'],
    'mpmath':     ['https://mpmath.org/doc/current/'],
    'networkx':   ['https://networkx.org/documentation/stable/'],
    'numpy':      ['https://numpy.org/doc/stable/'],
    'pexpect':    ['https://pexpect.readthedocs.io/en/stable/'],
    'pplpy':      [PPLPY_DOCS, 'https://www.sagemath.org/pplpy/'],
    'ptyprocess': ['https://ptyprocess.readthedocs.io/en/stable/'],
    'python':     ['https://docs.python.org/'],
    'rpy2':       ['https://rpy2.github.io/doc/latest/html/'],
    'scipy':      ['https://docs.scipy.org/doc/scipy/'],
    'sphinx':     ['https://www.sphinx-doc.org/en/master/'],
    'sympy':      ['https://docs.sympy.org/latest/'],
    'traitlets':  ['https://traitlets.readthedocs.io/en/stable/'],
}


def _intersphinx_mapping(key):
    inventories = []
    link_target = None
    for target in _intersphinx_targets[key]:
        if not target:
            pass
        elif target.startswith('http'):
            if not link_target:
                link_target = target
                if SAGE_DOC_REMOTE_INVENTORIES:
                    inventories.append(None)  # Try downloading inventory from link_target
        elif os.path.exists(target):
            if not link_target:
                link_target = target
            inventory = os.path.join(target, 'objects.inv')
            if os.path.exists(inventory):
                inventories.append(inventory)
                break
    else:
        vendored_inventory = os.path.join(_vendored_inventories_dir, key + '.inv')
        if os.path.exists(vendored_inventory):
            inventories.append(vendored_inventory)
        else:
            # To avoid docbuild failures when building Sage without internet
            # connection, we use the local python inventory file as a fallback for other
            # projects. Cross-references will not be resolved in that case, but the
            # docbuild will still succeed.
            python_inventory_file = os.path.join(_vendored_inventories_dir, "python.inv")
            inventories.append(python_inventory_file)
    assert link_target
    if len(inventories) == 1:
        return link_target, inventories[0]
    return link_target, tuple(inventories)


def set_intersphinx_mappings(app, config):
    """
    Add precompiled inventory (the objects.inv)
    """
    app.config.intersphinx_mapping = {}

    refpath = os.path.join(SAGE_DOC, "html", "en", "reference")
    invpath = os.path.join(SAGE_DOC, "inventory", "en", "reference")
    if app.config.multidoc_first_pass == 1 or not os.path.exists(invpath):
        return

    install_path = os.path.join(SAGE_DOC, "html", "en", "installation")
    install_inv = os.path.join(SAGE_DOC, "inventory", "en", "installation", "objects.inv")
    if os.path.exists(install_inv):
        app.config.intersphinx_mapping['installation'] = (install_path, install_inv)

    app.config.intersphinx_mapping.update({key: _intersphinx_mapping(key)
                                           for key in _intersphinx_targets})

    # Add master intersphinx mapping
    dst = os.path.join(invpath, 'objects.inv')
    if os.path.exists(dst):
        app.config.intersphinx_mapping['sagemath'] = (refpath, dst)

    # Add intersphinx mapping for subdirectories
    for directory in os.listdir(os.path.join(invpath)):
        if directory == 'jupyter_execute':
            # This directory is created by jupyter-sphinx extension for
            # internal use and should be ignored here. See Issue #33507.
            continue
        if os.path.isdir(os.path.join(invpath, directory)):
            src = os.path.join(refpath, directory)
            dst = os.path.join(invpath, directory, 'objects.inv')
            app.config.intersphinx_mapping[directory] = (src, dst)

    intersphinx.validate_intersphinx_mapping(app, config)


def prefer_python_inventory(app):
    r"""
    Give the Python inventory the last word on the builtins it shares.

    Intersphinx merges every inventory into a single one, project by
    project in the order of the project names, so that the project whose
    name sorts last silently claims each name that two projects define.
    Sphinx documents ``enumerate`` as a function of its own in the
    quickstart guide, for one, which is enough to divert every reference to
    the builtin away from the Python manual.

    Only the builtins are taken back, and only in the Python domain: a
    name that a project legitimately owns keeps its target, be it a
    glossary term, a document, a label or a function of the C API that
    the Python manual happens to name too.

    TESTS:

    Sphinx sorts after Python, hence wins every name the two share::

        sage: from types import SimpleNamespace
        sage: from sphinx.ext.intersphinx import InventoryAdapter
        sage: from sage_docbuild.ext.crossrefs import prefer_python_inventory
        sage: python = {'py:function': {'enumerate': 'python/enumerate'},
        ....:           'py:class': {'os.PathLike': 'python/pathlike'},
        ....:           'std:term': {'object': 'python/object'},
        ....:           'std:doc': {'tutorial/index': 'python/tutorial'},
        ....:           'std:label': {'glossary': 'python/glossary'},
        ....:           'c:function': {'PyType_GenericAlloc': 'python/alloc'}}
        sage: sphinx = {'py:function': {'enumerate': 'sphinx/enumerate'},
        ....:           'py:class': {'os.PathLike': 'sphinx/pathlike',
        ....:                        'sphinx.application.Sphinx': 'sphinx/app'},
        ....:           'std:term': {'object': 'sphinx/object'},
        ....:           'std:doc': {'tutorial/index': 'sphinx/tutorial'},
        ....:           'std:label': {'glossary': 'sphinx/glossary'},
        ....:           'c:function': {'PyType_GenericAlloc': 'sphinx/alloc'}}
        sage: app = SimpleNamespace(env=SimpleNamespace())
        sage: inventories = InventoryAdapter(app.env)
        sage: inventories.named_inventory.update(python=python, sphinx=sphinx)
        sage: for inventory in (python, sphinx):
        ....:     for objtype, objects in inventory.items():
        ....:         inventories.main_inventory.setdefault(objtype, {}).update(objects)
        sage: main = inventories.main_inventory
        sage: main['py:function']['enumerate']
        'sphinx/enumerate'
        sage: main['py:class']['os.PathLike']
        'sphinx/pathlike'

    A builtin of the Python domain goes back to the Python manual::

        sage: prefer_python_inventory(app)
        sage: main['py:function']['enumerate']
        'python/enumerate'

    Nothing else moves, neither the entries of the other domains, nor a
    Python object that is not a builtin, nor a name of Sphinx's own::

        sage: main['std:term']['object']
        'sphinx/object'
        sage: main['std:doc']['tutorial/index']
        'sphinx/tutorial'
        sage: main['std:label']['glossary']
        'sphinx/glossary'
        sage: main['c:function']['PyType_GenericAlloc']
        'sphinx/alloc'
        sage: main['py:class']['os.PathLike']
        'sphinx/pathlike'
        sage: main['py:class']['sphinx.application.Sphinx']
        'sphinx/app'
    """
    inventories = intersphinx.InventoryAdapter(app.env)
    for objtype, objects in inventories.named_inventory.get('python', {}).items():
        if not objtype.startswith('py:'):
            continue
        builtin_objects = {name: entry for name, entry in objects.items()
                           if hasattr(builtins, name)}
        if builtin_objects:
            inventories.main_inventory.setdefault(objtype, {}).update(builtin_objects)


dangling_debug = False


def debug_inf(app, message):
    if dangling_debug:
        logger.info(message)


def call_intersphinx(app, env, node, contnode):
    r"""
    Call intersphinx and make links between Sage manuals relative.

    TESTS:

    Check that the link from the thematic tutorials to the reference
    manual is relative, see :issue:`20118`::

        sage: from sage.env import SAGE_DOC
        sage: thematic_index = os.path.join(SAGE_DOC, "html", "en", "thematic_tutorials", "index.html")
        sage: for line in open(thematic_index).readlines():  # optional - sagemath_doc_html
        ....:     if "padics" in line:
        ....:         _ = sys.stdout.write(line)
        <li><p><a class="reference external" href="../reference/padics/sage/rings/padics/tutorial.html#sage-rings-padics-tutorial" title="(in $p$-adics v...)"><span>Introduction to the p-adics</span></a></p></li>
    """
    debug_inf(app, "???? Trying intersphinx for %s" % node['reftarget'])
    builder = app.builder
    res = intersphinx.missing_reference(
        app, env, node, contnode)
    if res:
        # Replace absolute links to $SAGE_DOC by relative links: this
        # allows to copy the whole documentation tree somewhere else
        # without breaking links, see Issue #20118.
        if res['refuri'].startswith(SAGE_DOC):
            here = os.path.dirname(os.path.join(builder.outdir,
                                                node['refdoc']))
            res['refuri'] = os.path.relpath(res['refuri'], here)
            debug_inf(app, "++++ Found at %s" % res['refuri'])
    else:
        debug_inf(app, "---- Intersphinx: %s not Found" % node['reftarget'])
    return res


# Objects that no importable module exposes under a documented name.  Only for
# cases that _public_alias() cannot derive; prefer fixing the reference itself.
_public_target_overrides = {
    # The workers of sage.parallel.map_reduce derive from the process class of
    # the 'fork' context, which multiprocessing documents as Process only.
    'multiprocessing.context.ForkProcess': 'multiprocessing.Process',
}


def _public_aliases(reftarget):
    r"""
    Return the public paths of the object with the dotted path ``reftarget``.

    Implementation modules are usually absent from the inventory of a project,
    while ``__module__`` points to them; the documented path is the one where a
    package re-exports the object.  Several packages on the way may re-export
    it and only one of them document it, so every alias is answered, the
    shortest path first.  Only modules that are imported already are inspected,
    and an alias must be the very same object, so nothing here is guessed.

    EXAMPLES::

        sage: import pexpect, unittest
        sage: from sage_docbuild.ext.crossrefs import _public_aliases
        sage: _public_aliases('unittest.case.TestCase')
        ['unittest.TestCase']
        sage: _public_aliases('pexpect.pty_spawn.spawn')
        ['pexpect.spawn']

    Paths that are public already, or that name no imported object, have no
    alias::

        sage: _public_aliases('unittest.TestCase')
        []
        sage: _public_aliases('unittest.case.NoSuchClass')
        []
        sage: _public_aliases('TestCase')
        []

    Every package that re-exports the object is answered for::

        sage: import sage.rings.polynomial.multi_polynomial_libsingular as libsing
        sage: aliases = _public_aliases(
        ....:     'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular')
        sage: aliases == sorted(aliases, key=len)
        True
    """
    if reftarget in _public_target_overrides:
        return [_public_target_overrides[reftarget]]
    modname, _, attrname = reftarget.rpartition('.')
    module = sys.modules.get(modname)
    if module is None:
        return []
    obj = safe_getattr(module, attrname, None)
    if obj is None:
        return []
    parts = modname.split('.')
    return ['.'.join(parts[:i] + [attrname])
            for i in range(1, len(parts))
            if (ancestor := sys.modules.get('.'.join(parts[:i]))) is not None
            and safe_getattr(ancestor, attrname, None) is obj]


def _retitle(node, contnode, reftarget, alias):
    r"""
    Rename the link ``contnode`` from ``reftarget`` to ``alias``.

    The text of a reference that names no title of its own is the path that was
    written, so a reference resolved under another path would read as the
    implementation module that no page documents.  A text that says anything
    else was written that way on purpose and is left alone.

    EXAMPLES::

        sage: from docutils import nodes
        sage: from sage_docbuild.ext.crossrefs import _retitle
        sage: xref = {'refexplicit': False}
        sage: content = nodes.literal('', 'unittest.case.TestCase')
        sage: _retitle(xref, content, 'unittest.case.TestCase', 'unittest.TestCase')
        sage: content.astext()
        'unittest.TestCase'

        sage: content = nodes.literal('', 'the base class')
        sage: _retitle(xref, content, 'unittest.case.TestCase', 'unittest.TestCase')
        sage: content.astext()
        'the base class'

    Even a title equal to its target is deliberate when it was written
    explicitly, so it is not replaced::

        sage: xref = {'refexplicit': True}
        sage: content = nodes.literal('', 'unittest.case.TestCase')
        sage: _retitle(xref, content, 'unittest.case.TestCase', 'unittest.TestCase')
        sage: content.astext()
        'unittest.case.TestCase'

    Sphinx appends parentheses to the automatic title of callable roles::

        sage: xref = {'refexplicit': False}
        sage: content = nodes.literal('', 'unittest.case.TestCase()')
        sage: _retitle(xref, content, 'unittest.case.TestCase', 'unittest.TestCase')
        sage: content.astext()
        'unittest.TestCase()'
    """
    from docutils import nodes

    if node.get('refexplicit') or contnode is None:
        return
    title = contnode.astext()
    if title == reftarget:
        replacement = alias
    elif title == reftarget + '()':
        replacement = alias + '()'
    else:
        return
    contnode.clear()
    contnode += nodes.Text(replacement)


def _type_parameter_names(modname):
    r"""
    Return the names of the type parameters declared in the module ``modname``.

    A class carries the type parameters of its bases too: the members that a
    subclass inherits are documented in the module of the subclass, so their
    annotations name type parameters declared elsewhere.

    EXAMPLES::

        sage: from sage_docbuild.ext.crossrefs import _type_parameter_names
        sage: sorted(_type_parameter_names('sage.categories.homset'))
        ['CodomainElementT', 'DomainElementT']

    ``VermaModuleHomset`` is a ``Homset``, and inherits its type parameters::

        sage: import sage.algebras.lie_algebras.verma_module
        sage: sorted(_type_parameter_names(
        ....:     'sage.algebras.lie_algebras.verma_module'))
        ['CodomainElementT', 'DomainElementT']

        sage: sorted(_type_parameter_names('no.such.module'))
        []

    The answer for a module that is not imported is not remembered: importing
    it later is what gives it type parameters to find::

        sage: name = 'sage.categories.homset'
        sage: _type_parameter_names(name) == _type_parameter_names(name)
        True

    """
    if modname not in sys.modules:
        # A later call, once something has imported it, has to be free to
        # answer differently.
        return frozenset()
    return _type_parameters_of_imported_module(modname)


@functools.lru_cache(maxsize=256)
def _type_parameters_of_imported_module(modname):
    r"""
    Do the work of :func:`_type_parameter_names` once per imported module.

    Module declarations are stable during a documentation build.  The bounded
    cache contains only module names and immutable answers; missing modules are
    filtered by :func:`_type_parameter_names` before reaching it.
    """
    from sage.misc.lazy_import import LazyImport

    module = sys.modules[modname]

    # Reading anything off a lazy import performs it, which loads a module for
    # no reason and warns when it is deprecated.  One stands for a name that
    # another module defines, so it declares no type parameter of this one.
    def type_params(obj):
        if isinstance(obj, LazyImport):
            return set()
        return {tp.__name__ for tp in safe_getattr(obj, '__type_params__', ())}

    names = set()
    seen_bases = set()
    for member in list(vars(module).values()):
        if isinstance(member, LazyImport):
            continue  # stands for a name of another module, and unread so far
        if safe_getattr(member, '__module__', None) != modname:
            continue  # imported from elsewhere
        names |= type_params(member)
        if isinstance(member, type):
            for base in safe_getattr(member, '__mro__', ()):
                identity = id(base)
                if identity in seen_bases:
                    continue
                seen_bases.add(identity)
                names |= type_params(base)
                for attr in list(vars(base).values()):
                    names |= type_params(attr)
    return frozenset(names)


def _is_type_parameter(node, reftarget):
    r"""
    Return whether ``reftarget`` names a type parameter in the context of ``node``.

    A :pep:`695` type parameter is documented nowhere: CPython creates it with
    ``__module__`` set to ``typing``, and it is meaningful only in the scope
    declaring it, so it is rendered either bare or under a ``typing.`` prefix.

    EXAMPLES::

        sage: from sage_docbuild.ext.crossrefs import _is_type_parameter
        sage: node = {'py:module': 'sage.categories.homset', 'py:class': 'Homset'}
        sage: _is_type_parameter(node, 'DomainElementT')
        True
        sage: _is_type_parameter(node, 'typing.DomainElementT')
        True
        sage: _is_type_parameter(node, 'typing.Any')
        False
        sage: _is_type_parameter(node, 'Parent')
        False
    """
    name = reftarget.removeprefix('typing.')
    if '.' in name or (name != reftarget and hasattr(typing, name)):
        return False
    return name in _type_parameter_names(node.get('py:module'))


def find_sage_dangling_links(app, env, node, contnode):
    r"""
    Try to find dangling link in local module imports or all.py.

    TESTS:

    Globally unique suffixes in inventories are deliberately not guessed.  A
    reference must resolve directly or in its Python module/class context::

        sage: from sage_docbuild.ext import crossrefs
        sage: node = {'refdoc': 'index', 'refdomain': 'py', 'reftype': 'obj',
        ....:         'reftarget': 'spawn'}
        sage: inventory = {'py:classmethod': {
        ....:     'ptyprocess.PtyProcess.spawn': object()}}
        sage: env = type('Env', (), {'intersphinx_inventory': inventory})()
        sage: from unittest.mock import patch
        sage: with patch.object(crossrefs, 'call_intersphinx', return_value=None):
        ....:     result = crossrefs.find_sage_dangling_links(None, env, node, None)
        sage: result is None
        True

    A builtin that Python documents as a function is retried under the function
    role, and the role the document asked for is restored afterwards.  Sphinx
    resolves nothing here and warns about nothing either, so without the retry
    the reference would silently render as plain text::

        sage: node = {'refdoc': 'index', 'refdomain': 'py', 'reftype': 'class',
        ....:         'reftarget': 'staticmethod'}
        sage: def as_func(app, env, node, contnode):
        ....:     if node['reftype'] == 'func':
        ....:         return {'refuri': 'library/functions.html#staticmethod'}
        sage: with patch.object(crossrefs, 'call_intersphinx', as_func):
        ....:     crossrefs.find_sage_dangling_links(None, env, node, None)
        {'refuri': 'library/functions.html#staticmethod'}
        sage: node['reftype']
        'class'
    """
    debug_inf(app, "==================== find_sage_dangling_links ")

    reftype = node['reftype']
    reftarget = node['reftarget']
    try:
        doc = node['refdoc']
    except KeyError:
        debug_inf(app, "-- no refdoc in node %s" % node)
        return None

    debug_inf(app, "Searching %s from %s" % (reftarget, doc))

    res = call_intersphinx(app, env, node, contnode)
    if res:
        debug_inf(app, "++ DONE %s" % (res['refuri']))
        return res

    # Python's inventory has changed over time for basic builtins.  Try the
    # alternate role after the role requested by the document.
    if reftarget in base_class_as_func and reftype == 'class':
        node['reftype'] = 'func'
        res = call_intersphinx(app, env, node, contnode)
        node['reftype'] = reftype
        if res:
            debug_inf(app, "++ DONE %s" % (res['refuri']))
            return res

    if reftarget in base_func_as_class and reftype == 'func':
        node['reftype'] = 'class'
        res = call_intersphinx(app, env, node, contnode)
        node['reftype'] = reftype
        if res:
            debug_inf(app, "++ DONE %s" % (res['refuri']))
            return res

    # Some inherited third-party signatures render ``typing.Any`` as ``t.Any``
    # and classify it as a class, while Python documents it as data.
    if reftarget == 't.Any' and reftype == 'class':
        node['reftarget'] = 'typing.Any'
        node['reftype'] = 'data'
        res = call_intersphinx(app, env, node, contnode)
        node['reftarget'] = reftarget
        node['reftype'] = reftype
        if res:
            debug_inf(app, "++ DONE %s" % (res['refuri']))
            return res

    if node.get('refdomain') != 'py':  # not a python file
        return None

    # A class is commonly re-exported by the package documenting it, while
    # ``__module__`` names the implementation module, which no inventory knows.
    # Any of the packages re-exporting it may be the one that documents it.
    for alias in _public_aliases(reftarget):
        node['reftarget'] = alias
        res = call_intersphinx(app, env, node, contnode)
        node['reftarget'] = reftarget
        if res:
            debug_inf(app, "++ DONE %s" % (res['refuri']))
            _retitle(node, contnode, reftarget, alias)
            return res

    # Type parameters have no target to link to; keep them as plain text.
    if _is_type_parameter(node, reftarget):
        debug_inf(app, "++ type parameter %s" % reftarget)
        return contnode

    try:
        module = node['py:module']
        cls = node['py:class']
    except KeyError:
        debug_inf(app, "-- no module or class for :%s:%s" % (reftype,
                                                             reftarget))
        return None

    def module_of(module, name):
        """The module defining the attribute ``name`` of ``module``, if any."""
        obj = safe_getattr(module, name, None)
        return None if obj is None else safe_getattr(obj, '__module__', None)

    basename = reftarget.split(".")[0]
    target_module = module_of(sys.modules['sage.all'], basename)
    if target_module is not None:
        debug_inf(app, "++ found %s using sage.all in %s" % (basename, target_module))
    else:
        this_module = sys.modules.get(node['py:module'])
        if this_module is not None:
            target_module = module_of(this_module, basename)
            if target_module is None:
                debug_inf(app, "-- %s not found in sage.all or this module" % (basename))
                return None
            debug_inf(app, "++ found %s in this module" % (basename,))
    if target_module is None:
        target_module = ""
        debug_inf(app, "?? found in None !!!")

    newtarget = target_module+'.'+reftarget
    node['reftarget'] = newtarget

    # adapted  from sphinx/domains/python.py
    builder = app.builder
    searchmode = node.hasattr('refspecific') and 1 or 0
    matches = builder.env.domains['py'].find_obj(
        builder.env, module, cls, newtarget, reftype, searchmode)
    if not matches:
        debug_inf(app, "?? no matching doc for %s" % newtarget)
        return call_intersphinx(app, env, node, contnode)
    if len(matches) > 1:
        logger.warning('more than one target found for cross-reference %r: %s',
                       newtarget, ', '.join(match[0] for match in matches),
                       location=node)
    name, obj = matches[0]
    debug_inf(app, "++ match = %s %s" % (name, obj))

    from docutils import nodes
    newnode = nodes.reference('', '', internal=True)
    if name == target_module:
        newnode['refid'] = name
    else:
        newnode['refuri'] = builder.get_relative_uri(node['refdoc'], obj[0])
        newnode['refuri'] += '#' + name
        debug_inf(app, "++ DONE at URI %s" % (newnode['refuri']))
    newnode['reftitle'] = name
    newnode.append(contnode)
    return newnode


# Basic Python classes documented as functions in some old inventories, and
# builtins which are classes that Python documents as functions to this day.
# The latter need this list to be linked at all: Sphinx's own
# ``builtin_resolver()`` renders a ``:class:`` reference to any builtin class as
# plain text and suppresses the warning that nitpicky mode would otherwise
# raise, so a reference that no inventory resolves is invisible in the output.
base_class_as_func = [
    'bool', 'classmethod', 'complex', 'dict', 'enumerate', 'filter',
    'float', 'frozenset', 'int', 'list', 'map', 'object', 'reversed',
    'set', 'slice', 'staticmethod', 'str', 'tuple', 'type', 'zip']

# Basic Python functions documented as classes in modern inventories.
base_func_as_class = [
    'bool', 'complex', 'dict', 'float', 'frozenset', 'int', 'list',
    'object', 'range', 'set', 'slice', 'str', 'tuple', 'type']


def setup(app):
    """
    Register this extension with Sphinx.
    """
    # When building the standard docs, app.srcdir is set to SAGE_DOC_SRC +
    # 'LANGUAGE/DOCNAME'.
    if not app.srcdir.is_relative_to(SAGE_DOC_SRC):
        return {'parallel_read_safe': True}

    app.add_config_value('intersphinx_resolve_self', 'sagemath', False)
    app.add_config_value('intersphinx_mapping', {}, False)
    app.add_config_value('intersphinx_cache_limit', 5, False)
    app.add_config_value('intersphinx_disabled_reftypes', [], False)
    app.add_config_value('intersphinx_timeout', None, False)
    app.connect('config-inited', set_intersphinx_mappings)
    app.connect('builder-inited', intersphinx.load_mappings)
    app.connect('builder-inited', prefer_python_inventory, priority=600)
    # We do *not* fully initialize intersphinx since we call it by hand
    # in find_sage_dangling_links.
    #   app.connect('missing-reference', missing_reference)
    app.connect('missing-reference', find_sage_dangling_links)
    return {'parallel_read_safe': True}
