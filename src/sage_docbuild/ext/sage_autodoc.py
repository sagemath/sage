"""
Sage autodoc extension

This is :mod:`sphinx.ext.autodoc` extension modified for Sage objects.

The original header of :mod:`sphinx.ext.autodoc`:

    Extension to create automatic documentation from code docstrings.

    Automatically insert docstrings for functions, classes or whole modules into
    the doctree, thus avoiding duplication between docstrings and documentation
    for those who like elaborate docstrings.

This module is currently based on the legacy class-based implementation of
:mod:`sphinx.ext.autodoc` from Sphinx version 9.1.0.  Compare (do diff) with
the upstream source file
`sphinx/ext/autodoc/_legacy_class_based/_documenters.py
<https://github.com/sphinx-doc/sphinx/blob/v9.1.0/sphinx/ext/autodoc/_legacy_class_based/_documenters.py>`_.

In the source file of this module, major modifications are delimited by a pair
of comment dividers. To lessen maintenance burden, we aim at reducing the modifications.
What Sphinx keeps outside that file -- the directive options, the sentinels and
the signature regular expressions -- is imported from :mod:`sphinx.ext.autodoc`,
which exports them publicly, so that a rebase only has to consider the parts
that Sage really changes.

.. NOTE::

    Since Sphinx 9.1.0, the class-based documenters that this module derives
    from are the opt-in ``autodoc_use_legacy_class_based`` implementation:
    Sphinx builds its own ``auto*`` directives on a different implementation by
    default.  Sage does not run into the two of them at once, because it loads
    this module in place of :mod:`sphinx.ext.autodoc`.

AUTHORS:

- ? (before 2011): initial version derived from sphinx.ext.autodoc

- Johan S. R. Nielsen: support for _sage_argspec_

- Simon King (2011-04): use sageinspect; include public cython attributes
  only in the documentation if they have a docstring

- Kwankyu Lee (2018-12-26, 2022-11-08): rebased on the latest sphinx.ext.autodoc

- Kwankyu Lee (2024-02-14): rebased on Sphinx 7.2.6

- François Bissey (2024-08-24): rebased on Sphinx 8.0.2

- François Bissey (2024-09-10): Tweaks to support python 3.9 (and older sphinx) as well

- François Bissey (2024-11-12): rebased on Sphinx 8.1.3 (while trying to keep python 3.9 compatibility)

- François Bissey (2025-02-24): Remove python 3.9 support hacks, making us closer to upstream

- François Bissey (2025-03-18): rebased on sphinx 8.2.3

- Chenxin Zhong (2026-07-07): rebased on Sphinx 9.1.0, where the class-based
  implementation this module derives from became ``autodoc_use_legacy_class_based``

- Chenxin Zhong (2026-07-25): import the directive options and the sentinels
  from :mod:`sphinx.ext.autodoc` instead of copying them
"""

from __future__ import annotations

import ast
import functools
import importlib.util
import operator
import re
import sys
import tokenize
from inspect import Parameter, Signature
from types import MappingProxyType, ModuleType
from typing import TYPE_CHECKING, Any, NewType, TypeVar

import sphinx
from docutils.statemachine import StringList
from sphinx.config import ENUM
from sphinx.errors import PycodeError

# Sphinx keeps these in modules of its own; import them rather than carrying a
# copy, so that only the parts that Sage really changes have to be rebased.
from sphinx.ext.autodoc import (
    ALL,
    INSTANCEATTR,
    SLOTSATTR,
    SUPPRESS,
    UNINITIALIZED_ATTR,
    annotation_option,
    bool_option,
    class_doc_from_option,
    exclude_members_option,
    identity,
    inherited_members_option,
    member_order_option,
    members_option,
    merge_members_option,
    py_ext_sig_re,
    special_member_re,
)
from sphinx.ext.autodoc.importer import get_class_members, import_module, import_object
from sphinx.ext.autodoc.mock import ismock, mock, undecorate
from sphinx.locale import _, __
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import inspect, logging
from sphinx.util.docstrings import prepare_docstring, separate_metadata
from sphinx.util.inspect import (
    TypeAliasForwardRef,
    TypeAliasNamespace,
    evaluate_signature,
    object_description,
    safe_getattr,
)
from sphinx.util.inspect import (
    stringify_signature as _sphinx_stringify_signature,
)
from sphinx.util.typing import (
    get_type_hints,
)
from sphinx.util.typing import (
    restify as _sphinx_restify,
)
from sphinx.util.typing import (
    stringify_annotation as _sphinx_stringify_annotation,
)

# ------------------------------------------------------------------
from sage.misc.sageinspect import (
    is_function_or_cython_function,
    isclassinstance,
    sage_formatargspec,
    sage_getargspec,
    sage_getdoc_original,
)


def getdoc(obj, *args, **kwargs):
    """Read a docstring the way Sage does, in place of Sphinx's own reader.

    This accepts the calling convention of
    ``sphinx.util.inspect.getdoc`` rather than repeating its signature,
    so that the documenters below can call it either way round.  Everything
    past the object is discarded: Sage resolves inheritance and Cython
    docstrings itself, so the ``attrgetter``, ``allow_inherited``, ``cls``
    and ``name`` that Sphinx would use have nothing left to decide.
    """
    return sage_getdoc_original(obj)


# ------------------------------------------------------------------

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator, Mapping, Sequence
    from typing import ClassVar, Literal

    from sphinx.application import Sphinx
    from sphinx.environment import BuildEnvironment, _CurrentDocument
    from sphinx.events import EventManager
    from sphinx.registry import SphinxComponentRegistry
    from sphinx.util.typing import ExtensionMetadata, _RestifyMode

logger = logging.getLogger(__name__)


def _get_render_mode(
    typehints_format: Literal['fully-qualified', 'short'],
) -> _RestifyMode:
    if typehints_format == 'short':
        return 'smart'
    return 'fully-qualified-except-typing'


# ------------------------------------------------------------------
# Annotations commonly refer to names that their module imports only under
# ``if TYPE_CHECKING:``.  Sphinx evaluates the annotations of a signature in
# one go, so a single such name makes the *whole* signature fall back to its
# unevaluated source text, where every name loses its module and no longer
# cross-references.  We therefore collect the deferred imports of a module and
# pass them to Sphinx as :confval:`autodoc_type_aliases` entries.

_NO_TYPE_ALIASES: Mapping[str, str] = MappingProxyType({})


class _ModuleBindingVisitor(ast.NodeVisitor):
    """Find names that a statement may bind in its surrounding module."""

    def __init__(self) -> None:
        self.names: set[str] = set()
        self.has_star_import = False

    def visit_Name(self, node: ast.Name) -> None:
        if isinstance(node.ctx, (ast.Store, ast.Del)):
            self.names.add(node.id)

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            self.names.add(alias.asname or alias.name.partition('.')[0])

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        for alias in node.names:
            if alias.name == '*':
                self.has_star_import = True
            else:
                self.names.add(alias.asname or alias.name)

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self.names.add(node.name)
        for decorator in node.decorator_list:
            self.visit(decorator)
        for default in (*node.args.defaults,
                        *(value for value in node.args.kw_defaults
                          if value is not None)):
            self.visit(default)
        arguments = (*node.args.posonlyargs, *node.args.args,
                     *node.args.kwonlyargs)
        for argument in arguments:
            if argument.annotation is not None:
                self.visit(argument.annotation)
        for argument in (node.args.vararg, node.args.kwarg):
            if argument is not None and argument.annotation is not None:
                self.visit(argument.annotation)
        if node.returns is not None:
            self.visit(node.returns)
        for parameter in getattr(node, 'type_params', ()):
            for attribute in ('bound', 'default_value'):
                value = getattr(parameter, attribute, None)
                if value is not None:
                    self.visit(value)

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self.names.add(node.name)
        for expression in (*node.decorator_list, *node.bases,
                           *(keyword.value for keyword in node.keywords)):
            self.visit(expression)
        for parameter in getattr(node, 'type_params', ()):
            for attribute in ('bound', 'default_value'):
                value = getattr(parameter, attribute, None)
                if value is not None:
                    self.visit(value)

    def visit_Lambda(self, node: ast.Lambda) -> None:
        for default in (*node.args.defaults,
                        *(value for value in node.args.kw_defaults
                          if value is not None)):
            self.visit(default)

    def _visit_comprehension(self, node, values) -> None:
        # Iteration targets live in the comprehension's implicit function
        # scope on Python 3.  Its iterable, filters, and result expressions can
        # still contain assignment expressions that bind in the surrounding
        # module, so visit those without visiting ``generator.target``.
        for generator in node.generators:
            self.visit(generator.iter)
            for condition in generator.ifs:
                self.visit(condition)
        for value in values:
            self.visit(value)

    def visit_ListComp(self, node: ast.ListComp) -> None:
        self._visit_comprehension(node, (node.elt,))

    visit_SetComp = visit_ListComp
    visit_GeneratorExp = visit_ListComp

    def visit_DictComp(self, node: ast.DictComp) -> None:
        self._visit_comprehension(node, (node.key, node.value))

    def visit_ExceptHandler(self, node: ast.ExceptHandler) -> None:
        if node.name:
            self.names.add(node.name)
        self.generic_visit(node)

    def visit_MatchAs(self, node: ast.MatchAs) -> None:
        if node.name:
            self.names.add(node.name)
        self.generic_visit(node)

    def visit_MatchStar(self, node: ast.MatchStar) -> None:
        if node.name:
            self.names.add(node.name)


def _bound_module_names(node: ast.AST) -> tuple[set[str], bool]:
    """Return names that *node* may bind, without entering nested scopes."""
    visitor = _ModuleBindingVisitor()
    visitor.visit(node)
    return visitor.names, visitor.has_star_import


def _type_checking_value(test: ast.expr, flags: set[str],
                         modules: set[str]) -> str | None:
    """Classify *test* as the typing flag or module, if it names either."""
    if isinstance(test, ast.Name):
        if test.id in flags:
            return 'flag'
        if test.id in modules:
            return 'module'
    if (isinstance(test, ast.Attribute) and test.attr == 'TYPE_CHECKING'
            and isinstance(test.value, ast.Name) and test.value.id in modules):
        return 'flag'
    return None


def _update_type_checking_names(node: ast.stmt, flags: set[str],
                                modules: set[str]) -> None:
    """Update known typing bindings after the module executes *node*."""
    def forget(names: set[str]) -> None:
        flags.difference_update(names)
        modules.difference_update(names)

    if (isinstance(node, ast.If)
            and _tests_type_checking(node.test, flags, modules)):
        # typing.TYPE_CHECKING is false while the module executes.  Bindings
        # in the positive body therefore do not replace the runtime aliases
        # used by a later guard; an ``else`` body does execute.
        for statement in node.orelse:
            _update_type_checking_names(statement, flags, modules)
        return

    if isinstance(node, ast.Import):
        for alias in node.names:
            name = alias.asname or alias.name.partition('.')[0]
            forget({name})
            if (alias.name == 'typing'
                    or (alias.asname is None
                        and alias.name.startswith('typing.'))):
                modules.add(name)
        return

    if isinstance(node, ast.ImportFrom):
        for alias in node.names:
            if alias.name == '*':
                # A star import may replace any tracked alias (``Any`` is a
                # particularly plausible one).  Only the canonical flag that
                # typing itself exports can be known afterwards.
                flags.clear()
                modules.clear()
                if node.level == 0 and node.module == 'typing':
                    flags.add('TYPE_CHECKING')
                continue
            name = alias.asname or alias.name
            forget({name})
            if (node.level == 0 and node.module == 'typing'
                    and alias.name == 'TYPE_CHECKING'):
                flags.add(name)
        return

    if isinstance(node, (ast.Assign, ast.AnnAssign)):
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if isinstance(node, ast.AnnAssign) and node.value is None:
            # A bare annotation does not bind its target at runtime.  Its
            # annotation expression can still contain a module-level named
            # expression, so invalidate only names bound there.
            forget(_bound_module_names(node.annotation)[0])
            return
        names = {name for target in targets
                 for name in _bound_module_names(target)[0]}
        if isinstance(node, ast.AnnAssign):
            names.update(_bound_module_names(node.annotation)[0])
        value = node.value
        if value is not None:
            # Assignment expressions in the value bind in the surrounding
            # scope too: ``holder = (TC := False)`` invalidates ``TC``.
            names.update(_bound_module_names(value)[0])
        kind = _type_checking_value(value, flags, modules) if value is not None else None
        forget(names)
        # Propagate a simple alias, but not a destructuring assignment whose
        # runtime value cannot be recovered from the syntax alone.
        simple_names = [target.id for target in targets if isinstance(target, ast.Name)]
        if kind == 'flag':
            flags.update(simple_names)
        elif kind == 'module':
            modules.update(simple_names)
        return

    names, has_star_import = _bound_module_names(node)
    if has_star_import:
        flags.clear()
        modules.clear()
    else:
        forget(names)


def _type_checking_guards(tree: ast.Module) -> Iterator[ast.If]:
    r"""
    Yield positive ``TYPE_CHECKING`` guards in their binding context.

    Imports take effect only after their statement, and bindings in nested
    scopes never leak into the module::

        sage: import ast
        sage: from sage_docbuild.ext.sage_autodoc import _type_checking_guards
        sage: source = '''\
        ....: from typing import TYPE_CHECKING as TC
        ....: if TC: pass
        ....: TC = False
        ....: if TC: pass
        ....: if LATE: pass
        ....: from typing import TYPE_CHECKING as LATE
        ....: def f():
        ....:     from typing import TYPE_CHECKING as LOCAL
        ....: if LOCAL: pass
        ....: import typing as t
        ....: if t.TYPE_CHECKING: pass
        ....: '''
        sage: [ast.unparse(guard.test)
        ....:  for guard in _type_checking_guards(ast.parse(source))]
        ['TC', 't.TYPE_CHECKING']

    Comprehension iteration variables do not rebind the module's flag::

        sage: tree = ast.parse(
        ....:     'from typing import TYPE_CHECKING as TC\n'
        ....:     'dummy = [TC for TC in ()]\nif TC: pass')
        sage: [ast.unparse(guard.test) for guard in _type_checking_guards(tree)]
        ['TC']

    Expressions evaluated while defining a nested scope can rebind it::

        sage: tree = ast.parse(
        ....:     'from typing import TYPE_CHECKING as TC\n'
        ....:     'def f(value=(TC := False)): pass\nif TC: pass')
        sage: list(_type_checking_guards(tree))
        []

    A bare module annotation does not replace the value it annotates::

        sage: tree = ast.parse(
        ....:     'from typing import TYPE_CHECKING as TC\n'
        ....:     'TC: bool\nif TC: pass')
        sage: [ast.unparse(guard.test) for guard in _type_checking_guards(tree)]
        ['TC']

    Nor does an assignment in the non-executed body of an earlier typing
    guard replace the flag for a later guard::

        sage: tree = ast.parse(
        ....:     'from typing import TYPE_CHECKING as TC\n'
        ....:     'if TC:\n    TC = False\nif TC: pass')
        sage: [ast.unparse(guard.test) for guard in _type_checking_guards(tree)]
        ['TC', 'TC']

    A star import from :mod:`typing` binds the flag too::

        sage: tree = ast.parse('from typing import *\nif TYPE_CHECKING: pass')
        sage: [ast.unparse(guard.test) for guard in _type_checking_guards(tree)]
        ['TYPE_CHECKING']
        sage: tree = ast.parse(
        ....:     'from typing import TYPE_CHECKING as Any\n'
        ....:     'from typing import *\nif Any: pass')
        sage: list(_type_checking_guards(tree))
        []
        sage: tree = ast.parse(
        ....:     'from .typing import TYPE_CHECKING\nif TYPE_CHECKING: pass')
        sage: list(_type_checking_guards(tree))
        []

    Importing a :mod:`typing` submodule without an alias binds the top-level
    module, just as Python does::

        sage: tree = ast.parse(
        ....:     'import typing.io\nif typing.TYPE_CHECKING: pass')
        sage: [ast.unparse(guard.test) for guard in _type_checking_guards(tree)]
        ['typing.TYPE_CHECKING']
    """
    flags: set[str] = set()
    modules: set[str] = set()
    for node in tree.body:
        if (isinstance(node, ast.If)
                and _tests_type_checking(node.test, flags, modules)):
            yield node
        _update_type_checking_names(node, flags, modules)


def _tests_type_checking(test: ast.expr, flags: set[str],
                         modules: set[str]) -> bool:
    """
    Return whether the expression *test* is a guard on ``TYPE_CHECKING``.

    *flags* and *modules* are the bindings known immediately before *test*.
    Only the plain guard counts, and only on the flag of :mod:`typing`: a
    negated or a compound test binds its names under a condition of its own,
    and an unrelated ``TYPE_CHECKING`` - a setting of a framework, say - says
    nothing about what is deferred for typing.

    EXAMPLES::

        sage: import ast
        sage: from sage_docbuild.ext.sage_autodoc import _tests_type_checking
        sage: def guards(source, flags={'TYPE_CHECKING'}, modules={'typing'}):
        ....:     return _tests_type_checking(
        ....:         ast.parse(source, mode='eval').body, flags, modules)
        sage: guards('TYPE_CHECKING')
        True
        sage: guards('typing.TYPE_CHECKING')
        True
        sage: guards('sys.version_info')
        False
        sage: guards('not TYPE_CHECKING')
        False
        sage: guards('TYPE_CHECKING and sys.platform')
        False

    A module that never imported the flag does not test it::

        sage: guards('TYPE_CHECKING', flags=set())
        False
        sage: guards('settings.TYPE_CHECKING')
        False
    """
    if isinstance(test, ast.Name):
        return test.id in flags
    return (isinstance(test, ast.Attribute) and test.attr == 'TYPE_CHECKING'
            and isinstance(test.value, ast.Name) and test.value.id in modules)


def _type_checking_aliases(modname: str) -> Mapping[str, str]:
    r"""
    Return the names that the module *modname* imports under ``TYPE_CHECKING``.

    The names are mapped to their fully qualified targets, in the format
    expected of :confval:`autodoc_type_aliases`.  Names that the module also
    provides at runtime are left out, since those already evaluate.

    EXAMPLES::

        sage: from sage_docbuild.ext.sage_autodoc import _type_checking_aliases
        sage: aliases = _type_checking_aliases('sage_docbuild.ext.sage_autodoc')
        sage: aliases['Sphinx']
        'sphinx.application.Sphinx'
        sage: aliases['ExtensionMetadata']
        'sphinx.util.typing.ExtensionMetadata'

    Modules that are not imported, that have no Python source available, or
    that defer no import give an empty mapping::

        sage: sorted(_type_checking_aliases('no.such.module'))
        []
        sage: sorted(_type_checking_aliases('sage.rings.integer'))   # compiled
        []
        sage: sorted(_type_checking_aliases('sage.misc.sageinspect'))
        []

    Looking for source metadata does not invoke a module's PEP 562 hook::

        sage: from types import ModuleType
        sage: dynamic = ModuleType('_sage_docbuild_dynamic_alias_test')
        sage: calls = []
        sage: dynamic.__getattr__ = lambda name: calls.append(name)
        sage: sys.modules[dynamic.__name__] = dynamic
        sage: dict(_type_checking_aliases(dynamic.__name__)), calls
        ({}, [])
        sage: del sys.modules[dynamic.__name__]

    The answer for a module that is not imported is not remembered: importing
    it later is what gives it one::

        sage: import sage.geometry.polyhedron.base  # any module with aliases
        sage: name = 'sage.geometry.polyhedron.base'
        sage: _type_checking_aliases(name) == _type_checking_aliases(name)
        True

    Parsed source is cached by its stable module metadata, while names
    appearing in the live module namespace are filtered on every call::

        sage: from pathlib import Path
        sage: from sage.misc.temporary_file import tmp_filename
        sage: name = '_sage_docbuild_type_alias_test'
        sage: filename = tmp_filename(ext='.py')
        sage: source = '''\
        ....: from typing import TYPE_CHECKING
        ....: if TYPE_CHECKING:
        ....:     if MAYBE:
        ....:         from wrong import Conditional
        ....:     from example import Deferred
        ....:     Deferred: object
        ....: '''
        sage: _ = Path(filename).write_text(source)
        sage: module = ModuleType(name)
        sage: module.__file__ = filename
        sage: sys.modules[name] = module
        sage: dict(_type_checking_aliases(name))
        {'Deferred': 'example.Deferred'}
        sage: module.Deferred = object()
        sage: dict(_type_checking_aliases(name))
        {}
        sage: del sys.modules[name]
        sage: Path(filename).unlink()
    """
    module = sys.modules.get(modname)
    if module is None:
        # Reading it needs __file__, and a later call, once something has
        # imported it, has to be free to answer differently.
        return _NO_TYPE_ALIASES
    if not isinstance(module, ModuleType):
        return _NO_TYPE_ALIASES
    namespace = vars(module)
    filename = namespace.get('__file__') or ''
    if not filename.endswith('.py'):
        return _NO_TYPE_ALIASES  # compiled: no source to scan
    module_name = namespace.get('__name__') or modname
    package = namespace.get('__package__') or module_name.rpartition('.')[0]
    aliases = _aliases_of_imported_module(filename, package)
    if not aliases or not any(name in namespace for name in aliases):
        return aliases
    return MappingProxyType({name: target for name, target in aliases.items()
                             if name not in namespace})


@functools.lru_cache(maxsize=256)
def _aliases_of_imported_module(
    filename: str,
    package: str,
) -> Mapping[str, str]:
    """
    Parse deferred aliases once per source file.

    Source contents, module objects, and module namespaces are deliberately
    absent from the bounded cache key.  A running docbuild does not edit the
    source of an imported module, while runtime names are filtered cheaply by
    :func:`_type_checking_aliases` after this returns.
    """
    try:
        with tokenize.open(filename) as source_file:
            source = source_file.read()
        tree = ast.parse(source, filename=filename)
    except (OSError, SyntaxError, UnicodeError, ValueError):
        return _NO_TYPE_ALIASES

    aliases: dict[str, str] = {}
    for guard in _type_checking_guards(tree):
        # Only ``from ... import ...`` is handled: a plain ``import x`` binds a
        # module, which Sphinx cannot look attributes up on.
        # A nested conditional has terms of its own, so only imports that the
        # positive guard executes directly are unconditional typing aliases.
        for node in guard.body:
            if not isinstance(node, ast.ImportFrom):
                if isinstance(node, ast.AnnAssign) and node.value is None:
                    names, has_star_import = _bound_module_names(
                        node.annotation)
                else:
                    names, has_star_import = _bound_module_names(node)
                if has_star_import:
                    aliases.clear()
                else:
                    for name in names:
                        aliases.pop(name, None)
                continue
            if any(alias.name == '*' for alias in node.names):
                aliases.clear()
                continue
            if node.level:
                try:
                    target = importlib.util.resolve_name(
                        '.' * node.level + (node.module or ''), package
                    )
                except (ImportError, ValueError):
                    continue
            else:
                target = node.module or ''
            if not target:
                continue
            for alias in node.names:
                name = alias.asname or alias.name
                aliases.pop(name, None)
                if name == '*':
                    continue
                aliases[name] = f'{target}.{alias.name}'
    return MappingProxyType(aliases)


def _unwrap_type_alias_forward_refs(text: str) -> str:
    """
    Strip ``TypeAliasForwardRef`` wrappers from an already-rendered annotation.

    Sphinx resolves the values of :confval:`autodoc_type_aliases` to
    :class:`~sphinx.util.inspect.TypeAliasForwardRef` objects and unwraps them
    again for the annotation of a parameter or of a return value.  Inside a
    composite annotation such as ``Alias | None`` the wrapper survives, and
    Sphinx renders it as its :func:`repr`.

    EXAMPLES::

        sage: from sage_docbuild.ext.sage_autodoc import _unwrap_type_alias_forward_refs
        sage: _unwrap_type_alias_forward_refs(
        ....:     "TypeAliasForwardRef('sphinx.application.Sphinx') | None")
        'sphinx.application.Sphinx | None'

    Quoted literal values may force :func:`repr` to escape the quote that wraps
    the forward-reference name.  Decode that string literal rather than leaving
    its representation escapes in the displayed annotation::

        sage: from sphinx.util.inspect import TypeAliasForwardRef
        sage: name = '''typing.Literal['left"right']'''
        sage: rendered = f'{TypeAliasForwardRef(name)!r} | None'
        sage: _unwrap_type_alias_forward_refs(rendered) == f'{name} | None'
        True

    Wrapper-looking text that is itself quoted is left alone::

        sage: text = '''typing.Literal["TypeAliasForwardRef('value')"]'''
        sage: _unwrap_type_alias_forward_refs(text) == text
        True
    """
    pattern = r'''TypeAliasForwardRef\((?P<literal>'(?:\\.|[^'\\])*'|"(?:\\.|[^"\\])*")\)'''
    wrapper = re.compile(pattern)
    answer = []
    quote = None
    position = 0
    while position < len(text):
        character = text[position]
        if quote is not None:
            answer.append(character)
            if character == '\\' and position + 1 < len(text):
                position += 1
                answer.append(text[position])
            elif character == quote:
                quote = None
            position += 1
            continue
        if character in "'\"":
            quote = character
            answer.append(character)
            position += 1
            continue
        match = wrapper.match(text, position)
        if match is None:
            answer.append(character)
            position += 1
            continue
        try:
            name = ast.literal_eval(match['literal'])
        except (SyntaxError, ValueError):
            name = None
        answer.append(name if isinstance(name, str) else match[0])
        position = match.end()
    return ''.join(answer)


def _render_subscript_argument(argument: Any) -> str:
    """
    Render *argument* as it was written between the brackets of an annotation.

    :func:`stringify_annotation` renders a type; the brackets of an annotation
    hold other things too, and each of them has a source form of its own.

    EXAMPLES::

        sage: from sage_docbuild.ext.sage_autodoc import _render_subscript_argument
        sage: _render_subscript_argument(int)
        'int'
        sage: _render_subscript_argument('literal')            # as in Literal
        "'literal'"
        sage: _render_subscript_argument(())                   # as in Tuple[()]
        '()'
        sage: _render_subscript_argument((int,))
        '(int,)'
        sage: _render_subscript_argument([int, str])            # as in Callable
        '[int, str]'
        sage: _render_subscript_argument({'item': int})         # as in Annotated
        "{'item': int}"
        sage: _render_subscript_argument(slice(None, 3, None))  # as in A[:3]
        ':3'
        sage: _render_subscript_argument(Ellipsis)
        '...'
    """
    if isinstance(argument, TypeAliasForwardRef):
        return argument.name
    if isinstance(argument, (str, bytes)):
        # A Literal holds values, and a string among them keeps its quotes.
        return repr(argument)
    if isinstance(argument, tuple):
        inside = ', '.join(map(_render_subscript_argument, argument))
        if len(argument) == 1:
            inside += ','
        return f'({inside})'
    if isinstance(argument, list):
        return '[' + ', '.join(map(_render_subscript_argument, argument)) + ']'
    if isinstance(argument, dict):
        items = (f'{_render_subscript_argument(key)}: '
                 f'{_render_subscript_argument(value)}'
                 for key, value in argument.items())
        return '{' + ', '.join(items) + '}'
    if isinstance(argument, (set, frozenset)):
        items = sorted(map(_render_subscript_argument, argument))
        if isinstance(argument, frozenset):
            return ('frozenset({' + ', '.join(items) + '})'
                    if items else 'frozenset()')
        return '{' + ', '.join(items) + '}' if items else 'set()'
    if isinstance(argument, slice):
        parts = [argument.start, argument.stop, argument.step]
        while len(parts) > 2 and parts[-1] is None:
            parts.pop()
        return ':'.join('' if part is None
                        else _render_subscript_argument(part) for part in parts)
    return stringify_annotation(argument)


def _render_literal_argument(argument: Any) -> str:
    """Render a value in a deferred :class:`typing.Literal`."""
    from sphinx.util.inspect import isenumattribute

    if isinstance(argument, TypeAliasForwardRef):
        return argument.name
    if isenumattribute(argument):
        enum = type(argument)
        return f'{enum.__module__}.{enum.__qualname__}.{argument.name}'
    return repr(argument)


def _subscript_type_alias_forward_ref(self: TypeAliasForwardRef, item: Any) -> Any:
    """
    Return the subscripted forward reference, as another forward reference.

    An annotation such as ``list[ReferenceType[Expect]]`` subscripts the names
    it is made of, and the class that Sphinx substitutes for an entry of
    :confval:`autodoc_type_aliases` supports ``|`` but not ``[]``.  Evaluating
    such an annotation raises, and Sphinx then keeps the whole signature as the
    source text it was written as, where no name cross-references any more.

    EXAMPLES::

        sage: from sphinx.util.inspect import TypeAliasForwardRef
        sage: import sage_docbuild.ext.sage_autodoc
        sage: TypeAliasForwardRef('weakref.ReferenceType')[int]
        TypeAliasForwardRef('weakref.ReferenceType[int]')
        sage: TypeAliasForwardRef('collections.abc.Callable')[..., int]
        TypeAliasForwardRef('collections.abc.Callable[..., int]')

    The forms that a subscription can take are kept as they were written::

        sage: TypeAliasForwardRef('typing.Tuple')[()]
        TypeAliasForwardRef('typing.Tuple[()]')
        sage: TypeAliasForwardRef('typing.Literal')['left', 'right']
        TypeAliasForwardRef("typing.Literal['left', 'right']")
        sage: TypeAliasForwardRef('collections.abc.Callable')[[int, str], int]
        TypeAliasForwardRef('collections.abc.Callable[[int, str], int]')
        sage: TypeAliasForwardRef('example.Alias')[int,]
        TypeAliasForwardRef('example.Alias[int,]')
        sage: TypeAliasForwardRef('example.Alias')[((int,), str)]
        TypeAliasForwardRef('example.Alias[(int,), str]')

    Enum members keep the qualified source form required by
    :class:`typing.Literal`::

        sage: from enum import Enum
        sage: class Color(Enum):
        ....:     red = 1
        sage: deferred = TypeAliasForwardRef('typing.Literal')[Color.red]
        sage: deferred.name.endswith('.Color.red]') and '<Color.red:' not in deferred.name
        True
    """
    render = (_render_literal_argument
              if self.name.rpartition('.')[2] == 'Literal'
              else _render_subscript_argument)
    if isinstance(item, tuple):
        if not item:
            # The empty tuple is the item of A[()], not an empty subscription.
            inside = '()'
        else:
            inside = ', '.join(map(render, item))
            # __class_getitem__ receives a tuple, rather than its sole member,
            # for a one-item comma-separated subscription.
            if (len(item) == 1
                    and not (isinstance(item[0], TypeAliasForwardRef)
                             and item[0].name.startswith('*'))):
                inside += ','
    else:
        inside = render(item)
    return TypeAliasForwardRef(f'{self.name}[{inside}]')


def _iterate_type_alias_forward_ref(
    self: TypeAliasForwardRef,
) -> Iterator[TypeAliasForwardRef]:
    r"""
    Yield the unpacked form of a forward reference exactly once.

    Merely adding ``__getitem__`` would otherwise activate Python's legacy
    sequence iteration, which keeps asking for indices forever because every
    integer is a valid annotation parameter.  A one-item iterator also models
    the behavior of :class:`typing.TypeVarTuple` in ``Alias[*Parameters]``::

        sage: from sphinx.util.inspect import TypeAliasForwardRef
        sage: import sage_docbuild.ext.sage_autodoc
        sage: parameters = TypeAliasForwardRef('example.Parameters')
        sage: list(parameters)
        [TypeAliasForwardRef('*example.Parameters')]
        sage: TypeAliasForwardRef('example.Alias')[*parameters]
        TypeAliasForwardRef('example.Alias[*example.Parameters]')
    """
    yield TypeAliasForwardRef(f'*{self.name}')


if not hasattr(TypeAliasForwardRef, '__getitem__'):
    TypeAliasForwardRef.__getitem__ = _subscript_type_alias_forward_ref
if not hasattr(TypeAliasForwardRef, '__iter__'):
    TypeAliasForwardRef.__iter__ = _iterate_type_alias_forward_ref


def restify(cls: Any, mode: _RestifyMode = 'fully-qualified-except-typing') -> str:
    """
    Convert a type-like object to a reST reference.
    """
    return _unwrap_type_alias_forward_refs(_sphinx_restify(cls, mode=mode))


def stringify_annotation(
    annotation: Any,
    /,
    mode: _RestifyMode = 'fully-qualified-except-typing',
    *,
    short_literals: bool = False,
) -> str:
    """
    Stringify a type annotation.
    """
    return _unwrap_type_alias_forward_refs(
        _sphinx_stringify_annotation(
            annotation, mode=mode, short_literals=short_literals
        )
    )


def stringify_signature(
    sig: Signature,
    show_annotation: bool = True,
    show_return_annotation: bool = True,
    unqualified_typehints: bool = False,
    short_literals: bool = False,
) -> str:
    """
    Stringify a signature.
    """
    return _unwrap_type_alias_forward_refs(
        _sphinx_stringify_signature(
            sig,
            show_annotation=show_annotation,
            show_return_annotation=show_return_annotation,
            unqualified_typehints=unqualified_typehints,
            short_literals=short_literals,
        )
    )


# ------------------------------------------------------------------




























# Some useful event listener factories for autodoc-process-docstring.








class ObjectMember:
    """A member of object.

    This is used for the result of `Documenter.get_module_members()` to
    represent each member of the object.
    """

    __slots__ = '__name__', 'object', 'docstring', 'class_', 'skipped'

    __name__: str
    object: Any
    docstring: str | None
    class_: Any
    skipped: bool

    def __init__(
        self,
        name: str,
        obj: Any,
        *,
        docstring: str | None = None,
        class_: Any = None,
        skipped: bool = False,
    ) -> None:
        self.__name__ = name
        self.object = obj
        self.docstring = docstring
        self.class_ = class_
        self.skipped = skipped

    def __repr__(self) -> str:
        return (
            f'ObjectMember('
            f'name={self.__name__!r}, '
            f'obj={self.object!r}, '
            f'docstring={self.docstring!r}, '
            f'class_={self.class_!r}, '
            f'skipped={self.skipped!r}'
            f')'
        )


class Documenter:
    """A Documenter knows how to autodocument a single object type.  When
    registered with the AutoDirective, it will be used to document objects
    of that type when needed by autodoc.

    Its *objtype* attribute selects what auto directive it is assigned to
    (the directive name is 'auto' + objtype), and what directive it generates
    by default, though that can be overridden by an attribute called
    *directivetype*.

    A Documenter has an *option_spec* that works like a docutils directive's;
    in fact, it will be used to parse an auto directive's options that matches
    the Documenter.
    """

    #: name by which the directive is called (auto...) and the default
    #: generated directive name
    objtype = 'object'
    #: indentation by which to indent the directive content
    content_indent = '   '
    #: priority if multiple documenters return True from can_document_member
    priority = 0
    #: order if autodoc_member_order is set to 'groupwise'
    member_order = 0
    #: true if the generated content may contain titles
    titles_allowed = True

    option_spec: ClassVar[dict[str, Any]] = {
        'no-index': bool_option,
        'no-index-entry': bool_option,
        'noindex': bool_option,
    }

    def get_attr(self, obj: Any, name: str, *defargs: Any) -> Any:
        """getattr() override for types such as Zope interfaces."""
        return autodoc_attrgetter(obj, name, *defargs, registry=self.env._registry)

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        """Called to see if a member can be documented by this Documenter."""
        msg = 'must be implemented in subclasses'
        raise NotImplementedError(msg)

    def __init__(
        self, directive: Any, name: str, indent: str = ''
    ) -> None:
        self.directive = directive
        self.config: sphinx.config.Config = directive.env.config
        self.env: BuildEnvironment = directive.env
        self._current_document: _CurrentDocument = directive.env.current_document
        self._events: EventManager = directive.env.events
        self.options = directive.genopt
        self.name = name
        self.indent = indent
        # the module and object path within the module, and the fully
        # qualified name (all set after resolve_name succeeds)
        self.modname: str = ''
        self.module: ModuleType | None = None
        self.objpath: list[str] = []
        self.fullname = ''
        # extra signature items (arguments and return annotation,
        # also set after resolve_name succeeds)
        self.args: str | None = None
        self.retann: str = ''
        # the object to document (set after import_object succeeds)
        self.object: Any = None
        self.object_name = ''
        # the parent/owner of the object to document
        self.parent: Any = None
        # the module analyzer to get at attribute docs, or None
        self.analyzer: ModuleAnalyzer | None = None

    @property
    def documenters(self) -> dict[str, type[Documenter]]:
        """Returns registered Documenter classes"""
        return self.env._registry.documenters

    def add_line(self, line: str, source: str, *lineno: int) -> None:
        """Append one line of generated reST to the output."""
        if line.strip():  # not a blank line
            self.directive.result.append(self.indent + line, source, *lineno)
        else:
            self.directive.result.append('', source, *lineno)

    def resolve_name(
        self, modname: str | None, parents: Any, path: str, base: str
    ) -> tuple[str | None, list[str]]:
        """Resolve the module and name of the object to document given by the
        arguments and the current module/class.

        Must return a pair of the module name and a chain of attributes; for
        example, it would return ``('zipfile', ['ZipFile', 'open'])`` for the
        ``zipfile.ZipFile.open`` method.
        """
        msg = 'must be implemented in subclasses'
        raise NotImplementedError(msg)

    def parse_name(self) -> bool:
        """Determine what module to import and what attribute to document.

        Returns True and sets *self.modname*, *self.objpath*, *self.fullname*,
        *self.args* and *self.retann* if parsing and resolving was successful.
        """
        # first, parse the definition -- auto directives for classes and
        # functions can contain a signature which is then used instead of
        # an autogenerated one
        matched = py_ext_sig_re.match(self.name)
        if matched is None:
            logger.warning(
                __('invalid signature for auto%s (%r)'),
                self.objtype,
                self.name,
                type='autodoc',
            )
            return False
        explicit_modname, path, base, tp_list, args, retann = matched.groups()

        # support explicit module and class name separation via ::
        if explicit_modname is not None:
            modname = explicit_modname[:-2]
            parents = path.rstrip('.').split('.') if path else []
        else:
            modname = None
            parents = []

        with mock(self.config.autodoc_mock_imports):
            modname, self.objpath = self.resolve_name(modname, parents, path, base)

        if not modname:
            return False

        self.modname = modname
        self.args = args
        self.retann = retann
        self.fullname = '.'.join((self.modname or '', *self.objpath))
        return True

    def import_object(self, raiseerror: bool = False) -> bool:
        """Import the object given by *self.modname* and *self.objpath* and set
        it as *self.object*.

        Returns True if successful, False if an error occurred.
        """
        with mock(self.config.autodoc_mock_imports):
            try:
                ret = import_object(
                    self.modname, self.objpath, self.objtype, attrgetter=self.get_attr
                )
                self.module, self.parent, self.object_name, self.object = ret
                if ismock(self.object):
                    self.object = undecorate(self.object)
                return True
            except ImportError as exc:
                if raiseerror:
                    raise
                logger.warning(exc.args[0], type='autodoc', subtype='import_object')
                self.env.note_reread()
                return False

    def get_real_modname(self) -> str:
        """Get the real module name of an object to document.

        It can differ from the name of the module through which the object was
        imported.
        """
        return self.get_attr(self.object, '__module__', None) or self.modname

    # ------------------------------------------------------------------
    def _type_aliases(self, obj: Any = None) -> Mapping[str, str]:
        """Get the type aliases to evaluate the annotations of *obj* with.

        These are the configured :confval:`autodoc_type_aliases`, augmented
        with the names that the module defining *obj* imports under
        ``if TYPE_CHECKING:``.  Without them, a single deferred name makes
        Sphinx keep the whole signature unevaluated.  *obj* defaults to the
        object being documented; pass a module or a class to evaluate the
        annotations of its members.

        This is for :func:`~sphinx.util.inspect.signature`, which builds a
        namespace of its own; see :meth:`_type_alias_namespace` for the rest.
        """
        modname = self.get_attr(obj if obj is not None else self.object,
                                '__module__', None) or self.modname
        aliases = _type_checking_aliases(modname)
        if not aliases:
            return self.config.autodoc_type_aliases
        # Configured aliases take precedence over the deferred imports.
        return {**aliases, **self.config.autodoc_type_aliases}

    def _type_alias_namespace(self, obj: Any = None) -> TypeAliasNamespace:
        """Get the local namespace to evaluate the annotations of *obj* in.

        The aliases of :meth:`_type_aliases` are plain strings, which
        :func:`~sphinx.util.typing.get_type_hints` and
        :func:`~sphinx.util.inspect.evaluate_signature` would evaluate a
        second time; that fails for a target whose module is not imported.
        """
        return TypeAliasNamespace(self._type_aliases(obj))

    # ------------------------------------------------------------------

    def check_module(self) -> bool:
        """Check if *self.object* is really defined in the module given by
        *self.modname*.
        """
        if self.options.imported_members:
            return True

        subject = inspect.unpartial(self.object)
        modname = self.get_attr(subject, '__module__', None)
        return not modname or modname == self.modname

    def format_args(self, **kwargs: Any) -> str:
        """Format the argument signature of *self.object*.

        Should return None if the object does not have a signature.
        """
        return ''

    def format_name(self) -> str:
        """Format the name of *self.object*.

        This normally should be something that can be parsed by the generated
        directive, but doesn't need to be (Sphinx will display it unparsed
        then).
        """
        # normally the name doesn't contain the module (except for module
        # directives of course)
        return '.'.join(self.objpath) or self.modname

    def _call_format_args(self, **kwargs: Any) -> str:
        if kwargs:
            try:
                return self.format_args(**kwargs)
            except TypeError:
                # avoid chaining exceptions, by putting nothing here
                pass

        # retry without arguments for old documenters
        return self.format_args()

    def format_signature(self, **kwargs: Any) -> str:
        """Format the signature (arguments and return annotation) of the object.

        Let the user process it via the ``autodoc-process-signature`` event.
        """
        if self.args is not None:
            # signature given explicitly
            args = f'({self.args})'
            retann = self.retann
        else:
            # try to introspect the signature
            try:
                retann = None
                args = self._call_format_args(**kwargs)
                if args:
                    matched = re.match(r'^(\(.*\))\s+->\s+(.*)$', args)
                    if matched:
                        args = matched.group(1)
                        retann = matched.group(2)
            except Exception as exc:
                logger.warning(
                    __('error while formatting arguments for %s: %s'),
                    self.fullname,
                    exc,
                    type='autodoc',
                )
                args = None

        result = self._events.emit_firstresult(
            'autodoc-process-signature',
            self.objtype,
            self.fullname,
            self.object,
            self.options,
            args,
            retann,
        )
        if result:
            args, retann = result

        if args is not None:
            return args + ((' -> %s' % retann) if retann else '')
        return ''

    def add_directive_header(self, sig: str) -> None:
        """Add the directive header and options to the generated content."""
        domain = getattr(self, 'domain', 'py')
        directive = getattr(self, 'directivetype', self.objtype)
        name = self.format_name()
        sourcename = self.get_sourcename()

        # one signature per line, indented by column
        prefix = f'.. {domain}:{directive}:: '
        for i, sig_line in enumerate(sig.split('\n')):
            self.add_line(f'{prefix}{name}{sig_line}', sourcename)
            if i == 0:
                prefix = ' ' * len(prefix)

        if self.options.no_index or self.options.noindex:
            self.add_line('   :no-index:', sourcename)
        if self.options.no_index_entry:
            self.add_line('   :no-index-entry:', sourcename)
        if self.objpath:
            # Be explicit about the module, this is necessary since .. class::
            # etc. don't support a prepended module name
            self.add_line('   :module: %s' % self.modname, sourcename)

    def get_doc(self) -> list[list[str]] | None:
        """Decode and return lines of the docstring(s) for the object.

        When it returns None, autodoc-process-docstring will not be called for this
        object.
        """
        docstring = getdoc(
            self.object,
            self.get_attr,
            self.config.autodoc_inherit_docstrings,
            self.parent,
            self.object_name,
        )
        if docstring:
            tab_width = self.directive.state.document.settings.tab_width
            return [prepare_docstring(docstring, tab_width)]
        return []

    def process_doc(self, docstrings: list[list[str]]) -> Iterator[str]:
        """Let the user process the docstrings before adding them."""
        for docstringlines in docstrings:
            if self._events is not None:
                # let extensions preprocess docstrings
                self._events.emit(
                    'autodoc-process-docstring',
                    self.objtype,
                    self.fullname,
                    self.object,
                    self.options,
                    docstringlines,
                )

                if docstringlines and docstringlines[-1]:
                    # append a blank line to the end of the docstring
                    docstringlines.append('')

            yield from docstringlines

    def get_sourcename(self) -> str:
        obj_module = inspect.safe_getattr(self.object, '__module__', None)
        obj_qualname = inspect.safe_getattr(self.object, '__qualname__', None)
        if obj_module and obj_qualname:
            # Get the correct location of docstring from self.object
            # to support inherited methods
            fullname = f'{self.object.__module__}.{self.object.__qualname__}'
        else:
            fullname = self.fullname

        if self.analyzer:
            return f'{self.analyzer.srcname}:docstring of {fullname}'
        return 'docstring of %s' % fullname

    def add_content(self, more_content: StringList | None) -> None:
        """Add content from docstrings, attribute documentation and user."""
        docstring = True

        # set sourcename and add content from attribute documentation
        sourcename = self.get_sourcename()
        if self.analyzer:
            attr_docs = self.analyzer.find_attr_docs()
            if self.objpath:
                key = ('.'.join(self.objpath[:-1]), self.objpath[-1])
                if key in attr_docs:
                    docstring = False
                    # make a copy of docstring for attributes to avoid cache
                    # the change of autodoc-process-docstring event.
                    attribute_docstrings = [list(attr_docs[key])]

                    for i, line in enumerate(self.process_doc(attribute_docstrings)):
                        self.add_line(line, sourcename, i)

        # add content from docstrings
        if docstring:
            docstrings = self.get_doc()
            if docstrings is None:
                # Do not call autodoc-process-docstring on get_doc() returns None.
                pass
            else:
                if not docstrings:
                    # append at least a dummy docstring, so that the event
                    # autodoc-process-docstring is fired and can add some
                    # content if desired
                    docstrings.append([])
                for i, line in enumerate(self.process_doc(docstrings)):
                    self.add_line(line, sourcename, i)

        # add additional content (e.g. from document), if present
        if more_content:
            for line, src in zip(more_content.data, more_content.items, strict=True):
                self.add_line(line, src[0], src[1])

    def get_object_members(self, want_all: bool) -> tuple[bool, list[ObjectMember]]:
        """Return `(members_check_module, members)` where `members` is a
        list of `(membername, member)` pairs of the members of *self.object*.

        If *want_all* is True, return all members.  Else, only return those
        members given by *self.options.members* (which may also be None).
        """
        msg = 'must be implemented in subclasses'
        raise NotImplementedError(msg)

    def filter_members(
        self, members: list[ObjectMember], want_all: bool
    ) -> list[tuple[str, Any, bool]]:
        """Filter the given member list.

        Members are skipped if

        - they are private (except if given explicitly or the private-members
          option is set)
        - they are special methods (except if given explicitly or the
          special-members option is set)
        - they are undocumented (except if the undoc-members option is set)

        The user can override the skipping decision by connecting to the
        ``autodoc-skip-member`` event.
        """

        def is_filtered_inherited_member(name: str, obj: Any) -> bool:
            inherited_members = self.options.inherited_members or set()
            seen = set()

            if inspect.isclass(self.object):
                for cls in self.object.__mro__:
                    if name in cls.__dict__:
                        seen.add(cls)
                    if (
                        cls.__name__ in inherited_members
                        and cls != self.object
                        and any(
                            issubclass(potential_child, cls) for potential_child in seen
                        )
                    ):
                        # given member is a member of specified *super class*
                        return True
                    if name in cls.__dict__:
                        return False
                    if name in self.get_attr(cls, '__annotations__', {}):
                        return False
                    if isinstance(obj, ObjectMember) and obj.class_ is cls:
                        return False

            return False

        ret = []

        # search for members in source code too
        namespace = '.'.join(self.objpath)  # will be empty for modules

        if self.analyzer:
            attr_docs = self.analyzer.find_attr_docs()
        else:
            attr_docs = {}

        # process members and determine which to skip
        for obj in members:
            membername = obj.__name__
            member = obj.object

            # ---------------------------------------------------
            # Issue #17455: Immediately skip lazy imports to avoid
            # deprecation messages.
            from sage.misc.lazy_import import LazyImport
            if isinstance(member, LazyImport):
                continue
            # ---------------------------------------------------

            # if isattr is True, the member is documented as an attribute
            isattr = member is INSTANCEATTR or (namespace, membername) in attr_docs

            try:
                doc = getdoc(
                    member,
                    self.get_attr,
                    self.config.autodoc_inherit_docstrings,
                    self.object,
                    membername,
                )
                if not isinstance(doc, str):
                    # Ignore non-string __doc__
                    doc = None

                # if the member __doc__ is the same as self's __doc__, it's just
                # inherited and therefore not the member's doc
                cls = self.get_attr(member, '__class__', None)
                if cls:
                    cls_doc = self.get_attr(cls, '__doc__', None)
                    if cls_doc == doc:
                        doc = None

                if isinstance(obj, ObjectMember) and obj.docstring:
                    # hack for ClassDocumenter to inject docstring via ObjectMember
                    doc = obj.docstring

                doc, metadata = separate_metadata(doc)
                has_doc = bool(doc)

                if 'private' in metadata:
                    # consider a member private if docstring has "private" metadata
                    isprivate = True
                elif 'public' in metadata:
                    # consider a member public if docstring has "public" metadata
                    isprivate = False
                else:
                    isprivate = membername.startswith('_')

                keep = False
                if ismock(member) and (namespace, membername) not in attr_docs:
                    # mocked module or object
                    pass
                elif (
                    self.options.exclude_members
                    and membername in self.options.exclude_members
                ):
                    # remove members given by exclude-members
                    keep = False
                elif want_all and special_member_re.match(membername):
                    # special __methods__
                    if (
                        self.options.special_members
                        and membername in self.options.special_members
                    ):
                        if membername == '__doc__':  # NoQA: SIM114
                            keep = False
                        elif is_filtered_inherited_member(membername, obj):
                            keep = False
                        else:
                            keep = has_doc or self.options.undoc_members
                    else:
                        keep = False
                elif (namespace, membername) in attr_docs:
                    if want_all and isprivate:
                        if self.options.private_members is None:
                            keep = False
                        else:
                            keep = membername in self.options.private_members
                    else:
                        # keep documented attributes
                        keep = True
                elif want_all and isprivate:
                    if has_doc or self.options.undoc_members:
                        if self.options.private_members is None:  # NoQA: SIM114
                            keep = False
                        elif is_filtered_inherited_member(membername, obj):
                            keep = False
                        else:
                            keep = membername in self.options.private_members
                    else:
                        keep = False
                elif self.options.members is ALL and is_filtered_inherited_member(
                    membername, obj
                ):
                    keep = False
                else:
                    # ignore undocumented members if :undoc-members: is not given
                    keep = has_doc or self.options.undoc_members

                if isinstance(obj, ObjectMember) and obj.skipped:
                    # forcedly skipped member (ex. a module attribute not defined in __all__)
                    keep = False

                # give the user a chance to decide whether this member
                # should be skipped
                if self._events is not None:
                    # let extensions preprocess docstrings
                    skip_user = self._events.emit_firstresult(
                        'autodoc-skip-member',
                        self.objtype,
                        membername,
                        member,
                        not keep,
                        self.options,
                    )
                    if skip_user is not None:
                        keep = not skip_user
            except Exception as exc:
                logger.warning(
                    __(
                        'autodoc: failed to determine %s.%s (%r) to be documented, '
                        'the following exception was raised:\n%s'
                    ),
                    self.name,
                    membername,
                    member,
                    exc,
                    type='autodoc',
                )
                keep = False

            if keep:
                ret.append((membername, member, isattr))

        return ret

    def document_members(self, all_members: bool = False) -> None:
        """Generate reST for member documentation.

        If *all_members* is True, document all members, else those given by
        *self.options.members*.
        """
        # set current namespace for finding members
        self._current_document.autodoc_module = self.modname
        if self.objpath:
            self._current_document.autodoc_class = self.objpath[0]

        want_all = (
            all_members or self.options.inherited_members or self.options.members is ALL
        )
        # find out which members are documentable
        members_check_module, members = self.get_object_members(want_all)

        # document non-skipped members
        member_documenters: list[tuple[Documenter, bool]] = []
        for mname, member, isattr in self.filter_members(members, want_all):
            classes = [
                cls
                for cls in self.documenters.values()
                if cls.can_document_member(member, mname, isattr, self)
            ]
            if not classes:
                # don't know how to document this member
                continue
            # prefer the documenter with the highest priority
            classes.sort(key=lambda cls: cls.priority)
            # give explicitly separated module name, so that members
            # of inner classes can be documented
            full_mname = f'{self.modname}::' + '.'.join((*self.objpath, mname))
            documenter = classes[-1](self.directive, full_mname, self.indent)
            member_documenters.append((documenter, isattr))

        member_order = self.options.member_order or self.config.autodoc_member_order
        # We now try to import all objects before ordering them. This is to
        # avoid possible circular imports if we were to import objects after
        # their associated documenters have been sorted.
        member_documenters = [
            (documenter, isattr)
            for documenter, isattr in member_documenters
            if documenter.parse_name() and documenter.import_object()
        ]
        member_documenters = self.sort_members(member_documenters, member_order)

        for documenter, isattr in member_documenters:
            assert documenter.modname
            # We can directly call ._generate() since the documenters
            # already called parse_name() and import_object() before.
            #
            # Note that those two methods above do not emit events, so
            # whatever objects we deduced should not have changed.
            documenter._generate(
                all_members=True,
                real_modname=self.real_modname,
                check_module=members_check_module and not isattr,
            )

        # reset current objects
        self._current_document.autodoc_module = ''
        self._current_document.autodoc_class = ''

    def sort_members(
        self, documenters: list[tuple[Documenter, bool]], order: str
    ) -> list[tuple[Documenter, bool]]:
        """Sort the given member list."""
        if order == 'groupwise':
            # sort by group; alphabetically within groups
            documenters.sort(key=lambda e: (e[0].member_order, e[0].name))
        elif order == 'bysource':
            # By default, member discovery order matches source order,
            # as dicts are insertion-ordered.
            if self.analyzer:
                # sort by source order, by virtue of the module analyzer
                tagorder = self.analyzer.tagorder

                def keyfunc(entry: tuple[Documenter, bool]) -> int:
                    fullname = entry[0].name.split('::')[1]
                    return tagorder.get(fullname, len(tagorder))

                documenters.sort(key=keyfunc)
        else:  # alphabetical
            documenters.sort(key=lambda e: e[0].name)

        return documenters

    def generate(
        self,
        more_content: StringList | None = None,
        real_modname: str | None = None,
        check_module: bool = False,
        all_members: bool = False,
    ) -> None:
        """Generate reST for the object given by *self.name*, and possibly for
        its members.

        If *more_content* is given, include that content. If *real_modname* is
        given, use that module name to find attribute docs. If *check_module* is
        True, only generate if the object is defined in the module name it is
        imported from. If *all_members* is True, document all members.
        """
        if not self.parse_name():
            # need a module to import
            logger.warning(
                __(
                    "don't know which module to import for autodocumenting "
                    '%r (try placing a "module" or "currentmodule" directive '
                    'in the document, or giving an explicit module name)'
                ),
                self.name,
                type='autodoc',
            )
            return

        # now, import the module and get object to document
        if not self.import_object():
            return

        self._generate(more_content, real_modname, check_module, all_members)

    def _generate(
        self,
        more_content: StringList | None = None,
        real_modname: str | None = None,
        check_module: bool = False,
        all_members: bool = False,
    ) -> None:
        # If there is no real module defined, figure out which to use.
        # The real module is used in the module analyzer to look up the module
        # where the attribute documentation would actually be found in.
        # This is used for situations where you have a module that collects the
        # functions and classes of internal submodules.
        guess_modname = self.get_real_modname()
        self.real_modname: str = real_modname or guess_modname

        # try to also get a source code analyzer for attribute docs
        try:
            self.analyzer = ModuleAnalyzer.for_module(self.real_modname)
            # parse right now, to get PycodeErrors on parsing (results will
            # be cached anyway)
            self.analyzer.find_attr_docs()
        except PycodeError as exc:
            logger.debug('[autodoc] module analyzer failed: %s', exc)
            # no source file -- e.g. for builtin and C modules
            self.analyzer = None
            # at least add the module.__file__ as a dependency
            if module___file__ := getattr(self.module, '__file__', ''):
                self.directive.record_dependencies.add(module___file__)
        else:
            self.directive.record_dependencies.add(self.analyzer.srcname)

        if self.real_modname != guess_modname:
            # Add module to dependency list if target object is defined in other module.
            try:
                analyzer = ModuleAnalyzer.for_module(guess_modname)
                self.directive.record_dependencies.add(analyzer.srcname)
            except PycodeError:
                pass

        docstrings: list[str] = functools.reduce(
            operator.iadd, self.get_doc() or [], []
        )
        if ismock(self.object) and not docstrings:
            logger.warning(
                __('A mocked object is detected: %r'),
                self.name,
                type='autodoc',
                subtype='mocked_object',
            )

        # check __module__ of object (for members not given explicitly)
        if check_module:
            if not self.check_module():
                return

        sourcename = self.get_sourcename()

        # make sure that the result starts with an empty line.  This is
        # necessary for some situations where another directive preprocesses
        # reST and no starting newline is present
        self.add_line('', sourcename)

        # format the object's signature, if any
        try:
            sig = self.format_signature()
        except Exception as exc:
            logger.warning(
                __('error while formatting signature for %s: %s'),
                self.fullname,
                exc,
                type='autodoc',
            )
            return

        # generate the directive header and options, if applicable
        self.add_directive_header(sig)
        self.add_line('', sourcename)

        # e.g. the module directive doesn't have content
        self.indent += self.content_indent

        # add all content (from docstrings, attribute docs etc.)
        self.add_content(more_content)

        # document members, if possible
        self.document_members(all_members)


class ModuleDocumenter(Documenter):
    """Specialized Documenter subclass for modules."""

    objtype = 'module'
    content_indent = ''
    _extra_indent = '   '

    option_spec: ClassVar[dict[str, Any]] = {
        'members': members_option,
        'undoc-members': bool_option,
        'no-index': bool_option,
        'no-index-entry': bool_option,
        'inherited-members': inherited_members_option,
        'show-inheritance': bool_option,
        'synopsis': identity,
        'platform': identity,
        'deprecated': bool_option,
        'member-order': member_order_option,
        'exclude-members': exclude_members_option,
        'private-members': members_option,
        'special-members': members_option,
        'imported-members': bool_option,
        'ignore-module-all': bool_option,
        'no-value': bool_option,
        'noindex': bool_option,
    }

    def __init__(self, *args: Any) -> None:
        super().__init__(*args)
        merge_members_option(self.options)
        self.__all__: Sequence[str] | None = None

    def add_content(self, more_content: StringList | None) -> None:
        old_indent = self.indent
        self.indent += self._extra_indent
        super().add_content(None)
        self.indent = old_indent
        if more_content:
            for line, src in zip(more_content.data, more_content.items, strict=True):
                self.add_line(line, src[0], src[1])

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        # don't document submodules automatically
        return False

    def resolve_name(
        self, modname: str | None, parents: Any, path: str, base: str
    ) -> tuple[str | None, list[str]]:
        if modname is not None:
            logger.warning(
                __('"::" in automodule name doesn\'t make sense'), type='autodoc'
            )
        return (path or '') + base, []

    def parse_name(self) -> bool:
        ret = super().parse_name()
        if self.args or self.retann:
            logger.warning(
                __('signature arguments or return annotation given for automodule %s'),
                self.fullname,
                type='autodoc',
            )
        return ret

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)

        try:
            if not self.options.ignore_module_all:
                self.__all__ = inspect.getall(self.object)
        except ValueError as exc:
            # invalid __all__ found.
            logger.warning(
                __(
                    '__all__ should be a list of strings, not %r '
                    '(in module %s) -- ignoring __all__'
                ),
                exc.args[0],
                self.fullname,
                type='autodoc',
            )

        return ret

    def add_directive_header(self, sig: str) -> None:
        Documenter.add_directive_header(self, sig)

        sourcename = self.get_sourcename()

        # add some module-specific options
        if self.options.synopsis:
            self.add_line('   :synopsis: ' + self.options.synopsis, sourcename)
        if self.options.platform:
            self.add_line('   :platform: ' + self.options.platform, sourcename)
        if self.options.deprecated:
            self.add_line('   :deprecated:', sourcename)
        if self.options.no_index_entry:
            self.add_line('   :no-index-entry:', sourcename)

    def get_module_members(self) -> dict[str, ObjectMember]:
        """Get members of target module."""
        if self.analyzer:
            attr_docs = self.analyzer.attr_docs
        else:
            attr_docs = {}

        members: dict[str, ObjectMember] = {}
        for name in dir(self.object):
            try:
                value = safe_getattr(self.object, name, None)
                if ismock(value):
                    value = undecorate(value)
                docstring = attr_docs.get(('', name), [])
                members[name] = ObjectMember(
                    name, value, docstring='\n'.join(docstring)
                )
            except AttributeError:
                continue

        # annotation only member (ex. attr: int)
        for name in inspect.getannotations(self.object):
            if name not in members:
                docstring = attr_docs.get(('', name), [])
                members[name] = ObjectMember(
                    name, INSTANCEATTR, docstring='\n'.join(docstring)
                )

        return members

    def get_object_members(self, want_all: bool) -> tuple[bool, list[ObjectMember]]:
        members = self.get_module_members()
        if want_all:
            if self.__all__ is None:
                # for implicit module members, check __module__ to avoid
                # documenting imported objects
                return True, list(members.values())
            for member in members.values():
                if member.__name__ not in self.__all__:
                    member.skipped = True
            return False, list(members.values())
        memberlist = self.options.members or []
        ret = []
        for name in memberlist:
            if name in members:
                ret.append(members[name])
            else:
                logger.warning(
                    __(
                        'missing attribute mentioned in :members: option: '
                        'module %s, attribute %s'
                    ),
                    safe_getattr(self.object, '__name__', '???'),
                    name,
                    type='autodoc',
                )
        return False, ret

    def sort_members(
        self, documenters: list[tuple[Documenter, bool]], order: str
    ) -> list[tuple[Documenter, bool]]:
        if order == 'bysource' and self.__all__:
            assert self.__all__ is not None
            module_all = self.__all__
            module_all_set = set(module_all)
            module_all_len = len(module_all)

            # Sort alphabetically first (for members not listed on the __all__)
            documenters.sort(key=lambda e: e[0].name)

            # Sort by __all__
            def keyfunc(entry: tuple[Documenter, bool]) -> int:
                name = entry[0].name.split('::')[1]
                if name in module_all_set:
                    return module_all.index(name)
                return module_all_len

            documenters.sort(key=keyfunc)

            return documenters
        return super().sort_members(documenters, order)


class ModuleLevelDocumenter(Documenter):
    """Specialized Documenter subclass for objects on module level (functions,
    classes, data/constants).
    """

    def resolve_name(
        self, modname: str | None, parents: Any, path: str, base: str
    ) -> tuple[str | None, list[str]]:
        if modname is not None:
            return modname, [*parents, base]
        if path:
            modname = path.rstrip('.')
            return modname, [*parents, base]

        # if documenting a toplevel object without explicit module,
        # it can be contained in another auto directive ...
        modname = self._current_document.autodoc_module
        # ... or in the scope of a module directive
        if not modname:
            modname = self.env.ref_context.get('py:module')
        # ... else, it stays None, which means invalid
        return modname, [*parents, base]


class ClassLevelDocumenter(Documenter):
    """Specialized Documenter subclass for objects on class level (methods,
    attributes).
    """

    def resolve_name(
        self, modname: str | None, parents: Any, path: str, base: str
    ) -> tuple[str | None, list[str]]:
        if modname is not None:
            return modname, [*parents, base]

        if path:
            mod_cls = path.rstrip('.')
        else:
            # if documenting a class-level object without path,
            # there must be a current class, either from a parent
            # auto directive ...
            mod_cls = self._current_document.autodoc_class
            # ... or from a class directive
            if not mod_cls:
                mod_cls = self.env.ref_context.get('py:class', '')
                # ... if still falsy, there's no way to know
                if not mod_cls:
                    return None, []
        modname, sep, cls = mod_cls.rpartition('.')
        parents = [cls]
        # if the module name is still missing, get it like above
        if not modname:
            modname = self._current_document.autodoc_module
        if not modname:
            modname = self.env.ref_context.get('py:module')
        # ... else, it stays None, which means invalid
        return modname, [*parents, base]


class DocstringSignatureMixin:
    """Mixin for FunctionDocumenter and MethodDocumenter to provide the
    feature of reading the signature from the docstring.
    """

    _new_docstrings: list[list[str]] | None = None
    _signatures: list[str] = []

    def _find_signature(self) -> tuple[str | None, str | None] | None:
        # candidates of the object name
        valid_names = [self.objpath[-1]]  # type: ignore[attr-defined]
        if isinstance(self, ClassDocumenter):
            valid_names.append('__init__')
            if hasattr(self.object, '__mro__'):
                valid_names.extend(cls.__name__ for cls in self.object.__mro__)

        docstrings = self.get_doc()
        if docstrings is None:
            return None, None
        self._new_docstrings = docstrings[:]
        self._signatures = []
        result = None
        for i, doclines in enumerate(docstrings):
            for j, line in enumerate(doclines):
                if not line:
                    # no lines in docstring, no match
                    break

                if line.endswith('\\'):
                    line = line.rstrip('\\').rstrip()

                # match first line of docstring against signature RE
                match = py_ext_sig_re.match(line)
                if not match:
                    break
                exmod, path, base, tp_list, args, retann = match.groups()

                # the base name must match ours
                if base not in valid_names:
                    break

                # re-prepare docstring to ignore more leading indentation
                directive = self.directive  # type: ignore[attr-defined]
                tab_width = directive.state.document.settings.tab_width
                self._new_docstrings[i] = prepare_docstring(
                    '\n'.join(doclines[j + 1 :]), tab_width
                )

                if result is None:
                    # first signature
                    result = args, retann
                else:
                    # subsequent signatures
                    self._signatures.append(f'({args}) -> {retann}')

            if result is not None:
                # finish the loop when signature found
                break

        return result

    def get_doc(self) -> list[list[str]] | None:
        if self._new_docstrings is not None:
            return self._new_docstrings
        return super().get_doc()  # type: ignore[misc]

    def format_signature(self, **kwargs: Any) -> str:
        self.args: str | None
        if self.args is None and self.config.autodoc_docstring_signature:  # type: ignore[attr-defined]
            # only act if a signature is not explicitly given already, and if
            # the feature is enabled
            result = self._find_signature()
            if result is not None:
                self.args, self.retann = result
        sig = super().format_signature(**kwargs)  # type: ignore[misc]
        if self._signatures:
            return '\n'.join((sig, *self._signatures))
        return sig


class DocstringStripSignatureMixin(DocstringSignatureMixin):
    """Mixin for AttributeDocumenter to provide the
    feature of stripping any function signature from the docstring.
    """

    def format_signature(self, **kwargs: Any) -> str:
        if self.args is None and self.config.autodoc_docstring_signature:  # type: ignore[attr-defined]
            # only act if a signature is not explicitly given already, and if
            # the feature is enabled
            result = self._find_signature()
            if result is not None:
                # Discarding _args is a only difference with
                # DocstringSignatureMixin.format_signature.
                # Documenter.format_signature use self.args value to format.
                _args, self.retann = result
        return super().format_signature(**kwargs)


class FunctionDocumenter(DocstringSignatureMixin, ModuleLevelDocumenter):  # type: ignore[misc]
    """Specialized Documenter subclass for functions."""

    objtype = 'function'
    member_order = 30

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        # --------------------------------------------------------------------
        # supports functions, builtins but, unlike Sphinx' autodoc,
        # does not support bound methods exported at the module level
        if is_function_or_cython_function(member) or inspect.isbuiltin(member):
            return True
        # Issue #9976: It can be documented if it is a genuine function.
        # Often, a class instance has the same documentation as its class,
        # and then we typically want to document the class and not the
        # instance.  However, there is an exception: CachedFunction(f) returns
        # a class instance, whose doc string coincides with that of f and is
        # thus different from that of the class CachedFunction. In that
        # situation, we want that f is documented.
        return (isclassinstance(member) and
                sage_getdoc_original(member) != sage_getdoc_original(member.__class__))
        # --------------------------------------------------------------------

    def format_args(self, **kwargs: Any) -> str:
        if self.config.autodoc_typehints in {'none', 'description'}:
            kwargs.setdefault('show_annotation', False)
        if self.config.autodoc_typehints_format == 'short':
            kwargs.setdefault('unqualified_typehints', True)
        if self.config.python_display_short_literal_types:
            kwargs.setdefault('short_literals', True)

        try:
            self._events.emit('autodoc-before-process-signature', self.object, False)
            # ----------------------------------------------------------------
            # Issue #9976: Support the _sage_argspec_ attribute which makes it
            # possible to get argument specification of decorated callables in
            # documentation correct. See e.g. sage.misc.decorators.sage_wraps
            obj = self.object

            if hasattr(obj, "_sage_argspec_"):
                argspec = obj._sage_argspec_()
            if inspect.isbuiltin(obj) or \
                   inspect.ismethoddescriptor(obj):
                # cannot introspect arguments of a C function or method
                # unless a function to do so is supplied
                argspec = sage_getargspec(obj)
            argspec = sage_getargspec(obj)
            if isclassinstance(obj) or inspect.isclass(obj):
                # if a class should be documented as function, we try
                # to use the constructor signature as function
                # signature without the first argument.
                if argspec is not None and argspec[0]:
                    del argspec[0][0]

            if argspec is None:
                return None
            args = sage_formatargspec(*argspec)
            # ----------------------------------------------------------------
        except TypeError as exc:
            logger.warning(
                __('Failed to get a function signature for %s: %s'), self.fullname, exc
            )
            return ''
        except ValueError:
            args = ''

        if self.config.strip_signature_backslash:
            # escape backslashes for reST
            args = args.replace('\\', '\\\\')
        return args

    def document_members(self, all_members: bool = False) -> None:
        pass

    def add_directive_header(self, sig: str) -> None:
        sourcename = self.get_sourcename()
        super().add_directive_header(sig)

        is_coro = inspect.iscoroutinefunction(self.object)
        is_acoro = inspect.isasyncgenfunction(self.object)
        if is_coro or is_acoro:
            self.add_line('   :async:', sourcename)

    def format_signature(self, **kwargs: Any) -> str:
        if self.config.autodoc_typehints_format == 'short':
            kwargs.setdefault('unqualified_typehints', True)
        if self.config.python_display_short_literal_types:
            kwargs.setdefault('short_literals', True)

        sigs = []
        if (
            self.analyzer
            and '.'.join(self.objpath) in self.analyzer.overloads
            and self.config.autodoc_typehints != 'none'
        ):
            # Use signatures for overloaded functions instead of the implementation function.
            overloaded = True
        else:
            overloaded = False
            sig = super().format_signature(**kwargs)
            sigs.append(sig)

        if inspect.is_singledispatch_function(self.object):
            # append signature of singledispatch'ed functions
            for typ, func in self.object.registry.items():
                if typ is object:
                    pass  # default implementation. skipped.
                else:
                    dispatchfunc = self.annotate_to_first_argument(func, typ)
                    if dispatchfunc:
                        documenter = FunctionDocumenter(self.directive, '')
                        documenter.object = dispatchfunc
                        documenter.objpath = ['']
                        sigs.append(documenter.format_signature())
        if overloaded and self.analyzer is not None:
            actual = inspect.signature(
                self.object, type_aliases=self._type_aliases()
            )
            __globals__ = safe_getattr(self.object, '__globals__', {})
            for overload in self.analyzer.overloads['.'.join(self.objpath)]:
                overload = self.merge_default_value(actual, overload)
                overload = evaluate_signature(
                    overload, __globals__, self._type_alias_namespace()
                )

                sig = stringify_signature(overload, **kwargs)
                sigs.append(sig)

        return '\n'.join(sigs)

    def merge_default_value(self, actual: Signature, overload: Signature) -> Signature:
        """Merge default values of actual implementation to the overload variants."""
        parameters = list(overload.parameters.values())
        for i, param in enumerate(parameters):
            actual_param = actual.parameters.get(param.name)
            if actual_param and param.default == '...':
                parameters[i] = param.replace(default=actual_param.default)

        return overload.replace(parameters=parameters)

    def annotate_to_first_argument(
        self, func: Callable[..., Any], typ: type
    ) -> Callable[..., Any] | None:
        """Annotate type hint to the first argument of function if needed."""
        try:
            sig = inspect.signature(func, type_aliases=self._type_aliases(func))
        except TypeError as exc:
            logger.warning(
                __('Failed to get a function signature for %s: %s'), self.fullname, exc
            )
            return None
        except ValueError:
            return None

        if len(sig.parameters) == 0:
            return None

        def dummy():  # type: ignore[no-untyped-def]  # NoQA: ANN202
            pass

        params = list(sig.parameters.values())
        if params[0].annotation is Parameter.empty:
            params[0] = params[0].replace(annotation=typ)
            try:
                dummy.__signature__ = sig.replace(parameters=params)  # type: ignore[attr-defined]
                return dummy
            except (AttributeError, TypeError):
                # failed to update signature (ex. built-in or extension types)
                return None

        return func


class DecoratorDocumenter(FunctionDocumenter):
    """Specialized Documenter subclass for decorator functions."""

    objtype = 'decorator'

    # must be lower than FunctionDocumenter
    priority = -1

    def format_args(self, **kwargs: Any) -> str:
        args = super().format_args(**kwargs)
        if ',' in args:
            return args
        return ''


# Types which have confusing metaclass signatures it would be best not to show.
# These are listed by name, rather than storing the objects themselves, to avoid
# needing to import the modules.
_METACLASS_CALL_BLACKLIST = frozenset({
    'enum.EnumType.__call__',
})


# Types whose __new__ signature is a pass-through.
_CLASS_NEW_BLACKLIST = frozenset({
    'typing.Generic.__new__',
})


class ClassDocumenter(DocstringSignatureMixin, ModuleLevelDocumenter):  # type: ignore[misc]
    """Specialized Documenter subclass for classes."""

    objtype = 'class'
    member_order = 20
    option_spec: ClassVar[dict[str, Any]] = {
        'members': members_option,
        'undoc-members': bool_option,
        'no-index': bool_option,
        'no-index-entry': bool_option,
        'inherited-members': inherited_members_option,
        'show-inheritance': bool_option,
        'member-order': member_order_option,
        'exclude-members': exclude_members_option,
        'private-members': members_option,
        'special-members': members_option,
        'class-doc-from': class_doc_from_option,
        'noindex': bool_option,
    }

    # Must be higher than FunctionDocumenter, ClassDocumenter, and
    # AttributeDocumenter as NewType can be an attribute and is a class
    # after Python 3.10.
    priority = 15

    _signature_class: Any = None
    _signature_method_name: str = ''

    def __init__(self, *args: Any) -> None:
        super().__init__(*args)

        if self.config.autodoc_class_signature == 'separated':
            self.options = self.options.copy()

            # show __init__() method
            if self.options.special_members is None:
                self.options['special-members'] = ['__new__', '__init__']
            else:
                self.options.special_members.append('__new__')
                self.options.special_members.append('__init__')

        merge_members_option(self.options)

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        return isinstance(member, type) or (
            isattr and isinstance(member, NewType | TypeVar)
        )

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)
        # if the class is documented under another name, document it
        # as data/attribute
        if ret:
            if hasattr(self.object, '__name__'):
                self.doc_as_attr = self.objpath[-1] != self.object.__name__
            # -------------------------------------------------------------------
            # Issue #27692, #7448: The original goal of this was that if some
            # class is aliased, the alias is generated as a link rather than
            # duplicated. For example in
            #
            #     class A:
            #         pass
            #     B = A
            #
            # Then B is an alias of A, and should be generated as such.
            #
            # The way it was solved is to compare the name under which the
            # current class is found and the actual name of the class (stored in
            # the attribute __name__):
            #
            #   if hasattr(self.object, '__name__'):
            #       self.doc_as_attr = (self.objpath[-1] != self.object.__name__)
            #   else:
            #       self.doc_as_attr = True
            #
            # Now, to work around a pickling bug of nested class in Python,
            # by using the metaclass NestedMetaclass, we change the attribute
            # __name__ of the nested class. For example, in
            #
            #     class A(metaclass=NestedMetaclass):
            #        class B():
            #            pass
            #
            # the class B get its name changed to 'A.B'. Such dots '.' in names
            # are not supposed to occur in normal python name. I use it to check
            # if the class is a nested one and to compare its __name__ with its
            # path.
            #
            # The original implementation as well as the new one here doesn't work
            # if a class is aliased from a different place under the same name. For
            # example, in the following,
            #
            #     class A:
            #         pass
            #     class Container:
            #         A = A
            #
            # The nested copy Container.A is also documented. Actually, it seems
            # that there is no way to solve this by introspection. I'll submbit
            # this problem on sphinx trac.
            #
            # References: trac #5986, file sage/misc/nested_class.py
            import sys
            module = getattr(self.object, '__module__', False)
            name = getattr(self.object, '__name__', False)
            qualname = getattr(self.object, '__qualname__', name)
            if qualname and module:
                # walk the standard attribute lookup path for this object
                qualname_parts = qualname.split('.')
                cls = getattr(sys.modules[module], qualname_parts[0], None)
                for part in qualname_parts[1:]:
                    if cls is None:
                        break
                    cls = getattr(cls, part, None)
                self.doc_as_attr = (self.objpath != qualname_parts and
                                    self.object is cls)
            # -------------------------------------------------------------------
            else:
                self.doc_as_attr = True
            if isinstance(self.object, NewType | TypeVar):
                modname = getattr(self.object, '__module__', self.modname)
                if modname != self.modname and self.modname.startswith(modname):
                    bases = self.modname[len(modname) :].strip('.').split('.')
                    self.objpath = bases + self.objpath
                    self.modname = modname
        return ret

    def _get_signature(self) -> tuple[Any | None, str | None, Signature | None]:
        if isinstance(self.object, NewType | TypeVar):
            # Suppress signature
            return None, None, None

        def get_user_defined_function_or_method(obj: Any, attr: str) -> Any:
            """Get the `attr` function or method from `obj`, if it is user-defined."""
            if inspect.is_builtin_class_method(obj, attr):
                return None
            attr = self.get_attr(obj, attr, None)
            if not (inspect.ismethod(attr) or inspect.isfunction(attr)):
                return None
            return attr

        # This sequence is copied from inspect._signature_from_callable.
        # ValueError means that no signature could be found, so we keep going.

        # First, we check if obj has a __signature__ attribute
        if hasattr(self.object, '__signature__'):
            object_sig = self.object.__signature__
            if isinstance(object_sig, Signature):
                return None, None, object_sig
            if sys.version_info[:2] <= (3, 14) and callable(object_sig):
                # Support for enum.Enum.__signature__ in Python 3.12, 3.13, and 3.14.
                if isinstance(object_sig_str := object_sig(), str):
                    return None, None, inspect.signature_from_str(object_sig_str)

        # Next, let's see if it has an overloaded __call__ defined
        # in its metaclass
        call = get_user_defined_function_or_method(type(self.object), '__call__')

        if call is not None:
            if f'{call.__module__}.{call.__qualname__}' in _METACLASS_CALL_BLACKLIST:
                call = None

        if call is not None:
            self._events.emit('autodoc-before-process-signature', call, True)
            try:
                sig = inspect.signature(
                    call,
                    bound_method=True,
                    type_aliases=self._type_aliases(call),
                )
                return type(self.object), '__call__', sig
            except ValueError:
                pass

        # Now we check if the 'obj' class has a '__new__' method
        new = get_user_defined_function_or_method(self.object, '__new__')

        if new is not None:
            if f'{new.__module__}.{new.__qualname__}' in _CLASS_NEW_BLACKLIST:
                new = None

        if new is not None:
            self._events.emit('autodoc-before-process-signature', new, True)
            try:
                sig = inspect.signature(
                    new,
                    bound_method=True,
                    type_aliases=self._type_aliases(new),
                )
                return self.object, '__new__', sig
            except ValueError:
                pass

        # Finally, we should have at least __init__ implemented
        init = get_user_defined_function_or_method(self.object, '__init__')
        if init is not None:
            self._events.emit('autodoc-before-process-signature', init, True)
            try:
                sig = inspect.signature(
                    init,
                    bound_method=True,
                    type_aliases=self._type_aliases(init),
                )
                return self.object, '__init__', sig
            except ValueError:
                pass

        # None of the attributes are user-defined, so fall back to let inspect
        # handle it.
        # We don't know the exact method that inspect.signature will read
        # the signature from, so just pass the object itself to our hook.
        self._events.emit('autodoc-before-process-signature', self.object, False)
        try:
            sig = inspect.signature(
                self.object,
                bound_method=False,
                type_aliases=self._type_aliases(),
            )
            return None, None, sig
        except ValueError:
            pass

        # Still no signature: happens e.g. for old-style classes
        # with __init__ in C and no `__text_signature__`.
        return None, None, None

    def format_args(self, **kwargs: Any) -> str:
        if self.config.autodoc_typehints in {'none', 'description'}:
            kwargs.setdefault('show_annotation', False)
        if self.config.autodoc_typehints_format == 'short':
            kwargs.setdefault('unqualified_typehints', True)
        if self.config.python_display_short_literal_types:
            kwargs.setdefault('short_literals', True)

        try:
            self._signature_class, _signature_method_name, sig = self._get_signature()
        except TypeError as exc:
            # __signature__ attribute contained junk
            logger.warning(
                __('Failed to get a constructor signature for %s: %s'),
                self.fullname,
                exc,
            )
            return ''
        self._signature_method_name = _signature_method_name or ''

        if sig is None:
            return ''

        return stringify_signature(sig, show_return_annotation=False, **kwargs)

    def _find_signature(self) -> tuple[str | None, str | None] | None:
        result = super()._find_signature()
        if result is not None:
            # Strip a return value from signature of constructor in docstring (first entry)
            result = (result[0], None)

        for i, sig in enumerate(self._signatures):
            if sig.endswith(' -> None'):
                # Strip a return value from signatures of constructor in docstring (subsequent
                # entries)
                self._signatures[i] = sig[:-8]

        return result

    def format_signature(self, **kwargs: Any) -> str:
        if self.doc_as_attr:
            return ''
        if self.config.autodoc_class_signature == 'separated':
            # do not show signatures
            return ''

        if self.config.autodoc_typehints_format == 'short':
            kwargs.setdefault('unqualified_typehints', True)
        if self.config.python_display_short_literal_types:
            kwargs.setdefault('short_literals', True)

        sig = super().format_signature()
        sigs = []

        overloads = self.get_overloaded_signatures()
        if overloads and self.config.autodoc_typehints != 'none':
            # Use signatures for overloaded methods instead of the implementation method.
            method = safe_getattr(
                self._signature_class, self._signature_method_name, None
            )
            __globals__ = safe_getattr(method, '__globals__', {})
            for overload in overloads:
                overload = evaluate_signature(
                    overload, __globals__, self._type_alias_namespace(method)
                )

                parameters = list(overload.parameters.values())
                overload = overload.replace(
                    parameters=parameters[1:], return_annotation=Parameter.empty
                )
                sig = stringify_signature(overload, **kwargs)
                sigs.append(sig)
        else:
            sigs.append(sig)

        return '\n'.join(sigs)

    def get_overloaded_signatures(self) -> list[Signature]:
        if self._signature_class and self._signature_method_name:
            for cls in self._signature_class.__mro__:
                try:
                    analyzer = ModuleAnalyzer.for_module(cls.__module__)
                    analyzer.analyze()
                    qualname = f'{cls.__qualname__}.{self._signature_method_name}'
                    if qualname in analyzer.overloads:
                        return analyzer.overloads.get(qualname, [])
                    if qualname in analyzer.tagorder:
                        # the constructor is defined in the class, but not overridden.
                        return []
                except PycodeError:
                    pass

        return []

    def get_canonical_fullname(self) -> str | None:
        __modname__ = safe_getattr(self.object, '__module__', self.modname)
        __qualname__ = safe_getattr(self.object, '__qualname__', None)
        if __qualname__ is None:
            __qualname__ = safe_getattr(self.object, '__name__', None)
        if __qualname__ and '<locals>' in __qualname__:
            # No valid qualname found if the object is defined as locals
            __qualname__ = None

        if __modname__ and __qualname__:
            return f'{__modname__}.{__qualname__}'
        return None

    def add_directive_header(self, sig: str) -> None:
        sourcename = self.get_sourcename()

        if self.doc_as_attr:
            self.directivetype = 'attribute'
        super().add_directive_header(sig)

        if isinstance(self.object, NewType | TypeVar):
            return

        if self.analyzer and '.'.join(self.objpath) in self.analyzer.finals:
            self.add_line('   :final:', sourcename)

        canonical_fullname = self.get_canonical_fullname()
        if (
            not self.doc_as_attr
            and not isinstance(self.object, NewType)
            and canonical_fullname
            and self.fullname != canonical_fullname
        ):
            self.add_line('   :canonical: %s' % canonical_fullname, sourcename)

        # add inheritance info, if wanted
        if not self.doc_as_attr and self.options.show_inheritance:
            if inspect.getorigbases(self.object):
                # A subclass of generic types
                # refs: PEP-560 <https://peps.python.org/pep-0560/>
                bases = list(self.object.__orig_bases__)
            elif hasattr(self.object, '__bases__') and len(self.object.__bases__):
                # A normal class
                bases = list(self.object.__bases__)
            else:
                bases = []

            self._events.emit(
                'autodoc-process-bases', self.fullname, self.object, self.options, bases
            )

            mode = _get_render_mode(self.config.autodoc_typehints_format)
            base_classes = [self._restify_base_class(cls, mode=mode) for cls in bases]

            sourcename = self.get_sourcename()
            self.add_line('', sourcename)
            self.add_line('   ' + _('Bases: %s') % ', '.join(base_classes), sourcename)

    def _restify_base_class(self, cls: Any, mode: str) -> str:
        """
        Return the reStructuredText representation of a base class.

        Python's inventory documents ``functools.partial`` as a function. Use
        that role for the inheritance display as well.

        Sage dynamically creates category helper classes such as
        ``FiniteEnumeratedSets.parent_class``.  These classes are not
        importable Python-domain objects, but the underlying category class is.
        Link the inheritance display to that stable category target.

        TESTS::

            sage: import functools
            sage: from sage_docbuild.ext.sage_autodoc import ClassDocumenter
            sage: ClassDocumenter._restify_base_class(None, functools.partial, 'smart')
            ':py:func:`~functools.partial`'
            sage: ClassDocumenter._restify_base_class(
            ....:     None, functools.partial, 'fully-qualified-except-typing')
            ':py:func:`functools.partial`'
        """
        if cls is functools.partial:
            # Sphinx identifies this callable type as a class, but Python's
            # inventory documents it as a function.  Keep the established
            # restify target/title (including smart mode's ``~``) and change
            # only the role that is wrong for this one inheritance display.
            return restify(cls, mode=mode).replace(':py:class:', ':py:func:', 1)
        if isinstance(cls, type):
            module = getattr(cls, '__module__', '')
            qualname = getattr(cls, '__qualname__', '')
            helper_names = ('parent_class', 'element_class', 'subcategory_class')
            category_name, sep, helper_name = qualname.rpartition('.')
            if (
                module.startswith('sage.categories.')
                and sep
                and helper_name in helper_names
            ):
                target = f'{module}.{category_name}'
                return f':class:`{qualname} <{target}>`'
        return restify(cls, mode=mode)

    def get_object_members(self, want_all: bool) -> tuple[bool, list[ObjectMember]]:
        members = get_class_members(
            self.object,
            self.objpath,
            self.get_attr,
            self.config.autodoc_inherit_docstrings,
        )
        if not want_all:
            if not self.options.members:
                return False, []
            # specific members given
            selected = []
            for name in self.options.members:
                if name in members:
                    selected.append(members[name])
                else:
                    logger.warning(
                        __('missing attribute %s in object %s'),
                        name,
                        self.fullname,
                        type='autodoc',
                    )
            return False, selected
        if self.options.inherited_members:
            return False, list(members.values())
        return False, [m for m in members.values() if m.class_ == self.object]

    def get_doc(self) -> list[list[str]] | None:
        if isinstance(self.object, TypeVar):
            if self.object.__doc__ == TypeVar.__doc__:
                return []
        if self.doc_as_attr:
            # Don't show the docstring of the class when it is an alias.
            if self.get_variable_comment():
                return []
            return None

        lines = getattr(self, '_new_docstrings', None)
        if lines is not None:
            return lines

        classdoc_from = self.options.get(
            'class-doc-from', self.config.autoclass_content
        )

        docstrings = []
        attrdocstring = getdoc(self.object, self.get_attr)
        if attrdocstring:
            docstrings.append(attrdocstring)

        # for classes, what the "docstring" is can be controlled via a
        # config value; the default is only the class docstring
        if classdoc_from in {'both', 'init'}:
            __init__ = self.get_attr(self.object, '__init__', None)
            initdocstring = getdoc(
                __init__,
                self.get_attr,
                self.config.autodoc_inherit_docstrings,
                self.object,
                '__init__',
            )
            # for new-style classes, no __init__ means default __init__
            if initdocstring is not None and (
                initdocstring == object.__init__.__doc__  # for pypy
                or initdocstring.strip() == object.__init__.__doc__  # for !pypy
            ):
                initdocstring = None
            if not initdocstring:
                # try __new__
                __new__ = self.get_attr(self.object, '__new__', None)
                initdocstring = getdoc(
                    __new__,
                    self.get_attr,
                    self.config.autodoc_inherit_docstrings,
                    self.object,
                    '__new__',
                )
                # for new-style classes, no __new__ means default __new__
                if initdocstring is not None and (
                    initdocstring == object.__new__.__doc__  # for pypy
                    or initdocstring.strip() == object.__new__.__doc__  # for !pypy
                ):
                    initdocstring = None
            if initdocstring:
                if classdoc_from == 'init':
                    docstrings = [initdocstring]
                else:
                    docstrings.append(initdocstring)

        tab_width = self.directive.state.document.settings.tab_width
        return [prepare_docstring(docstring, tab_width) for docstring in docstrings]

    def get_variable_comment(self) -> list[str] | None:
        try:
            key = ('', '.'.join(self.objpath))
            if self.doc_as_attr:
                analyzer = ModuleAnalyzer.for_module(self.modname)
            else:
                analyzer = ModuleAnalyzer.for_module(self.get_real_modname())
            analyzer.analyze()
            return list(analyzer.attr_docs.get(key, []))
        except PycodeError:
            return None

    def add_content(self, more_content: StringList | None) -> None:
        mode = _get_render_mode(self.config.autodoc_typehints_format)
        short_literals = self.config.python_display_short_literal_types

        if isinstance(self.object, NewType):
            supertype = restify(self.object.__supertype__, mode=mode)

            more_content = StringList([_('alias of %s') % supertype, ''], source='')
        if isinstance(self.object, TypeVar):
            attrs = [repr(self.object.__name__)]
            attrs.extend(
                stringify_annotation(constraint, mode, short_literals=short_literals)
                for constraint in self.object.__constraints__
            )
            if self.object.__bound__:
                bound = restify(self.object.__bound__, mode=mode)
                attrs.append(r'bound=\ ' + bound)
            if self.object.__covariant__:
                attrs.append('covariant=True')
            if self.object.__contravariant__:
                attrs.append('contravariant=True')

            more_content = StringList(
                [_('alias of TypeVar(%s)') % ', '.join(attrs), ''], source=''
            )
        if self.doc_as_attr and self.modname != self.get_real_modname():
            try:
                # override analyzer to obtain doccomment around its definition.
                self.analyzer = ModuleAnalyzer.for_module(self.modname)
                self.analyzer.analyze()
            except PycodeError:
                pass

        if self.doc_as_attr and not self.get_variable_comment():
            try:
                alias = restify(self.object, mode=mode)
                more_content = StringList([_('alias of %s') % alias], source='')
            except AttributeError:
                pass  # Invalid class object is passed.

        super().add_content(more_content)

    def document_members(self, all_members: bool = False) -> None:
        if self.doc_as_attr:
            return
        super().document_members(all_members)

    def generate(
        self,
        more_content: StringList | None = None,
        real_modname: str | None = None,
        check_module: bool = False,
        all_members: bool = False,
    ) -> None:
        # Do not pass real_modname and use the name from the __module__
        # attribute of the class.
        # If a class gets imported into the module real_modname
        # the analyzer won't find the source of the class, if
        # it looks in real_modname.
        return super().generate(
            more_content=more_content,
            check_module=check_module,
            all_members=all_members,
        )


class ExceptionDocumenter(ClassDocumenter):
    """Specialized ClassDocumenter subclass for exceptions."""

    objtype = 'exception'
    member_order = 10

    # needs a higher priority than ClassDocumenter
    priority = ClassDocumenter.priority + 5

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        try:
            return isinstance(member, type) and issubclass(member, BaseException)
        except TypeError as exc:
            # It's possible for a member to be considered a type, but fail
            # issubclass checks due to not being a class. For example:
            # https://github.com/sphinx-doc/sphinx/issues/11654#issuecomment-1696790436
            msg = (
                f'{cls.__name__} failed to discern if member {member} with'
                f' membername {membername} is a BaseException subclass.'
            )
            raise ValueError(msg) from exc


class DataDocumenterMixinBase:
    # define types of instance variables
    config: sphinx.config.Config
    env: BuildEnvironment
    modname: str
    parent: Any
    object: Any
    objpath: list[str]

    def should_suppress_directive_header(self) -> bool:
        """Check directive header should be suppressed."""
        return False

    def should_suppress_value_header(self) -> bool:
        """Check :value: header should be suppressed."""
        return False

    def update_content(self, more_content: StringList) -> None:
        """Update docstring, for example with TypeVar variance."""
        pass


class GenericAliasMixin(DataDocumenterMixinBase):
    """Mixin for DataDocumenter and AttributeDocumenter to provide the feature for
    supporting GenericAliases.
    """

    def should_suppress_directive_header(self) -> bool:
        return (
            inspect.isgenericalias(self.object)
            or super().should_suppress_directive_header()
        )

    def update_content(self, more_content: StringList) -> None:
        if inspect.isgenericalias(self.object):
            mode = _get_render_mode(self.config.autodoc_typehints_format)
            alias = restify(self.object, mode=mode)

            more_content.append(_('alias of %s') % alias, '')
            more_content.append('', '')

        super().update_content(more_content)


class UninitializedGlobalVariableMixin(DataDocumenterMixinBase):
    """Mixin for DataDocumenter to provide the feature for supporting uninitialized
    (type annotation only) global variables.
    """

    def import_object(self, raiseerror: bool = False) -> bool:
        try:
            return super().import_object(raiseerror=True)  # type: ignore[misc]
        except ImportError as exc:
            # annotation only instance variable (PEP-526)
            try:
                with mock(self.config.autodoc_mock_imports):
                    parent = import_module(self.modname)
                    annotations = get_type_hints(
                        parent,
                        None,
                        self._type_alias_namespace(parent),
                        include_extras=True,
                    )
                    if self.objpath[-1] in annotations:
                        self.object = UNINITIALIZED_ATTR
                        self.parent = parent
                        return True
            except ImportError:
                pass

            if raiseerror:
                raise
            logger.warning(exc.args[0], type='autodoc', subtype='import_object')
            self.env.note_reread()
            return False

    def should_suppress_value_header(self) -> bool:
        return (
            self.object is UNINITIALIZED_ATTR or super().should_suppress_value_header()
        )

    def get_doc(self) -> list[list[str]] | None:
        if self.object is UNINITIALIZED_ATTR:
            return []
        return super().get_doc()  # type: ignore[misc]


class DataDocumenter(
    GenericAliasMixin, UninitializedGlobalVariableMixin, ModuleLevelDocumenter
):
    """Specialized Documenter subclass for data items."""

    objtype = 'data'
    member_order = 40
    priority = -10
    option_spec: ClassVar[dict[str, Any]] = dict(ModuleLevelDocumenter.option_spec)
    option_spec['annotation'] = annotation_option
    option_spec['no-value'] = bool_option

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        return isinstance(parent, ModuleDocumenter) and isattr

    def update_annotations(self, parent: Any) -> None:
        """Update __annotations__ to support type_comment and so on."""
        annotations = dict(inspect.getannotations(parent))
        parent.__annotations__ = annotations

        try:
            analyzer = ModuleAnalyzer.for_module(self.modname)
            analyzer.analyze()
            for (classname, attrname), annotation in analyzer.annotations.items():
                if not classname and attrname not in annotations:
                    annotations[attrname] = annotation
        except PycodeError:
            pass

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)
        if self.parent:
            self.update_annotations(self.parent)

        return ret

    def should_suppress_value_header(self) -> bool:
        if super().should_suppress_value_header():
            return True
        doc = self.get_doc() or []
        docstring, metadata = separate_metadata(
            '\n'.join(functools.reduce(operator.iadd, doc, []))
        )
        if 'hide-value' in metadata:
            return True

        return False

    def add_directive_header(self, sig: str) -> None:
        super().add_directive_header(sig)
        sourcename = self.get_sourcename()
        if (
            self.options.annotation is SUPPRESS
            or self.should_suppress_directive_header()
        ):
            pass
        elif self.options.annotation:
            self.add_line('   :annotation: %s' % self.options.annotation, sourcename)
        else:
            if self.config.autodoc_typehints != 'none':
                # obtain annotation for this data
                annotations = get_type_hints(
                    self.parent,
                    None,
                    self._type_alias_namespace(self.parent),
                    include_extras=True,
                )
                if self.objpath[-1] in annotations:
                    mode = _get_render_mode(self.config.autodoc_typehints_format)
                    short_literals = self.config.python_display_short_literal_types
                    objrepr = stringify_annotation(
                        annotations.get(self.objpath[-1]),
                        mode,
                        short_literals=short_literals,
                    )
                    self.add_line('   :type: ' + objrepr, sourcename)

            try:
                if (
                    self.options.no_value
                    or self.should_suppress_value_header()
                    or ismock(self.object)
                ):
                    pass
                else:
                    objrepr = object_description(self.object)
                    self.add_line('   :value: ' + objrepr, sourcename)
            except ValueError:
                pass

    def document_members(self, all_members: bool = False) -> None:
        pass

    def get_real_modname(self) -> str:
        real_modname = self.get_attr(self.parent or self.object, '__module__', None)
        return real_modname or self.modname

    def get_module_comment(self, attrname: str) -> list[str] | None:
        try:
            analyzer = ModuleAnalyzer.for_module(self.modname)
            analyzer.analyze()
            key = ('', attrname)
            if key in analyzer.attr_docs:
                return list(analyzer.attr_docs[key])
        except PycodeError:
            pass

        return None

    def get_doc(self) -> list[list[str]] | None:
        # Check the variable has a docstring-comment
        comment = self.get_module_comment(self.objpath[-1])
        if comment:
            return [comment]
        return super().get_doc()

    def add_content(self, more_content: StringList | None) -> None:
        # Disable analyzing variable comment on Documenter.add_content() to control it on
        # DataDocumenter.add_content()
        self.analyzer = None

        if not more_content:
            more_content = StringList()

        self.update_content(more_content)
        super().add_content(more_content)


class MethodDocumenter(DocstringSignatureMixin, ClassLevelDocumenter):  # type: ignore[misc]
    """Specialized Documenter subclass for methods (normal, static and class)."""

    objtype = 'method'
    directivetype = 'method'
    member_order = 50
    priority = 1  # must be more than FunctionDocumenter

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        return inspect.isroutine(member) and not isinstance(parent, ModuleDocumenter)

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)
        if not ret:
            return ret

        # to distinguish classmethod/staticmethod
        obj = self.parent.__dict__.get(self.object_name, self.object)
        if inspect.isstaticmethod(obj, cls=self.parent, name=self.object_name):
            # document static members before regular methods
            self.member_order -= 1
        elif inspect.isclassmethod(obj):
            # document class methods before static methods as
            # they usually behave as alternative constructors
            self.member_order -= 2
        return ret

    def format_args(self, **kwargs: Any) -> str:
        if self.config.autodoc_typehints in {'none', 'description'}:
            kwargs.setdefault('show_annotation', False)
        if self.config.autodoc_typehints_format == 'short':
            kwargs.setdefault('unqualified_typehints', True)
        if self.config.python_display_short_literal_types:
            kwargs.setdefault('short_literals', True)

        # -----------------------------------------------------------------
        # Issue #9976: Support the _sage_argspec_ attribute which makes it
        # possible to get argument specification of decorated callables in
        # documentation correct.  See e.g. sage.misc.decorators.sage_wraps.
        #
        # Note, however, that sage.misc.sageinspect.sage_getargspec already
        # uses a method _sage_argspec_, that only works on objects, not on
        # classes, though.
        obj = self.object
        if obj == object.__init__ and self.parent != object:  # NoQA: E721
            # Classes not having own __init__() method are shown as no arguments.
            #
            # Note: The signature of object.__init__() is
            # (self, /, *args, **kwargs).  But it makes users confused.
            args = '()'
        else:
            if hasattr(obj, "_sage_argspec_"):
                argspec = obj._sage_argspec_()
            elif inspect.isbuiltin(obj) or inspect.ismethoddescriptor(obj):
                # can never get arguments of a C function or method unless
                # a function to do so is supplied
                argspec = sage_getargspec(obj)
            else:
                # The check above misses ordinary Python methods in Cython
                # files.
                argspec = sage_getargspec(obj)

            if argspec is not None and argspec[0] and argspec[0][0] in ('cls', 'self'):
                del argspec[0][0]
            if argspec is None:
                return None
            args = sage_formatargspec(*argspec)
        # -----------------------------------------------------------------

        if self.config.strip_signature_backslash:
            # escape backslashes for reST
            args = args.replace('\\', '\\\\')
        return args

    def add_directive_header(self, sig: str) -> None:
        super().add_directive_header(sig)

        sourcename = self.get_sourcename()
        obj = self.parent.__dict__.get(self.object_name, self.object)
        if inspect.isabstractmethod(obj):
            self.add_line('   :abstractmethod:', sourcename)
        if inspect.iscoroutinefunction(obj) or inspect.isasyncgenfunction(obj):
            self.add_line('   :async:', sourcename)
        if (
            inspect.is_classmethod_like(obj)
            or inspect.is_singledispatch_method(obj)
            and inspect.is_classmethod_like(obj.func)
        ):
            self.add_line('   :classmethod:', sourcename)
        if inspect.isstaticmethod(obj, cls=self.parent, name=self.object_name):
            self.add_line('   :staticmethod:', sourcename)
        if self.analyzer and '.'.join(self.objpath) in self.analyzer.finals:
            self.add_line('   :final:', sourcename)

    def document_members(self, all_members: bool = False) -> None:
        pass

    # ------------------------------------------------------------------------
    # Issue #34730: The format_signature() of the class MethodDocumenter
    # supports overloaded methods via inspect.signature(), which does not work
    # with Sage yet. Hence the method was removed from here.
    # ------------------------------------------------------------------------

    def merge_default_value(self, actual: Signature, overload: Signature) -> Signature:
        """Merge default values of actual implementation to the overload variants."""
        parameters = list(overload.parameters.values())
        for i, param in enumerate(parameters):
            actual_param = actual.parameters.get(param.name)
            if actual_param and param.default == '...':
                parameters[i] = param.replace(default=actual_param.default)

        return overload.replace(parameters=parameters)

    def annotate_to_first_argument(
        self, func: Callable[..., Any], typ: type
    ) -> Callable[..., Any] | None:
        """Annotate type hint to the first argument of function if needed."""
        try:
            sig = inspect.signature(func, type_aliases=self._type_aliases(func))
        except TypeError as exc:
            logger.warning(
                __('Failed to get a method signature for %s: %s'), self.fullname, exc
            )
            return None
        except ValueError:
            return None

        if len(sig.parameters) == 1:
            return None

        def dummy():  # type: ignore[no-untyped-def]  # NoQA: ANN202
            pass

        params = list(sig.parameters.values())
        if params[1].annotation is Parameter.empty:
            params[1] = params[1].replace(annotation=typ)
            try:
                dummy.__signature__ = sig.replace(  # type: ignore[attr-defined]
                    parameters=params
                )
                return dummy
            except (AttributeError, TypeError):
                # failed to update signature (ex. built-in or extension types)
                return None

        return func

    def get_doc(self) -> list[list[str]] | None:
        if self._new_docstrings is not None:
            # docstring already returned previously, then modified by
            # `DocstringSignatureMixin`.  Just return the previously-computed
            # result, so that we don't lose the processing done by
            # `DocstringSignatureMixin`.
            return self._new_docstrings
        if self.objpath[-1] == '__init__':
            docstring = getdoc(
                self.object,
                self.get_attr,
                self.config.autodoc_inherit_docstrings,
                self.parent,
                self.object_name,
            )
            if docstring is not None and (
                docstring == object.__init__.__doc__  # for pypy
                or docstring.strip() == object.__init__.__doc__  # for !pypy
            ):
                docstring = None
            if docstring:
                tab_width = self.directive.state.document.settings.tab_width
                return [prepare_docstring(docstring, tabsize=tab_width)]
            return []
        if self.objpath[-1] == '__new__':
            docstring = getdoc(
                self.object,
                self.get_attr,
                self.config.autodoc_inherit_docstrings,
                self.parent,
                self.object_name,
            )
            if docstring is not None and (
                docstring == object.__new__.__doc__  # for pypy
                or docstring.strip() == object.__new__.__doc__  # for !pypy
            ):
                docstring = None
            if docstring:
                tab_width = self.directive.state.document.settings.tab_width
                return [prepare_docstring(docstring, tabsize=tab_width)]
            return []
        return super().get_doc()


class NonDataDescriptorMixin(DataDocumenterMixinBase):
    """Mixin for AttributeDocumenter to provide the feature for supporting non
    data-descriptors.

    .. note:: This mix-in must be inherited after other mix-ins.  Otherwise, docstring
              and :value: header will be suppressed unexpectedly.
    """

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)  # type: ignore[misc]
        if ret and not inspect.isattributedescriptor(self.object):
            self.non_data_descriptor = True
        else:
            self.non_data_descriptor = False

        return ret

    def should_suppress_value_header(self) -> bool:
        return (
            not getattr(self, 'non_data_descriptor', False)
            or super().should_suppress_directive_header()
        )

    def get_doc(self) -> list[list[str]] | None:
        if getattr(self, 'non_data_descriptor', False):
            # the docstring of non datadescriptor is very probably the wrong thing
            # to display
            return None
        return super().get_doc()  # type: ignore[misc]


class SlotsMixin(DataDocumenterMixinBase):
    """Mixin for AttributeDocumenter to provide the feature for supporting __slots__."""

    def isslotsattribute(self) -> bool:
        """Check the subject is an attribute in __slots__."""
        try:
            if parent___slots__ := inspect.getslots(self.parent):
                return self.objpath[-1] in parent___slots__
            return False
        except (ValueError, TypeError):
            return False

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)  # type: ignore[misc]
        if self.isslotsattribute():
            self.object = SLOTSATTR

        return ret

    def should_suppress_value_header(self) -> bool:
        if self.object is SLOTSATTR:
            return True
        return super().should_suppress_value_header()

    def get_doc(self) -> list[list[str]] | None:
        if self.object is SLOTSATTR:
            try:
                parent___slots__ = inspect.getslots(self.parent)
                if parent___slots__ and (
                    docstring := parent___slots__.get(self.objpath[-1])
                ):
                    docstring = prepare_docstring(docstring)
                    return [docstring]
                return []
            except ValueError as exc:
                logger.warning(
                    __('Invalid __slots__ found on %s. Ignored.'),
                    (self.parent.__qualname__, exc),
                    type='autodoc',
                )
                return []
        else:
            return super().get_doc()  # type: ignore[misc]


class RuntimeInstanceAttributeMixin(DataDocumenterMixinBase):
    """Mixin for AttributeDocumenter to provide the feature for supporting runtime
    instance attributes (that are defined in __init__() methods with doc-comments).

    Example::

        class Foo:
            def __init__(self):
                self.attr = None  #: This is a target of this mix-in.
    """

    RUNTIME_INSTANCE_ATTRIBUTE = object()

    def is_runtime_instance_attribute(self, parent: Any) -> bool:
        """Check the subject is an attribute defined in __init__()."""
        # An instance variable defined in __init__().
        if self.get_attribute_comment(parent, self.objpath[-1]):  # type: ignore[attr-defined]
            return True
        return self.is_runtime_instance_attribute_not_commented(parent)

    def is_runtime_instance_attribute_not_commented(self, parent: Any) -> bool:
        """Check the subject is an attribute defined in __init__() without comment."""
        for cls in inspect.getmro(parent):
            try:
                module = safe_getattr(cls, '__module__')
                qualname = safe_getattr(cls, '__qualname__')

                analyzer = ModuleAnalyzer.for_module(module)
                analyzer.analyze()
                if qualname and self.objpath:
                    key = f'{qualname}.{self.objpath[-1]}'
                    if key in analyzer.tagorder:
                        return True
            except (AttributeError, PycodeError):
                pass

        return False

    def import_object(self, raiseerror: bool = False) -> bool:
        """Check the existence of runtime instance attribute after failing to import the
        attribute.
        """
        try:
            return super().import_object(raiseerror=True)  # type: ignore[misc]
        except ImportError as exc:
            try:
                with mock(self.config.autodoc_mock_imports):
                    ret = import_object(
                        self.modname,
                        self.objpath[:-1],
                        'class',
                        attrgetter=self.get_attr,  # type: ignore[attr-defined]
                    )
                    parent = ret[3]
                    if self.is_runtime_instance_attribute(parent):
                        self.object = self.RUNTIME_INSTANCE_ATTRIBUTE
                        self.parent = parent
                        return True
            except ImportError:
                pass

            if raiseerror:
                raise
            logger.warning(exc.args[0], type='autodoc', subtype='import_object')
            self.env.note_reread()
            return False

    def should_suppress_value_header(self) -> bool:
        return (
            self.object is self.RUNTIME_INSTANCE_ATTRIBUTE
            or super().should_suppress_value_header()
        )

    def get_doc(self) -> list[list[str]] | None:
        if (
            self.object is self.RUNTIME_INSTANCE_ATTRIBUTE
            and self.is_runtime_instance_attribute_not_commented(self.parent)
        ):
            return None
        return super().get_doc()  # type: ignore[misc]


class UninitializedInstanceAttributeMixin(DataDocumenterMixinBase):
    """Mixin for AttributeDocumenter to provide the feature for supporting uninitialized
    instance attributes (PEP-526 styled, annotation only attributes).

    Example::

        class Foo:
            attr: int  #: This is a target of this mix-in.
    """

    def is_uninitialized_instance_attribute(self, parent: Any) -> bool:
        """Check the subject is an annotation only attribute."""
        annotations = get_type_hints(
            parent, None, self._type_alias_namespace(parent), include_extras=True
        )
        return self.objpath[-1] in annotations

    def import_object(self, raiseerror: bool = False) -> bool:
        """Check the existence of uninitialized instance attribute when failed to import
        the attribute.
        """
        try:
            return super().import_object(raiseerror=True)  # type: ignore[misc]
        except ImportError as exc:
            try:
                ret = import_object(
                    self.modname,
                    self.objpath[:-1],
                    'class',
                    attrgetter=self.get_attr,  # type: ignore[attr-defined]
                )
                parent = ret[3]
                if self.is_uninitialized_instance_attribute(parent):
                    self.object = UNINITIALIZED_ATTR
                    self.parent = parent
                    return True
            except ImportError:
                pass

            if raiseerror:
                raise
            logger.warning(exc.args[0], type='autodoc', subtype='import_object')
            self.env.note_reread()
            return False

    def should_suppress_value_header(self) -> bool:
        return (
            self.object is UNINITIALIZED_ATTR or super().should_suppress_value_header()
        )

    def get_doc(self) -> list[list[str]] | None:
        if self.object is UNINITIALIZED_ATTR:
            return None
        return super().get_doc()  # type: ignore[misc]


class AttributeDocumenter(  # type: ignore[misc]
    GenericAliasMixin,
    SlotsMixin,
    RuntimeInstanceAttributeMixin,
    UninitializedInstanceAttributeMixin,
    NonDataDescriptorMixin,
    DocstringStripSignatureMixin,
    ClassLevelDocumenter,
):
    """Specialized Documenter subclass for attributes."""

    objtype = 'attribute'
    member_order = 60
    option_spec: ClassVar[dict[str, Any]] = dict(ModuleLevelDocumenter.option_spec)
    option_spec['annotation'] = annotation_option
    option_spec['no-value'] = bool_option

    # must be higher than the MethodDocumenter, else it will recognize
    # some non-data descriptors as methods
    priority = 10

    @staticmethod
    def is_function_or_method(obj: Any) -> bool:
        return (
            inspect.isfunction(obj) or inspect.isbuiltin(obj) or inspect.ismethod(obj)
        )

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        if isinstance(parent, ModuleDocumenter):
            return False
        # ---------------------------------------------------------------------
        # Issue #34730: Do not pass objects of the class CachedMethodCaller as
        # attributes.
        #
        #   sage: from sphinx.util import inspect
        #   sage: A = AlgebrasWithBasis(QQ).example()
        #   sage: member = A.one_basis
        #   sage: type(member)
        #   <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
        #   sage: inspect.isattributedescriptor(member)
        #   True
        #   sage: inspect.isroutine(member)
        #   True
        #
        if inspect.isattributedescriptor(member) and not inspect.isroutine(member):
            return True
        # Issue #26522: Pass objects of classes that inherit ClasscallMetaclass
        # as attributes rather than method descriptors.
        from sage.misc.classcall_metaclass import ClasscallMetaclass
        if isinstance(type(member), ClasscallMetaclass):
            return True
        # ---------------------------------------------------------------------
        return not inspect.isroutine(member) and not isinstance(member, type)

    def document_members(self, all_members: bool = False) -> None:
        pass

    def update_annotations(self, parent: Any) -> None:
        """Update __annotations__ to support type_comment and so on."""
        try:
            annotations = dict(inspect.getannotations(parent))
            parent.__annotations__ = annotations

            for cls in inspect.getmro(parent):
                try:
                    module = safe_getattr(cls, '__module__')
                    qualname = safe_getattr(cls, '__qualname__')

                    analyzer = ModuleAnalyzer.for_module(module)
                    analyzer.analyze()
                    anns = analyzer.annotations
                    for (classname, attrname), annotation in anns.items():
                        if classname == qualname and attrname not in annotations:
                            annotations[attrname] = annotation
                except (AttributeError, PycodeError):
                    pass
        except (AttributeError, TypeError):
            # Failed to set __annotations__ (built-in, extensions, etc.)
            pass

    def import_object(self, raiseerror: bool = False) -> bool:
        ret = super().import_object(raiseerror)
        if inspect.isenumattribute(self.object):
            self.object = self.object.value
        if self.parent:
            self.update_annotations(self.parent)

        return ret

    def get_real_modname(self) -> str:
        real_modname = self.get_attr(self.parent or self.object, '__module__', None)
        return real_modname or self.modname

    def should_suppress_value_header(self) -> bool:
        if super().should_suppress_value_header():
            return True
        doc = self.get_doc()
        if doc:
            docstring, metadata = separate_metadata(
                '\n'.join(functools.reduce(operator.iadd, doc, []))
            )
            if 'hide-value' in metadata:
                return True
        return False

    def add_directive_header(self, sig: str) -> None:
        super().add_directive_header(sig)
        sourcename = self.get_sourcename()
        if (
            self.options.annotation is SUPPRESS
            or self.should_suppress_directive_header()
        ):
            pass
        elif self.options.annotation:
            self.add_line('   :annotation: %s' % self.options.annotation, sourcename)
        else:
            if self.config.autodoc_typehints != 'none':
                # obtain type annotation for this attribute
                annotations = get_type_hints(
                    self.parent,
                    None,
                    self._type_alias_namespace(self.parent),
                    include_extras=True,
                )
                if self.objpath[-1] in annotations:
                    mode = _get_render_mode(self.config.autodoc_typehints_format)
                    short_literals = self.config.python_display_short_literal_types
                    objrepr = stringify_annotation(
                        annotations.get(self.objpath[-1]),
                        mode,
                        short_literals=short_literals,
                    )
                    self.add_line('   :type: ' + objrepr, sourcename)

            try:
                if (
                    self.options.no_value
                    or self.should_suppress_value_header()
                    or ismock(self.object)
                ):
                    pass
                else:
                    objrepr = object_description(self.object)
                    self.add_line('   :value: ' + objrepr, sourcename)
            except ValueError:
                pass

    def get_attribute_comment(self, parent: Any, attrname: str) -> list[str] | None:
        for cls in inspect.getmro(parent):
            try:
                module = safe_getattr(cls, '__module__')
                qualname = safe_getattr(cls, '__qualname__')

                analyzer = ModuleAnalyzer.for_module(module)
                analyzer.analyze()
                if qualname and self.objpath:
                    key = (qualname, attrname)
                    if key in analyzer.attr_docs:
                        return list(analyzer.attr_docs[key])
            except (AttributeError, PycodeError):
                pass

        return None

    def get_doc(self) -> list[list[str]] | None:
        # Check the attribute has a docstring-comment
        comment = self.get_attribute_comment(self.parent, self.objpath[-1])
        if comment:
            return [comment]

        try:
            # Disable `autodoc_inherit_docstring` temporarily to avoid to obtain
            # a docstring from the value which descriptor returns unexpectedly.
            # See: https://github.com/sphinx-doc/sphinx/issues/7805
            orig = self.config.autodoc_inherit_docstrings
            self.config.autodoc_inherit_docstrings = False
            return super().get_doc()
        finally:
            self.config.autodoc_inherit_docstrings = orig

    def add_content(self, more_content: StringList | None) -> None:
        # Disable analyzing attribute comment on Documenter.add_content() to control it on
        # AttributeDocumenter.add_content()
        self.analyzer = None

        if more_content is None:
            more_content = StringList()
        self.update_content(more_content)
        super().add_content(more_content)


class PropertyDocumenter(DocstringStripSignatureMixin, ClassLevelDocumenter):  # type: ignore[misc]
    """Specialized Documenter subclass for properties."""

    objtype = 'property'
    member_order = 60

    # before AttributeDocumenter
    priority = AttributeDocumenter.priority + 1

    @classmethod
    def can_document_member(
        cls: type[Documenter], member: Any, membername: str, isattr: bool, parent: Any
    ) -> bool:
        if isinstance(parent, ClassDocumenter):
            if inspect.isproperty(member):
                return True
            __dict__ = safe_getattr(parent.object, '__dict__', {})
            obj = __dict__.get(membername)
            return isinstance(obj, classmethod) and inspect.isproperty(obj.__func__)
        return False

    def import_object(self, raiseerror: bool = False) -> bool:
        """Check the existence of uninitialized instance attribute when failed to import
        the attribute.
        """
        ret = super().import_object(raiseerror)
        if ret and not inspect.isproperty(self.object):
            __dict__ = safe_getattr(self.parent, '__dict__', {})
            obj = __dict__.get(self.objpath[-1])
            if isinstance(obj, classmethod) and inspect.isproperty(obj.__func__):
                self.object = obj.__func__
                self.isclassmethod: bool = True
                return True
            return False

        self.isclassmethod = False
        return ret

    def format_args(self, **kwargs: Any) -> str:
        func = self._get_property_getter()
        if func is None:
            return ''

        # update the annotations of the property getter
        self._events.emit('autodoc-before-process-signature', func, False)
        # correctly format the arguments for a property
        return super().format_args(**kwargs)

    def document_members(self, all_members: bool = False) -> None:
        pass

    def get_real_modname(self) -> str:
        real_modname = self.get_attr(self.parent or self.object, '__module__', None)
        return real_modname or self.modname

    def add_directive_header(self, sig: str) -> None:
        super().add_directive_header(sig)
        sourcename = self.get_sourcename()
        if inspect.isabstractmethod(self.object):
            self.add_line('   :abstractmethod:', sourcename)
        if self.isclassmethod:
            self.add_line('   :classmethod:', sourcename)

        func = self._get_property_getter()
        if func is None or self.config.autodoc_typehints == 'none':
            return

        try:
            signature = inspect.signature(
                func, type_aliases=self._type_aliases(func)
            )
            if signature.return_annotation is not Parameter.empty:
                mode = _get_render_mode(self.config.autodoc_typehints_format)
                short_literals = self.config.python_display_short_literal_types
                objrepr = stringify_annotation(
                    signature.return_annotation, mode, short_literals=short_literals
                )
                self.add_line('   :type: ' + objrepr, sourcename)
        except TypeError as exc:
            logger.warning(
                __('Failed to get a function signature for %s: %s'), self.fullname, exc
            )
            pass
        except ValueError:
            pass

    def _get_property_getter(self) -> Callable[..., Any] | None:
        if safe_getattr(self.object, 'fget', None):  # property
            return self.object.fget
        if safe_getattr(self.object, 'func', None):  # cached_property
            return self.object.func
        return None


def autodoc_attrgetter(
    obj: Any, name: str, *defargs: Any, registry: SphinxComponentRegistry
) -> Any:
    """Alternative getattr() for types"""
    for typ, func in registry.autodoc_attrgetters.items():
        if isinstance(obj, typ):
            return func(obj, name, *defargs)

    return safe_getattr(obj, name, *defargs)


def setup(app: Sphinx) -> ExtensionMetadata:
    from sphinx.ext.autodoc.preserve_defaults import update_defvalue
    from sphinx.ext.autodoc.type_comment import update_annotations_using_type_comments
    from sphinx.ext.autodoc.typehints import _merge_typehints, record_typehints

    app.add_autodocumenter(ModuleDocumenter)
    app.add_autodocumenter(ClassDocumenter)
    app.add_autodocumenter(ExceptionDocumenter)
    app.add_autodocumenter(DataDocumenter)
    app.add_autodocumenter(FunctionDocumenter)
    app.add_autodocumenter(DecoratorDocumenter)
    app.add_autodocumenter(MethodDocumenter)
    app.add_autodocumenter(AttributeDocumenter)
    app.add_autodocumenter(PropertyDocumenter)

    app.add_config_value(
        'autoclass_content',
        'class',
        'env',
        types=ENUM('both', 'class', 'init'),
    )
    app.add_config_value(
        'autodoc_member_order',
        'alphabetical',
        'env',
        types=ENUM('alphabetical', 'bysource', 'groupwise'),
    )
    app.add_config_value(
        'autodoc_class_signature',
        'mixed',
        'env',
        types=ENUM('mixed', 'separated'),
    )
    app.add_config_value('autodoc_default_options', {}, 'env', types=frozenset({dict}))
    app.add_config_value(
        'autodoc_docstring_signature', True, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_mock_imports', [], 'env', types=frozenset({list, tuple})
    )
    app.add_config_value(
        'autodoc_typehints',
        'signature',
        'env',
        types=ENUM('signature', 'description', 'none', 'both'),
    )
    app.add_config_value(
        'autodoc_typehints_description_target',
        'all',
        'env',
        types=ENUM('all', 'documented', 'documented_params'),
    )
    app.add_config_value('autodoc_type_aliases', {}, 'env', types=frozenset({dict}))
    app.add_config_value(
        'autodoc_typehints_format',
        'short',
        'env',
        types=ENUM('fully-qualified', 'short'),
    )
    app.add_config_value('autodoc_warningiserror', True, 'env', types=frozenset({bool}))
    app.add_config_value(
        'autodoc_inherit_docstrings', True, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_preserve_defaults', False, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_use_type_comments', True, 'env', types=frozenset({bool})
    )
    # Sage's autodoc fork is the legacy class-based implementation.  Keep the
    # Sphinx 9 configuration value available, but this extension does not
    # implement the new dynamic autodoc branch.
    app.add_config_value(
        'autodoc_use_legacy_class_based', True, 'env', types=frozenset({bool})
    )
    app.add_event('autodoc-before-process-signature')
    app.add_event('autodoc-process-docstring')
    app.add_event('autodoc-process-signature')
    app.add_event('autodoc-skip-member')
    app.add_event('autodoc-process-bases')

    app.connect('object-description-transform', _merge_typehints)
    app.connect('autodoc-before-process-signature', update_defvalue)
    app.connect(
        'autodoc-before-process-signature', update_annotations_using_type_comments
    )
    app.connect('autodoc-process-signature', record_typehints)

    return {
        'version': sphinx.__display_version__,
        'parallel_read_safe': True,
    }
