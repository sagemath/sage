#!/usr/bin/env python3
# Copyright (C) 2026 Tobias Diez
# This file is dual-licensed under GPLv2+ and Apache-2.0 so it can be
# upstreamed to Cython if desired.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "cython",
# ]
# ///
"""Generate and validate typing stubs for Cython sources.

This script uses Cython's compiler frontend to parse ``.pyx`` and ``.pxd``
files into compiler ASTs, walks those trees to collect public symbols, and then
either renders matching ``.pyi`` files or checks existing stubs for drift.

At a high level, the workflow is:

1. Parse the Cython source into compiler nodes.
2. Extract public classes, functions, and module attributes.
3. Render a lightweight stub or compare an existing stub against the extracted
    symbol set using Python's :mod:`ast` module.

Relevant Cython documentation:

- https://cython.readthedocs.io/en/latest/src/userguide/language_basics.html
- https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html

Inspired by:

- https://github.com/Vizonex/CyStub
"""

from __future__ import annotations

import argparse
import ast
import io
import sys
from contextlib import redirect_stderr
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from Cython.Compiler.CythonScope import create_cython_scope
from Cython.Compiler.Errors import (
    CompileError,
    init_thread,
    open_listing_file,
)
from Cython.Compiler.ExprNodes import (
    AddNode,
    AttributeNode,
    BoolNode,
    DictNode,
    ExprNode,
    FloatNode,
    GeneralCallNode,
    ImportNode,
    IndexNode,
    IntBinopNode,
    IntNode,
    ListNode,
    NameNode,
    NoneNode,
    PowNode,
    SetNode,
    SimpleCallNode,
    TupleNode,
    TypecastNode,
    UnaryMinusNode,
    UnaryPlusNode,
    UnicodeNode,
)
from Cython.Compiler.Main import (
    CompilationOptions,
    Context,
    FileSourceDescriptor,
    default_options,
)
from Cython.Compiler.Nodes import (
    CClassDefNode,
    CFuncDefNode,
    DefNode,
    Node,
    PyClassDefNode,
    SingleAssignmentNode,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

    from Cython.Compiler.ModuleNode import ModuleNode


@dataclass(frozen=True, slots=True)
class DefaultValue:
    """Normalized representation of a default argument or attribute value."""

    repr: str | None = None

    @classmethod
    def from_cython(cls, node: ExprNode | None) -> DefaultValue | None:
        """Convert a Cython expression node into a normalized default value representation."""

        if node is None:
            return None
        if isinstance(node, BoolNode):
            return cls("True" if node.value else "False")
        if isinstance(node, NoneNode):
            return cls("None")
        if isinstance(node, IntNode):
            return cls(str(node.value))
        if isinstance(node, FloatNode):
            return cls(str(node.value))
        if isinstance(node, IntBinopNode):
            return cls(str(node))
        if isinstance(node, UnicodeNode):
            return cls(repr(node.value))
        if (
            isinstance(node, PowNode)
            or isinstance(node, AddNode)
            or isinstance(node, UnaryPlusNode)
            or isinstance(node, UnaryMinusNode)
        ):
            return cls("...")  # not handled in detail
        if (
            isinstance(node, SimpleCallNode)
            or isinstance(node, ImportNode)
            or isinstance(node, GeneralCallNode)
        ):
            return cls("...")  # e.g., object(), not handled in detail
        if isinstance(node, TypecastNode):
            return cls("...")  # not handled in detail
        if isinstance(node, TupleNode):
            return cls("...")  # not handled in detail
        if isinstance(node, DictNode):
            return cls("...")  # not handled in detail
        if isinstance(node, ListNode):
            return cls("...")  # not handled in detail
        if isinstance(node, SetNode):
            return cls("...")  # not handled in detail
        if isinstance(node, AttributeNode):
            return cls("...")  # not handled in detail
        if isinstance(node, IndexNode):
            return cls("...")  # not handled in detail
        if isinstance(node, NameNode):
            return cls(node.name)
        raise NotImplementedError(
            f"Unhandled Cython default value node: {node.__class__.__name__}"
        )

    @classmethod
    def from_python(cls, node: ast.AST | None) -> DefaultValue | None:
        """Convert a Python AST node into a normalized default value representation."""

        if node is None:
            return None
        try:
            text = ast.unparse(node)
        except Exception:
            text = node.__class__.__name__
        return cls(text)


@dataclass(frozen=True, slots=True)
class Param:
    """Description of one function parameter."""

    name: str
    kind: str  # pos, vararg, varkw, kwonly
    type: str | None
    default: DefaultValue | None


@dataclass(frozen=True, slots=True)
class Symbol:
    """Base record for a public symbol within a module or nested class path."""

    path: tuple[str, ...]
    name: str

    @property
    def dotted(self) -> str:
        """Return the fully qualified symbol name."""

        if not self.path:
            return self.name
        return ".".join((*self.path, self.name))


@dataclass(frozen=True, slots=True)
class FunctionSymbol(Symbol):
    """Public callable symbol together with its parameter list."""

    params: tuple[Param, ...]


@dataclass(frozen=True, slots=True)
class AttributeSymbol(Symbol):
    """Public attribute symbol together with any detected default value."""

    default: DefaultValue | None


@dataclass(frozen=True, slots=True)
class WriteResult:
    """Result of attempting to generate a stub file."""

    path: Path
    written: bool


def _iter_child_nodes(node: Node) -> Iterator[Node]:
    """Yield child nodes using Cython's ``child_attrs`` metadata."""

    for attr in node.child_attrs:
        child = getattr(node, attr, None)
        if child is None:
            continue

        if isinstance(child, list):
            for item in child:
                if isinstance(item, Node):
                    yield item
        elif isinstance(child, Node):
            yield child


def _build_context(source: Path) -> Context:
    """Create a Cython compilation context for the given source file."""

    src_root = next(
        (
            parent
            for parent in source.resolve().parents
            if parent.name == "src" and (parent / "sage").is_dir()
        ),
        Path.cwd() / "src",
    )
    options = CompilationOptions(
        defaults=default_options,
        include_path=[str(source.parent), str(src_root)],
    )
    return Context.from_options(options)


def parse_cython_module(source: Path) -> ModuleNode:
    """Parse a .pyx/.pxd file into a Cython AST."""

    ctx = _build_context(source)
    init_thread()
    module_name = ctx.extract_module_name(str(source), None)
    scope = create_cython_scope(ctx)
    errors = io.StringIO()
    try:
        with redirect_stderr(errors):
            open_listing_file(path=None, echo_to_stderr=True)
            parsed = ctx.parse(
                source_desc=FileSourceDescriptor(str(source)),
                pxd=source.suffix == ".pxd",
                scope=scope,
                full_module_name=module_name,
            )
            parsed.scope = scope
            return parsed
    except CompileError as e:
        print(f"Error parsing {source}:", file=sys.stderr)
        details = errors.getvalue().strip()
        if details:
            print(details, file=sys.stderr)
        elif str(e):
            print(e, file=sys.stderr)
        raise


def _is_public_name(name: str) -> bool:
    """Filter out private helpers (leading double underscore) but allow special methods like __xyz__."""
    return not name.startswith("__") or name.endswith("__")


def _params_from_def(node: DefNode) -> tuple[Param, ...]:
    """Extract parameters from a Cython ``def`` node."""

    params: list[Param] = []

    for arg in node.args:
        default_node = arg.default or arg.default_value
        default = (
            DefaultValue.from_cython(default_node) if default_node is not None else None
        )
        params.append(
            Param(
                arg.declarator.name, "kwonly" if arg.kw_only else "pos", None, default
            )
        )

    if node.star_arg:
        params.append(Param(node.star_arg.name, "vararg", None, None))

    if node.starstar_arg:
        params.append(Param(node.starstar_arg.name, "varkw", None, None))
    return tuple(params)


def _params_from_cfunc(node: CFuncDefNode) -> tuple[Param, ...]:
    """Extract parameters from a Cython ``cdef`` or ``cpdef``."""

    py_func = node.py_func
    if isinstance(py_func, DefNode):
        return _params_from_def(py_func)
    if node.declarator and node.declarator.args:
        params: list[Param] = []
        for arg in node.declarator.args:
            default = DefaultValue.from_cython(arg.default) if arg.default else None
            if getattr(arg.declarator, "name", None):
                typ = None  # TODO: extract type info also for Python stubs, then use the following
                # typ = "bool" if arg.base_type.name == "bint" else None
                name = arg.declarator.name
            else:
                name = arg.base_type.name
                typ = None

            params.append(Param(name, "pos", typ, default))
        return tuple(params)
    return tuple()


def _collect_public_symbols(module: ModuleNode) -> set[Symbol]:
    """Walk a parsed Cython module and collect exported symbols."""

    symbols: set[Symbol] = set()

    def visit(node: Node, cls_path: tuple[str, ...]) -> None:
        if isinstance(node, (CClassDefNode, PyClassDefNode)):
            cls_name = (node.class_name if hasattr(node, "class_name") else None) or (
                node.name if hasattr(node, "name") else None
            )
            if cls_name is None:
                return
            new_path = (*cls_path, cls_name)
            for child in _iter_child_nodes(node):
                visit(child, new_path)
            return

        if isinstance(node, DefNode):
            if _is_public_name(node.name):
                symbols.add(FunctionSymbol(cls_path, node.name, _params_from_def(node)))
            return

        if isinstance(node, CFuncDefNode):
            name = node.declared_name()
            if (node.overridable or node.visibility == "public") and _is_public_name(
                name
            ):
                symbols.add(FunctionSymbol(cls_path, name, _params_from_cfunc(node)))
            return

        if isinstance(node, SingleAssignmentNode):
            lhs = node.lhs
            if lhs and lhs.is_name and _is_public_name(lhs.name):
                default = DefaultValue.from_cython(node.rhs)
                symbols.add(AttributeSymbol(cls_path, lhs.name, default))
            return

        for child in _iter_child_nodes(node):
            visit(child, cls_path)

    visit(module, ())
    return symbols


def _build_symbol_tree(symbols: Iterable[Symbol]) -> dict:
    """Create a nested mapping for module/classes → methods and attributes."""

    root: dict = {"methods": [], "attrs": [], "classes": {}}
    for sym in symbols:
        node = root
        for cls in sym.path:
            node = node.setdefault("classes", {}).setdefault(
                cls, {"methods": [], "attrs": [], "classes": {}}
            )
        if isinstance(sym, FunctionSymbol):
            node.setdefault("methods", []).append(sym)
        elif isinstance(sym, AttributeSymbol):
            node.setdefault("attrs", []).append(sym)
    return root


def _format_params(params: Sequence[Param]) -> str:
    """Format extracted parameters as a ``.pyi`` function signature fragment."""

    parts: list[str] = []
    for param in params:
        name = param.name
        if param.kind == "vararg":
            name = f"*{name}"
        elif param.kind == "varkw":
            name = f"**{name}"
        elif param.kind == "kwonly":
            name = f"*, {name}"
        if param.default is not None:
            name += f"={param.default.repr}"
        parts.append(name)
    return ", ".join(parts)


def _render_tree(tree: dict, indent: int = 0) -> list[str]:
    """Render the nested symbol tree into stub lines."""

    lines: list[str] = []
    current_indent = " " * indent

    for cls_name in tree.get("classes", {}):
        lines.append(f"{current_indent}class {cls_name}:")
        child = tree["classes"][cls_name]
        body = _render_tree(child, indent + 4)
        if not body:
            lines.append(" " * (indent + 4) + "pass")
        else:
            lines.extend(body)
        lines.append("")

    for attr in tree.get("attrs", []):
        suffix = " = ..." if attr.default is not None else ""
        lines.append(f"{current_indent}{attr.name}: Any{suffix}")

    for method in tree.get("methods", []):
        sig = _format_params(method.params)
        lines.append(f"{current_indent}def {method.name}({sig}) -> Any: ...")

    return lines


def render_stub(symbols: set[Symbol], source: Path) -> str:
    """Render the stub file content for the collected symbols."""

    lines: list[str] = [f"# Stubs generated from {source.name}"]
    if symbols:
        lines.append("from typing import Any")
        lines.append("")

    tree = _build_symbol_tree(symbols)
    lines.extend(_render_tree(tree))

    # Ensure a trailing newline and remove excessive blank lines
    while lines and lines[-1] == "":
        lines.pop()
    lines.append("")
    return "\n".join(lines)


def write_stub(
    source: Path, output_dir: Path | None, *, force: bool = False
) -> WriteResult:
    """Generate and write the stub for the given Cython source file."""

    target_dir = output_dir if output_dir is not None else source.parent
    target = target_dir / f"{source.stem}.pyi"
    if target.exists() and not force:
        return WriteResult(target, written=False)

    module = parse_cython_module(source)
    symbols = _collect_public_symbols(module)
    target_dir.mkdir(parents=True, exist_ok=True)
    target.write_text(render_stub(symbols, source), encoding="utf-8")
    return WriteResult(target, written=True)


def _params_from_stub_func(node: ast.AST) -> tuple[Param, ...]:
    """Extract parameter metadata from a Python stub function definition."""

    assert isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    params: list[Param] = []

    pos_args = list(node.args.posonlyargs) + list(node.args.args)
    defaults = node.args.defaults
    num_defaults = len(defaults)
    first_default_idx = len(pos_args) - num_defaults
    for idx, arg in enumerate(pos_args):
        default_node = (
            defaults[idx - first_default_idx] if idx >= first_default_idx else None
        )
        default = (
            DefaultValue.from_python(default_node) if default_node is not None else None
        )
        params.append(Param(arg.arg, "pos", None, default))

    if node.args.vararg is not None:
        params.append(Param(node.args.vararg.arg, "vararg", None, None))

    for kw_arg, kw_default in zip(node.args.kwonlyargs, node.args.kw_defaults):
        default = (
            DefaultValue.from_python(kw_default) if kw_default is not None else None
        )
        params.append(Param(kw_arg.arg, "kwonly", None, default))
    if node.args.kwarg is not None:
        params.append(Param(node.args.kwarg.arg, "varkw", None, None))

    return tuple(params)


def _collect_stub_symbols(stub_path: Path) -> set[Symbol]:
    """Collect public symbols from an existing ``.pyi`` file."""

    content = stub_path.read_text(encoding="utf-8")
    tree = ast.parse(content, filename=str(stub_path))
    symbols: set[Symbol] = set()
    class_stack: list[str] = []
    function_depth = 0

    class Visitor(ast.NodeVisitor):
        def visit_ClassDef(self, node: ast.ClassDef) -> None:  # noqa: N802
            class_stack.append(node.name)
            self.generic_visit(node)
            class_stack.pop()

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:  # noqa: N802
            nonlocal function_depth
            if _is_public_name(node.name):
                symbols.add(
                    FunctionSymbol(
                        tuple(class_stack), node.name, _params_from_stub_func(node)
                    )
                )
            function_depth += 1
            self.generic_visit(node)
            function_depth -= 1

        def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:  # noqa: N802
            nonlocal function_depth
            if _is_public_name(node.name):
                symbols.add(
                    FunctionSymbol(
                        tuple(class_stack), node.name, _params_from_stub_func(node)
                    )
                )
            function_depth += 1
            self.generic_visit(node)
            function_depth -= 1

        def visit_Assign(self, node: ast.Assign) -> None:  # noqa: N802
            if function_depth:
                return
            if len(node.targets) != 1:
                return
            target = node.targets[0]
            if isinstance(target, ast.Name) and _is_public_name(target.id):
                default = (
                    DefaultValue.from_python(node.value)
                    if node.value is not None
                    else None
                )
                symbols.add(
                    AttributeSymbol(tuple(class_stack), target.id, default=default)
                )

        def visit_AnnAssign(self, node: ast.AnnAssign) -> None:  # noqa: N802
            if function_depth:
                return
            target = node.target
            if isinstance(target, ast.Name) and _is_public_name(target.id):
                default = (
                    DefaultValue.from_python(node.value)
                    if node.value is not None
                    else None
                )
                symbols.add(
                    AttributeSymbol(tuple(class_stack), target.id, default=default)
                )

    Visitor().visit(tree)
    return symbols


def check_stub(source: Path, output_dir: Path | None) -> bool:
    """Compare a generated symbol model against the existing stub on disk."""

    expected = _collect_public_symbols(parse_cython_module(source))
    stub_path = (
        output_dir if output_dir is not None else source.parent
    ) / f"{source.stem}.pyi"
    if not stub_path.exists():
        # print(f"Missing stub file: {stub_path}", file=sys.stderr)
        return True  # don't treat missing stub as failure here (at least for now)

    found = _collect_stub_symbols(stub_path)
    exp_funcs = {
        (sym.path, sym.name): sym for sym in expected if isinstance(sym, FunctionSymbol)
    }
    found_funcs = {
        (sym.path, sym.name): sym for sym in found if isinstance(sym, FunctionSymbol)
    }

    exp_attrs = {
        (sym.path, sym.name): sym
        for sym in expected
        if isinstance(sym, AttributeSymbol)
    }
    found_attrs = {
        (sym.path, sym.name): sym for sym in found if isinstance(sym, AttributeSymbol)
    }

    missing = exp_funcs.keys() - found_funcs.keys()
    extra = found_funcs.keys() - exp_funcs.keys()

    param_mismatches: list[tuple[str, FunctionSymbol, FunctionSymbol]] = []
    for key in set(exp_funcs).intersection(found_funcs):
        if exp_funcs[key].params != found_funcs[key].params:
            dotted = ".".join((*key[0], key[1])) if key[0] else key[1]
            param_mismatches.append((dotted, exp_funcs[key], found_funcs[key]))

    attr_mismatches: list[tuple[str, DefaultValue | None, DefaultValue | None]] = []
    for key in set(exp_attrs).intersection(found_attrs):
        expected_default = exp_attrs[key].default
        found_default = found_attrs[key].default
        if expected_default != found_default and (
            expected_default is not None and expected_default.repr != "..."
        ):
            dotted = ".".join((*key[0], key[1])) if key[0] else key[1]
            attr_mismatches.append((dotted, expected_default, found_default))

    ok = True
    if missing:
        ok = False
        print(f"Missing functions in {stub_path}:")
        for parent, name in missing:
            print(f"  - {'.'.join((*parent, name)) if parent else name}")

    if extra:
        ok = False
        print(f"Extra functions in {stub_path}:")
        for parent, name in extra:
            print(f"  - {'.'.join((*parent, name)) if parent else name}")

    if param_mismatches:
        ok = False
        print(f"Param mismatches in {stub_path}:")
        for name, exp_func, found_func in sorted(param_mismatches):
            print(f"  - {name}: expected vs found")
            exp_sig = _format_params(exp_func.params)
            found_sig = _format_params(found_func.params)
            print(f"      ({exp_sig})")
            print(f"      ({found_sig})")

    if attr_mismatches:
        ok = False
        print(f"Attribute mismatches in {stub_path}:")
        for name, exp_default, found_default in sorted(attr_mismatches):
            print(f"  - {name}: expected vs found")
            exp_def = exp_default.repr if exp_default is not None else "<no default>"
            found_def = (
                found_default.repr if found_default is not None else "<no default>"
            )
            print(f"      {exp_def}")
            print(f"      {found_def}")

    return ok


def _iter_source_paths(inputs: Sequence[Path]) -> list[Path]:
    """Expand file and directory inputs into the list of ``.pyx`` sources."""

    paths: list[Path] = []
    for inp in inputs:
        if inp.is_file() and inp.suffix == ".pyx":
            paths.append(inp)
        elif inp.is_dir():
            for child in sorted(inp.rglob("*.pyx")):
                paths.append(child)
        else:
            raise FileNotFoundError(f"Unsupported source path: {inp}")
    return paths


def _parse_args(argv: Sequence[str]) -> argparse.Namespace:
    """Parse command-line arguments for the ``write`` and ``check`` subcommands."""

    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "sources", nargs="+", type=Path, help=".pyx files or directories to process"
    )
    common.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated/checked .pyi files (defaults to alongside source)",
    )

    write_parser = sub.add_parser(
        "write", parents=[common], help="Generate .pyi stubs from Cython sources"
    )
    write_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing .pyi files",
    )
    sub.add_parser(
        "check",
        parents=[common],
        help="Verify .pyi files match public symbols and default values",
    )

    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""

    args = _parse_args(list(argv) if argv is not None else sys.argv[1:])
    sources = _iter_source_paths(args.sources)

    if args.command == "write":
        ok = True
        for src in sources:
            try:
                result = write_stub(src, args.output_dir, force=args.force)
            except CompileError:
                ok = False
                continue
            if result.written:
                print(f"Wrote {result.path}")
            else:
                print(f"Skipped existing {result.path}")
        return 0 if ok else 1

    if args.command == "check":
        ok = True
        for src in sources:
            try:
                ok &= check_stub(src, args.output_dir)
            except CompileError as e:
                print(f"Error checking {src}: {e}")
                ok = False
        if ok:
            print("All stubs are up to date.")
        return 0 if ok else 1

    raise AssertionError("Unhandled command")


if __name__ == "__main__":
    raise SystemExit(main())
