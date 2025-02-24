#!/usr/bin/env python3

import argparse
import os
from argparse import Namespace
from pathlib import Path

from mesonbuild import mlog
from mesonbuild.ast import (
    AstPrinter,
    AstVisitor,
)
from mesonbuild.ast.interpreter import MethodNode
from mesonbuild.mformat import (
    run as meson_format,
)
from mesonbuild.mparser import (
    AssignmentNode,
    BaseNode,
    DictNode,
    SymbolNode,
)
from mesonbuild.rewriter import (
    ArgumentNode,
    ArrayNode,
    FunctionNode,
    Rewriter,
    StringNode,
    Token,
)

# Get target directory from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "sourcedir", help="Source directory", nargs="?", default=".", type=Path
)
options = parser.parse_args()


class AstPython(AstVisitor):
    install_sources_calls: list[MethodNode] = []
    extension_data: list[AssignmentNode] = []

    def visit_MethodNode(self, node: MethodNode) -> None:
        if node.name.value == "install_sources":
            self.install_sources_calls += [node]
        return super().visit_MethodNode(node)

    def visit_AssignmentNode(self, node: AssignmentNode) -> None:
        if node.var_name.value in ["extension_data", "extension_data_cpp"]:
            self.extension_data += [node]
        return super().visit_AssignmentNode(node)


# Utility function to get a list of the sources from a node
def arg_list_from_node(n):
    args = []
    if isinstance(n, FunctionNode) or isinstance(n, MethodNode):
        args = list(n.args.arguments)
    # if 'func_name' in n and n.func_name.value in BUILD_TARGET_FUNCTIONS:
    #     args.pop(0)
    elif isinstance(n, ArrayNode):
        args = n.args.arguments
    elif isinstance(n, ArgumentNode):
        args = n.arguments
    return args


def _symbol(val: str) -> SymbolNode:
    return SymbolNode(Token("", "", 0, 0, 0, (0, 0), val))


def update_python_sources(self: Rewriter, visitor: AstPython):
    for target in visitor.install_sources_calls:
        # Generate the current source list
        src_list: list[str] = []
        for arg in arg_list_from_node(target):
            if isinstance(arg, StringNode):
                src_list += [arg.value]

        folder = Path(target.filename).parent
        python_files = sorted(
            list(folder.glob("*.py")) + list(folder.glob('*.pxd'))
        )  # + list(folder.glob('*.pxd')) + list(folder.glob('*.h')))

        to_append: list[StringNode] = []
        for file in python_files:
            file_name = file.name
            if file_name == "__init__.py":
                # We don't want to add __init__.py files
                continue
            if file_name in src_list:
                continue
            token = Token("string", target.filename, 0, 0, 0, None, file_name)
            to_append += [StringNode(token)]

        # Get all deleted files
        to_remove = []
        for src in src_list:
            if not folder.joinpath(src).exists():
                to_remove += [src]

        if not to_append and not to_remove:
            continue

        # Update the source list
        target.args.arguments = sorted(
            [
                arg
                for arg in target.args.arguments
                if not (isinstance(arg, StringNode) and arg.value in to_remove)
            ]
            + to_append,
            key=lambda x: x.value,
        )

        # Mark the node as modified
        if target not in self.modified_nodes:
            self.modified_nodes += [target]

    ext_data: dict[Path, list[str]] = {}
    for target in visitor.extension_data:
        folder = Path(target.filename).parent
        # Generate the current source dict
        src_list: dict[str, BaseNode] = {}
        if isinstance(target.value, DictNode):
            src_list.update({k.value: v for k, v in target.value.args.kwargs.items()})
        ext_data.setdefault(folder, [])
        ext_data[folder] += src_list.keys()

    for target in visitor.extension_data:
        if target.var_name.value != "extension_data":
            continue
        folder = Path(target.filename).parent
        src_list = ext_data[folder]

        cython_files = sorted(list(folder.glob("*.pyx")))
        # Some cython files are compiled in a special way, so we don't want to add them
        special_cython_files = {
            "bliss.pyx",
            "mcqd.pyx",
            "tdlib.pyx",
        }
        cython_files = [x for x in cython_files if x.name not in special_cython_files]
        # Add all cython files that are not in the source list
        for file in cython_files:
            file_name = file.stem
            if file_name in src_list:
                continue
            token = Token("string", target.filename, 0, 0, 0, None, file_name)
            arg = ArgumentNode(Token("", target.filename, 0, 0, 0, None, "[]"))
            arg.append(
                StringNode(Token("string", target.filename, 0, 0, 0, None, file.name))
            )
            func = FunctionNode(_symbol("files"), _symbol("("), arg, _symbol(")"))
            target.value.args.kwargs.update({StringNode(token): func})
            if target not in self.modified_nodes:
                self.modified_nodes += [target]


def apply_changes(self: Rewriter):
    assert all(
        hasattr(x, "lineno") and hasattr(x, "colno") and hasattr(x, "filename")
        for x in self.modified_nodes
    )
    assert all(
        hasattr(x, "lineno") and hasattr(x, "colno") and hasattr(x, "filename")
        for x in self.to_remove_nodes
    )
    assert all(
        isinstance(x, (ArrayNode, FunctionNode, MethodNode, AssignmentNode))
        for x in self.modified_nodes
    )
    assert all(
        isinstance(x, (ArrayNode, AssignmentNode, FunctionNode))
        for x in self.to_remove_nodes
    )
    # Sort based on line and column in reversed order
    work_nodes = [{"node": x, "action": "modify"} for x in self.modified_nodes]
    work_nodes += [{"node": x, "action": "rm"} for x in self.to_remove_nodes]
    work_nodes = sorted(
        work_nodes, key=lambda x: (x["node"].lineno, x["node"].colno), reverse=True
    )
    work_nodes += [{"node": x, "action": "add"} for x in self.to_add_nodes]

    # Generating the new replacement string
    str_list = []
    for i in work_nodes:
        new_data = ""
        if i["action"] == "modify" or i["action"] == "add":
            printer = AstPrinter()
            i["node"].accept(printer)
            printer.post_process()
            new_data = printer.result.strip()
        data = {
            "file": i["node"].filename,
            "str": new_data,
            "node": i["node"],
            "action": i["action"],
        }
        str_list += [data]

    # Load build files
    files = {}
    for i in str_list:
        if i["file"] in files:
            continue
        fpath = os.path.realpath(os.path.join(self.sourcedir, i["file"]))
        fdata = ""
        # Create an empty file if it does not exist
        if not os.path.exists(fpath):
            with open(fpath, "w", encoding="utf-8"):
                pass
        with open(fpath, encoding="utf-8") as fp:
            fdata = fp.read()

        # Generate line offsets numbers
        m_lines = fdata.splitlines(True)
        offset = 0
        line_offsets = []
        for j in m_lines:
            line_offsets += [offset]
            offset += len(j)

        files[i["file"]] = {"path": fpath, "raw": fdata, "offsets": line_offsets}

    # Replace in source code
    def remove_node(i):
        offsets = files[i["file"]]["offsets"]
        raw = files[i["file"]]["raw"]
        node = i["node"]
        line = node.lineno - 1
        col = node.colno
        if isinstance(node, MethodNode):
            # The new data contains the source object as well
            col = node.source_object.colno
        elif isinstance(node, AssignmentNode):
            col = node.var_name.colno
        start = offsets[line] + col
        end = start
        if isinstance(node, (ArrayNode, FunctionNode, MethodNode)):
            end = offsets[node.end_lineno - 1] + node.end_colno
        elif isinstance(node, AssignmentNode):
            end = offsets[node.value.end_lineno - 1] + node.value.end_colno

        # Only removal is supported for assignments
        elif isinstance(node, AssignmentNode) and i["action"] == "rm":
            if isinstance(node.value, (ArrayNode, FunctionNode, MethodNode)):
                remove_node(
                    {"file": i["file"], "str": "", "node": node.value, "action": "rm"}
                )
                raw = files[i["file"]]["raw"]
            while raw[end] != "=":
                end += 1
            end += 1  # Handle the '='
            while raw[end] in {" ", "\n", "\t"}:
                end += 1

        files[i["file"]]["raw"] = raw[:start] + i["str"] + raw[end:]

    for i in str_list:
        if i["action"] in {"modify", "rm"}:
            remove_node(i)
        elif i["action"] == "add":
            files[i["file"]]["raw"] += i["str"] + "\n"

    # Write the files back
    for key, val in files.items():
        mlog.log("Rewriting", mlog.yellow(key))
        with open(val["path"], "w", encoding="utf-8") as fp:
            fp.write(val["raw"])


# Monkey patch the apply_changes method until https://github.com/mesonbuild/meson/pull/12899 is merged
Rewriter.apply_changes = apply_changes
# Monkey patch the update_python_sources method until this is upstreamed
Rewriter.process_update_python_sources = update_python_sources

rewriter = Rewriter(options.sourcedir)
visitor = AstPython()
rewriter.interpreter.visitors += [visitor]
rewriter.analyze_meson()
rewriter.process_update_python_sources(visitor)
rewriter.apply_changes()
rewriter.print_info()

# Run meson format
meson_format(
    Namespace(
        sources=[options.sourcedir],
        inplace=True,
        recursive=True,
        output=None,
        configuration=None,
        editor_config=None,
    )
)
