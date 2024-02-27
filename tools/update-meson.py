import argparse
from pathlib import Path

from mesonbuild.ast import AstVisitor
from mesonbuild.ast.interpreter import MethodNode
from mesonbuild.mparser import AssignmentNode, BaseNode, DictNode, SymbolNode
from mesonbuild.rewriter import (
    ArgumentNode,
    ArrayNode,
    BaseStringNode,
    FunctionNode,
    Rewriter,
    StringNode,
    Token,
)

# Get target directory from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("sourcedir", help="Source directory")
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
    return SymbolNode(Token('', '', 0, 0, 0, (0, 0), val))


def add_python_sources(self: Rewriter, visitor: AstPython):
    for target in visitor.install_sources_calls:
        # Generate the current source list
        src_list: list[str] = []
        for arg in arg_list_from_node(target):
            if isinstance(arg, BaseStringNode):
                src_list += [arg.value]

        folder = Path(target.filename).parent
        python_files = sorted(
            list(folder.glob("*.py"))
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

        # Append to the AST at the right place
        target.args.arguments = sorted(
            target.args.arguments + to_append, key=lambda x: x.value
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


Rewriter.process_add_python_sources = add_python_sources
rewriter = Rewriter(options.sourcedir)
visitor = AstPython()
rewriter.interpreter.visitors += [visitor]
rewriter.analyze_meson()
rewriter.process_add_python_sources(visitor)
rewriter.apply_changes()
rewriter.print_info()
