"""
Live code blocks

With ``SAGE_LIVE_DOC=yes``, this extension turns the examples of the
documentation into code blocks that the reader can execute, using
jupyter-sphinx.  Otherwise it makes the jupyter-sphinx directives no-ops.
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

from docutils import nodes
from sphinx.transforms import SphinxTransform
from sphinx.util import logging as sphinx_logging
from sphinx.util.docutils import SphinxDirective

SAGE_LIVE_DOC = os.environ.get('SAGE_LIVE_DOC', 'no')

logger = sphinx_logging.getLogger(__name__)


class SagecodeTransform(SphinxTransform):
    """
    Transform a code block to a live code block enabled by jupyter-sphinx.

    Effectively a code block like::

        EXAMPLE::

            sage: 1 + 1
            2

    is transformed into::

        EXAMPLE::

            sage: 1 + 1
            2

        .. ONLY:: html

            .. JUPYTER-EXECUTE::
                :hide-code:
                :hide-output:
                :raises:
                :stderr:

                1 + 1

    enabling live execution of the code.
    """
    # lower than the priority of jupyer_sphinx.execute.ExecuteJupyterCells
    default_priority = 170

    def apply(self):
        if self.app.builder.tags.has('html') or self.app.builder.tags.has('inventory'):
            for node in list(self.document.findall(nodes.literal_block)):
                if node.get('language') is None and node.astext().startswith('sage:'):
                    from docutils.nodes import Text
                    from docutils.nodes import container as Container
                    from docutils.nodes import label as Label
                    from docutils.nodes import literal_block as LiteralBlock
                    from sphinx_inline_tabs._impl import TabContainer
                    parent = node.parent
                    index = parent.index(node)
                    prev_node = node.previous_sibling()
                    if isinstance(prev_node, TabContainer):
                        # Make sure not to merge inline tabs for adjacent literal blocks
                        parent.insert(index, nodes.paragraph())
                        prev_node = parent[index]
                        index += 1
                    parent.remove(node)
                    # Tab for Sage code
                    container = TabContainer("", type="tab", new_set=False)
                    textnodes = [Text('Sage')]
                    label = Label("", "", *textnodes)
                    container += label
                    content = Container("", is_div=True, classes=["tab-content"])
                    content += node
                    container += content
                    parent.insert(index, container)
                    index += 1
                    if isinstance(prev_node, nodes.paragraph):
                        prev_node['classes'].append('with-sage-tab')

                    # Tab for preparsed version
                    from sage.repl.preparse import preparse
                    container = TabContainer("", type="tab", new_set=False)
                    textnodes = [Text('Python')]
                    label = Label("", "", *textnodes)
                    container += label
                    content = Container("", is_div=True, classes=["tab-content"])
                    example_lines = []
                    preparsed_lines = ['>>> from sage.all import *']
                    for line in node.rawsource.splitlines() + ['']:  # one extra to process last example
                        newline = line.lstrip()
                        if newline.startswith('....: '):
                            example_lines.append(newline[6:])
                        else:
                            if example_lines:
                                preparsed_example = preparse('\n'.join(example_lines))
                                prompt = '>>> '
                                for preparsed_line in preparsed_example.splitlines():
                                    preparsed_lines.append(prompt + preparsed_line)
                                    prompt = '... '
                                example_lines = []
                            if newline.startswith('sage: '):
                                example_lines.append(newline[6:])
                            else:
                                preparsed_lines.append(line)
                    preparsed = '\n'.join(preparsed_lines)
                    preparsed_node = LiteralBlock(preparsed, preparsed, language='ipycon')
                    content += preparsed_node
                    container += content
                    parent.insert(index, container)
                    index += 1
                    if isinstance(prev_node, nodes.paragraph):
                        prev_node['classes'].append('with-python-tab')

                    if SAGE_LIVE_DOC == 'yes':
                        # Tab for Jupyter-sphinx cell
                        from jupyter_sphinx.ast import CellInputNode, JupyterCellNode
                        source = node.rawsource
                        lines = []
                        for line in source.splitlines():
                            newline = line.lstrip()
                            if newline.startswith('sage: ') or newline.startswith('....: '):
                                lines.append(newline[6:])
                        cell_node = JupyterCellNode(
                                    execute=False,
                                    hide_code=False,
                                    hide_output=True,
                                    emphasize_lines=[],
                                    raises=False,
                                    stderr=True,
                                    code_below=False,
                                    classes=["jupyter_cell"])
                        cell_input = CellInputNode(classes=['cell_input','live-doc'])
                        cell_input += nodes.literal_block(
                            text='\n'.join(lines),
                            linenos=False,
                            linenostart=1)
                        cell_node += cell_input
                        container = TabContainer("", type="tab", new_set=False)
                        textnodes = [Text('Sage Live')]
                        label = Label("", "", *textnodes)
                        container += label
                        content = Container("", is_div=True, classes=["tab-content"])
                        content += cell_node
                        container += content
                        parent.insert(index, container)
                        index += 1
                        if isinstance(prev_node, nodes.paragraph):
                            prev_node['classes'].append('with-sage-live-tab')


class Ignore(SphinxDirective):

    has_content = True

    def run(self):
        return []


#: The directives of jupyter-sphinx, which execute the code they contain.
_JUPYTER_DIRECTIVES = (
    'jupyter-execute',
    'jupyter-kernel',
    'jupyter-input',
    'jupyter-output',
    'thebe-button',
)


def ignore_jupyter_directives(app):
    """
    Make the directives of jupyter-sphinx do nothing.

    Sphinx sets an extension up in the order of ``extensions``, where
    jupyter-sphinx comes after this one, so registering the directives here
    would leave the ones of jupyter-sphinx in place; they are registered once
    every extension is set up instead.  Executing the examples then needs a
    Sage kernel, which a documentation build has no reason to require.
    """
    for name in _JUPYTER_DIRECTIVES:
        app.add_directive(name, Ignore, override=True)


def setup(app):
    """
    Register this extension with Sphinx.
    """
    app.add_transform(SagecodeTransform)
    if SAGE_LIVE_DOC != 'yes':
        app.connect('builder-inited', ignore_jupyter_directives)
    return {'parallel_read_safe': True}
