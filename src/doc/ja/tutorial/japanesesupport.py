# -*- coding: utf-8 -*-
import re
__RGX = re.compile(r'([^!-~])[\n\r\t]+([^!-~])')

def trunc_whitespace(app, doctree, docname):
    from docutils.nodes import Text, paragraph
    if not app.config.japanesesupport_trunc_whitespace:
        return
    for node in doctree.traverse(Text):
        if isinstance(node.parent, paragraph):
            newtext = node.astext()
            #↓「非ASCII」+「"\n\r\t"たち」+「非ASCII」
            # の場合だけ置換する…
            newtext = __RGX.sub(r"\1\2", newtext)
            #newtext = newtext.strip()
            node.parent.replace(node, Text(newtext))

def setup(app):
    app.add_config_value('japanesesupport_trunc_whitespace', True, True)
    app.connect("doctree-resolved", trunc_whitespace)
