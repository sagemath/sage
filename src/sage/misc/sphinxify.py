# -*- coding: utf-8 -*
r"""
Process docstrings with Sphinx

Processes docstrings with Sphinx. Can also be used as a commandline script:

``python sphinxify.py <text>``

AUTHORS:

- Tim Joseph Dumol (2009-09-29): initial version
"""

#*****************************************************************************
#       Copyright (C) 2009 Tim Dumol <tim@timdumol.com>
#
# Distributed under the terms of the BSD License
#*****************************************************************************
from __future__ import absolute_import, print_function

import os
import re
import shutil
from tempfile import mkdtemp
from sage.env import SAGE_DOC_SRC
from sphinx.application import Sphinx


def sphinxify(docstring, format='html'):
    r"""
    Runs Sphinx on a ``docstring``, and outputs the processed
    documentation.

    INPUT:

    - ``docstring`` -- string -- a ReST-formatted docstring

    - ``format`` -- string (optional, default 'html') -- either 'html' or
      'text'

    OUTPUT:

    - string -- Sphinx-processed documentation, in either HTML or
      plain text format, depending on the value of ``format``

    EXAMPLES::

        sage: from sage.misc.sphinxify import sphinxify
        sage: sphinxify('A test')
        '...<div class="docstring">\n    \n  <p>A test</p>\n\n\n</div>'
        sage: sphinxify('**Testing**\n`monospace`')
        '...<div class="docstring"...<strong>Testing</strong>\n<span class="math"...</p>\n\n\n</div>'
        sage: sphinxify('`x=y`')
        '...<div class="docstring">\n    \n  <p><span class="math">x=y</span></p>\n\n\n</div>'
        sage: sphinxify('`x=y`', format='text')
        'x=y\n'
        sage: sphinxify(':math:`x=y`', format='text')
        'x=y\n'

    TESTS::

        sage: n = len(sys.path)
        sage: _ = sphinxify('A test')
        sage: assert n == len(sys.path)
    """
    srcdir = mkdtemp()
    base_name = os.path.join(srcdir, 'docstring')
    rst_name = base_name + '.rst'

    if format == 'html':
        suffix = '.html'
    else:
        suffix = '.txt'
    output_name = base_name + suffix

    with open(rst_name, 'w') as filed:
        filed.write(docstring)

    # Sphinx constructor: Sphinx(srcdir, confdir, outdir, doctreedir,
    # buildername, confoverrides, status, warning, freshenv).
    confdir = os.path.join(SAGE_DOC_SRC, 'en', 'introspect')

    doctreedir = os.path.join(srcdir, 'doctrees')
    confoverrides = {'html_context': {}, 'master_doc': 'docstring'}

    import sys
    old_sys_path = list(sys.path)  # Sphinx modifies sys.path
    sphinx_app = Sphinx(srcdir, confdir, srcdir, doctreedir, format,
                        confoverrides, None, None, True)
    sphinx_app.build(None, [rst_name])
    sys.path = old_sys_path

    # We need to remove "_" from __builtin__ that the gettext module installs
    from six.moves import builtins
    builtins.__dict__.pop('_', None)

    if os.path.exists(output_name):
        output = open(output_name, 'r').read()
        output = output.replace('<pre>', '<pre class="literal-block">')

        # Translate URLs for media from something like
        #    "../../media/...path.../blah.png"
        # or
        #    "/media/...path.../blah.png"
        # to
        #    "/doc/static/reference/media/...path.../blah.png"
        output = re.sub("""src=['"](/?\.\.)*/?media/([^"']*)['"]""",
                          'src="/doc/static/reference/media/\\2"',
                          output)
        # Remove spurious \(, \), \[, \].
        output = output.replace('\\(', '').replace('\\)', '').replace('\\[', '').replace('\\]', '')
    else:
        print("BUG -- Sphinx error")
        if format == 'html':
            output = '<pre class="introspection">%s</pre>' % docstring
        else:
            output = docstring

    shutil.rmtree(srcdir, ignore_errors=True)

    return output


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        print(sphinxify(sys.argv[1]))
    else:
        print("""Usage:
%s 'docstring'

docstring -- docstring to be processed
""")
