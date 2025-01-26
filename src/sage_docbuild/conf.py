r"""
Sphinx build configuration

This file contains configuration needed to customize Sphinx input and output
behavior.
"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Load configuration shared with sage.misc.sphinxify
from sage.misc.sagedoc_conf import *

import sys
import os
import re
import importlib
import dateutil.parser
import sphinx
import sphinx.ext.intersphinx as intersphinx
from sphinx import highlighting
from sphinx.transforms import SphinxTransform
from sphinx.util.docutils import SphinxDirective
from IPython.lib.lexers import IPythonConsoleLexer, IPyLexer
from sage.misc.sagedoc import extlinks
from sage.env import SAGE_DOC_SRC, SAGE_DOC, PPLPY_DOCS, MATHJAX_DIR
from sage.misc.latex_macros import sage_mathjax_macros
from sage.features.sphinx import JupyterSphinx
from sage.features.all import all_features
import sage.version

# ---------------------
# General configuration
# ---------------------

SAGE_LIVE_DOC = os.environ.get('SAGE_LIVE_DOC', 'no')
SAGE_PREPARSED_DOC = os.environ.get('SAGE_PREPARSED_DOC', 'yes')

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sage_docbuild.ext.inventory_builder',
    'sage_docbuild.ext.multidocs',
    'sage_docbuild.ext.sage_autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.extlinks',
    'sphinx.ext.mathjax',
    'sphinx.ext.linkcode',
    'sphinx_copybutton',
    'sphinx_inline_tabs',
    'IPython.sphinxext.ipython_directive',
    'matplotlib.sphinxext.plot_directive',
]

if JupyterSphinx().is_present():
    extensions.append('jupyter_sphinx')

jupyter_execute_default_kernel = 'sagemath'

if SAGE_LIVE_DOC == 'yes':
    JupyterSphinx().require()
    SAGE_JUPYTER_SERVER = os.environ.get('SAGE_JUPYTER_SERVER', 'binder')
    if SAGE_JUPYTER_SERVER.startswith('binder'):
        # format: "binder" or
        #         "binder:sagemath/sage-binder-env" or
        #         "binder:sagemath/sage-binder-env/dev"
        if SAGE_JUPYTER_SERVER == 'binder':
            binder_repo = "sagemath/sage-binder-env/master"
        else:
            binder_repo = SAGE_JUPYTER_SERVER[7:]
        s = binder_repo.split('/', 2)
        if len(s) > 2:
            binder_options = {
                'repo': s[0] + '/' + s[1],
                'ref': s[2]
            }
        else:
            binder_options = {
                'repo': binder_repo
            }
        jupyter_sphinx_thebelab_config = {
            'requestKernel': False,
            'binderOptions': binder_options,
            'kernelOptions': {
                'name': "sagemath",
                'kernelName': "sagemath",
                'path': ".",
            },
            'selector': "div.live-doc"
        }
    else:  # local jupyter server
        SAGE_JUPYTER_SERVER_TOKEN = os.environ.get('SAGE_JUPYTER_SERVER_TOKEN', 'secret')
        jupyter_sphinx_thebelab_config = {
            'requestKernel': False,
            'kernelOptions': {
                'name': "sagemath",
                'kernelName': "sagemath",
                'path': ".",
                'serverSettings': {
                    'baseUrl': SAGE_JUPYTER_SERVER,
                    'token': SAGE_JUPYTER_SERVER_TOKEN
                },
            },
            'selector': "div.live-doc"
        }
    jupyter_sphinx_thebelab_config.update({
        'codeMirrorConfig': {
            'lineNumbers': True,
        }
    })

# This code is executed before each ".. PLOT::" directive in the Sphinx
# documentation. It defines a 'sphinx_plot' function that displays a Sage object
# through matplotlib, so that it will be displayed in the HTML doc
plot_html_show_source_link = False
plot_pre_code = r"""
# Set locale to prevent having commas in decimal numbers
# in tachyon input (see https://github.com/sagemath/sage/issues/28971)
import locale
locale.setlocale(locale.LC_NUMERIC, 'C')
def sphinx_plot(graphics, **kwds):
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from sage.misc.temporary_file import tmp_filename
    from sage.plot.graphics import _parse_figsize
    if os.environ.get('SAGE_SKIP_PLOT_DIRECTIVE', 'no') != 'yes':
        ## Option handling is taken from Graphics.save
        options = dict()
        if isinstance(graphics, sage.plot.graphics.Graphics):
            options.update(sage.plot.graphics.Graphics.SHOW_OPTIONS)
            options.update(graphics._extra_kwds)
            options.update(kwds)
        elif isinstance(graphics, sage.plot.multigraphics.MultiGraphics):
            options.update(kwds)
        else:
            graphics = graphics.plot(**kwds)
        dpi = options.pop('dpi', None)
        transparent = options.pop('transparent', None)
        fig_tight = options.pop('fig_tight', None)
        figsize = options.pop('figsize', None)
        if figsize is not None:
            figsize = _parse_figsize(figsize)
        plt.figure(figsize=figsize)
        figure = plt.gcf()
        if isinstance(graphics, (sage.plot.graphics.Graphics,
                                 sage.plot.multigraphics.MultiGraphics)):
            graphics.matplotlib(figure=figure, figsize=figsize, **options)
            if isinstance(graphics, (sage.plot.graphics.Graphics,
                                     sage.plot.multigraphics.GraphicsArray)):
                # for Graphics and GraphicsArray, tight_layout adjusts the
                # *subplot* parameters so ticks aren't cut off, etc.
                figure.tight_layout()
        else:
            # 3d graphics via png
            import matplotlib as mpl
            mpl.rcParams['image.interpolation'] = 'bilinear'
            mpl.rcParams['image.resample'] = False
            mpl.rcParams['figure.figsize'] = [8.0, 6.0]
            mpl.rcParams['figure.dpi'] = 80
            mpl.rcParams['savefig.dpi'] = 100
            fn = tmp_filename(ext=".png")
            graphics.save(fn)
            img = mpimg.imread(fn)
            plt.imshow(img)
            plt.axis("off")
        plt.margins(0)
        if not isinstance(graphics, sage.plot.multigraphics.MultiGraphics):
            plt.tight_layout(pad=0)

from sage.all_cmdline import *
"""

plot_html_show_formats = False
plot_formats = ['svg', 'pdf', 'png']

# We do *not* fully initialize intersphinx since we call it by hand
# in find_sage_dangling_links.
#, 'sphinx.ext.intersphinx']

# Add any paths that contain templates here, relative to this directory.
templates_path = [os.path.join(SAGE_DOC_SRC, 'common', 'templates'), 'templates']

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = ""
copyright = "2005--{}, The Sage Development Team".format(dateutil.parser.parse(sage.version.date).year)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
version = sage.version.version
release = sage.version.version

source_repository = 'https://github.com/sagemath/sage/'
source_branch = 'develop'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
# language = None

# The LaTeX engine to build the docs.
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-latex_engine
latex_engine = 'lualatex'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = '%B %d, %Y'

# List of glob-style patterns that should be excluded when looking for
# source files. [1] They are matched against the source file names
# relative to the source directory, using slashes as directory
# separators on all platforms.
exclude_patterns = ['.build']

# If true, '()' will be appended to :func: etc. cross-reference text.
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = True

# Default lexer to use when highlighting code blocks, using the IPython
# console lexers. 'ipycon' is the IPython console, which is what we want
# for most code blocks: anything with "sage:" prompts. For other IPython,
# like blocks which might appear in a notebook cell, use 'ipython'.
highlighting.lexers['ipycon'] = IPythonConsoleLexer(in1_regex=r'(sage:|>>>)', in2_regex=r'([.][.][.][.]:|[.][.][.])')
highlighting.lexers['ipython'] = IPyLexer()
highlight_language = 'ipycon'

# Create table of contents entries for domain objects (e.g. functions, classes,
# attributes, etc.). Default is True.
toc_object_entries = True

# A string that determines how domain objects (e.g. functions, classes,
# attributes, etc.) are displayed in their table of contents entry.
#
# Use "domain" to allow the domain to determine the appropriate number of parents
# to show. For example, the Python domain would show Class.method() and
# function(), leaving out the module. level of parents. This is the default
# setting.
#
# Use "hide" to only show the name of the element without any parents (i.e. method()).
#
# Use "all" to show the fully-qualified name for the object (i.e. module.Class.method()),
# displaying all parents.
toc_object_entries_show_parents = 'hide'

# -----------------------
# Extension configuration
# -----------------------

# include the todos
todo_include_todos = True

#
# intersphinx: Cross-links to other projects' online or installed documentation.
#
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
    'ipywidgets': ['https://ipywidgets.readthedocs.io/en/stable/'],
    'matplotlib': ['https://matplotlib.org/stable/'],
    'mpmath':     ['https://mpmath.org/doc/current/'],
    'networkx':   ['https://networkx.org/documentation/stable/'],
    'numpy':      ['https://numpy.org/doc/stable/'],
    'pplpy':      [PPLPY_DOCS, 'https://www.sagemath.org/pplpy/'],
    'python':     ['https://docs.python.org/'],
    'rpy2':       ['https://rpy2.github.io/doc/latest/html/'],
    'scipy':      ['https://docs.scipy.org/doc/scipy/'],
    'sympy':      ['https://docs.sympy.org/latest/'],
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
    if app.config.multidoc_first_pass == 1 or \
            not (os.path.exists(refpath) and os.path.exists(invpath)):
        return

    app.config.intersphinx_mapping = {key: _intersphinx_mapping(key)
                                      for key in _intersphinx_targets}

    # Add master intersphinx mapping
    dst = os.path.join(invpath, 'objects.inv')
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

    intersphinx.normalize_intersphinx_mapping(app, config)


# By default document is master.
multidocs_is_master = True

# https://sphinx-copybutton.readthedocs.io/en/latest/use.html
copybutton_prompt_text = r"sage: |[.][.][.][.]: |>>> |[.][.][.] |\$ "
copybutton_prompt_is_regexp = True
copybutton_exclude = '.linenos, .c1'  # exclude single comments (in particular, # optional!)
copybutton_only_copy_prompt_lines = True


# https://www.sphinx-doc.org/en/master/usage/extensions/linkcode.html
def linkcode_resolve(domain, info):
    import inspect
    from urllib.parse import quote
    from sage.misc.sageinspect import sage_getsourcelines
    if domain != 'py':
        return None
    if info['module']:
        m = importlib.import_module(info['module'])
        filename = quote(info['module'].replace('.', '/'))
        if m.__file__.endswith('py'):
            filename += '.py'
        else:
            filename += '.pyx'
        if 'fullname' in info:
            fullname = info['fullname']
            obj = m
            try:
                for attr in fullname.split('.'):
                    obj = getattr(obj, attr)
                lineno = sage_getsourcelines(obj)[-1]
            except Exception:  # catch all
                return None
            anchor = f'#L{lineno}'
        else:
            anchor = ''
        return f"{source_repository}blob/develop/src/{filename}{anchor}"
    return None


# -----------------------
# Options for HTML output
# -----------------------

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = [os.path.join(SAGE_DOC_SRC, "common", "themes")]

# Deprecated Sage classic theme:
#
#   html_theme = "sage-classic"
#   html_theme_options = {}
#
# See the directory doc/common/themes/sage-classic/ for theme files.

# Sphinx theme "furo" does not permit an extension. Do not attempt to make
# a "sage-furo" theme.
html_theme = "furo"

# Theme options are theme-specific and customize the look and feel of
# a theme further.  For a list of options available for each theme,
# see the documentation.
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#0f0fff",
        "color-brand-content": "#0f0fff",
    },
    "light_logo": "logo_sagemath_black.svg",
    "dark_logo": "logo_sagemath_white.svg",
    # Furo can add a small edit button to each document to allow visitors to
    # easily propose changes to that document using the repository’s source
    # control system.
    # https://pradyunsg.me/furo/customisation/edit-button/#adding-an-edit-button
    "source_repository": source_repository,
    "source_branch": source_branch,
    # "source_directory" is defined in conf.py customized for the doc
}

# Check the condition for announcement banner
github_ref = os.environ.get('GITHUB_REF', '')
if github_ref:
    match = re.search(r'refs/pull/(\d+)/merge', github_ref)
    if match:
        pr_number = match.group(1)
is_for_develop = github_ref.startswith('refs/heads/develop')
is_for_github_pr = github_ref and match and pr_number
is_stable_release = version.split('.')[-1].isnumeric()

if is_for_develop or is_for_github_pr or not is_stable_release:  # condition for announcement banner
    # This URL is hardcoded in the file .github/workflows/doc-publish.yml.
    # See NETLIFY_ALIAS of the "Deploy to Netlify" step.
    ver = f'<a href="https://doc-develop--sagemath.netlify.app/html/en/index.html">{version}</a>'
    if is_for_github_pr:
        pr_url = f'https://github.com/sagemath/sage/pull/{pr_number}'
        pr_sha = os.environ.get('PR_SHA', '')
        pr_commit = pr_url + f'/commits/{pr_sha}'
        ver += f' built with GitHub PR <a href="{pr_url}">#{pr_number}</a>' \
               f' on <a href="{pr_commit}">{pr_sha[:7]}</a>' \
               f' [<a href="/changes.html">changes</a>]'
    banner = f'This is documentation for Sage version {ver} for development purpose.'
    html_theme_options.update({ "announcement": banner })

# The name of the Pygments (syntax highlighting) style to use. This
# overrides a HTML theme's corresponding setting.
pygments_style = "sphinx"
pygments_dark_style = "monokai"

# Add siderbar/home.html to the default sidebar.
html_sidebars = {
    "**": [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/version-selector.html",
        "sidebar/search.html",
        "sidebar/home.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
        "sidebar/variant-selector.html",
    ]
}

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'custom-furo.css',
    'custom-jupyter-sphinx.css',
    'custom-codemirror-monokai.css',
    'custom-tabs.css',
]

html_js_files = [
    'jupyter-sphinx-furo.js',
]

# A list of paths that contain extra templates (or templates that overwrite
# builtin/theme-specific templates). Relative paths are taken as relative
# to the configuration directory.
templates_path = [os.path.join(SAGE_DOC_SRC, 'common', 'templates-furo')] + templates_path

# HTML style sheet. This overrides a HTML theme's corresponding setting.
# html_style = 'default.css'

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (within the static path) to place at the top of
# the sidebar.
# html_logo = 'sagelogo-word.ico'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = 'favicon.ico'

# html_static_path defined here and imported in the actual configuration file
# conf.py read by Sphinx was the cause of subtle bugs in builders (see #30418 for
# instance). Hence now html_common_static_path contains the common paths to static
# files, and is combined to html_static_path in each conf.py file read by Sphinx.
html_common_static_path = [os.path.join(SAGE_DOC_SRC, 'common', 'static'), 'static']

# Configure MathJax
# https://docs.mathjax.org/en/latest/options/input/tex.html
mathjax3_config = {
    "tex": {
        # Add custom sage macros
        # http://docs.mathjax.org/en/latest/input/tex/macros.html
        "macros": sage_mathjax_macros(),
        # Add $...$ as possible inline math
        # https://docs.mathjax.org/en/latest/input/tex/delimiters.html#tex-and-latex-math-delimiters
        "inlineMath": [["$", "$"], ["\\(", "\\)"]],
        # Increase the limit the size of the string to be processed
        # https://docs.mathjax.org/en/latest/options/input/tex.html#option-descriptions
        "maxBuffer": 50 * 1024,
        # Use colorv2 extension instead of built-in color extension
        # https://docs.mathjax.org/en/latest/input/tex/extensions/autoload.html#tex-autoload-options
        # https://docs.mathjax.org/en/latest/input/tex/extensions/colorv2.html#tex-colorv2
        "autoload": {"color": [], "colorv2": ["color"]},
    },
}

if os.environ.get('SAGE_USE_CDNS', 'no') == 'yes':
    mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"
else:
    mathjax_path = os.path.join(MATHJAX_DIR, 'tex-chtml.js')

# A list of glob-style patterns that should be excluded when looking for source
# files. They are matched against the source file names relative to the
# source directory, using slashes as directory separators on all platforms.
exclude_patterns = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_use_modindex = True

# A list of prefixes that are ignored for sorting the Python module index ( if
# this is set to ['foo.'], then foo.bar is shown under B, not F). Works only
# for the HTML builder currently.
modindex_common_prefix = ['sage.']

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = True

# If true, the reST sources are included in the HTML build as _sources/<name>.
# html_copy_source = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = ''

# Output file base name for HTML help builder.
# htmlhelp_basename = ''

# ------------------------
# Options for LaTeX output
# ------------------------

# See http://sphinx-doc.org/config.html#confval-latex_elements
latex_elements = {}

# The paper size ('letterpaper' or 'a4paper').
#latex_elements['papersize'] = 'letterpaper'

# The font size ('10pt', '11pt' or '12pt').
#latex_elements['pointsize'] = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = []

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = 'sagelogo-word.png'

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# Additional stuff for the LaTeX preamble.
latex_elements['preamble'] = r"""
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{textcomp}
\usepackage{mathrsfs}
\usepackage{iftex}

\let\textLaTeX\LaTeX
\AtBeginDocument{\renewcommand*{\LaTeX}{\hbox{\textLaTeX}}}

% Workaround for a LaTeX bug -- see Issue #31397 and
% https://tex.stackexchange.com/questions/583391/mactex-2020-error-with-report-hyperref-mathbf-in-chapter.
\makeatletter
\pdfstringdefDisableCommands{%
  \let\mathbf\@firstofone
}
\makeatother
"""

# Enable "hard wrapping" long code lines (only applies if breaking
# long codelines at spaces or other suitable places failed, typically
# this is for long decimal expansions or possibly long string identifiers)
latex_elements['sphinxsetup'] = "verbatimforcewraps=true"

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_use_modindex = True

# -------------------------
# add LaTeX macros for Sage
# -------------------------

from sage.misc.latex_macros import sage_latex_macros

try:
    pngmath_latex_preamble  # check whether this is already defined
except NameError:
    pngmath_latex_preamble = ""

for macro in sage_latex_macros():
    # used when building latex and pdf versions
    latex_elements['preamble'] += macro + '\n'
    # used when building html version
    pngmath_latex_preamble += macro + '\n'


# ------------------------------------------
# add custom context variables for templates
# ------------------------------------------

def add_page_context(app, pagename, templatename, context, doctree):
    # # The template function
    # def template_function(arg):
    #     return "Your string is " + arg
    # # Add it to the page's context
    # context['template_function'] = template_function
    path1 = os.path.dirname(app.builder.get_outfilename(pagename))
    path2 = os.path.join(SAGE_DOC, 'html', 'en')
    relpath = os.path.relpath(path2, path1)
    context['release'] = release
    context['documentation_title'] = f'Version {release} Documentation'
    context['documentation_root'] = os.path.join(relpath, 'index.html')
    if 'website' in path1:
        context['title'] = 'Documentation'
        context['website'] = True
        context['documentation_root'] = 'index.html'
    if 'reference' in path1 and not path1.endswith('reference'):
        path2 = os.path.join(SAGE_DOC, 'html', 'en', 'reference')
        relpath = os.path.relpath(path2, path1)
        context['reference_title'] = f'Version {release} Reference Manual'
        context['reference_root'] = os.path.join(relpath, 'index.html')
        context['refsub'] = True
        if pagename.startswith('sage/'):
            # This is for adding small view/edit buttons using Furo's feature:
            # https://pradyunsg.me/furo/customisation/top-of-page-buttons/
            # This works well if the source file is '.rst' file. But the '.rst'
            # files in the directory 'sage/' are generated by the Sphinx
            # autodoc from the Python or Cython source files. Hence we tweak
            # here template context variables so that links to the correct
            # source files are generated.
            suffix = '.py' if importlib.import_module(pagename.replace('/','.')).__file__.endswith('.py') else '.pyx'
            context['page_source_suffix'] = suffix
            context['theme_source_view_link'] = os.path.join(source_repository, f'blob/develop/src', '{filename}')
            context['theme_source_edit_link'] = os.path.join(source_repository, f'edit/develop/src', '{filename}')


dangling_debug = False


def debug_inf(app, message):
    if dangling_debug:
        app.info(message)


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


def find_sage_dangling_links(app, env, node, contnode):
    r"""
    Try to find dangling link in local module imports or all.py.
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

    # Workaround: in Python's doc 'object', 'list', ... are documented as a
    # function rather than a class
    if reftarget in base_class_as_func and reftype == 'class':
        node['reftype'] = 'func'

    res = call_intersphinx(app, env, node, contnode)
    if res:
        debug_inf(app, "++ DONE %s" % (res['refuri']))
        return res

    if node.get('refdomain') != 'py':  # not a python file
        return None

    try:
        module = node['py:module']
        cls = node['py:class']
    except KeyError:
        debug_inf(app, "-- no module or class for :%s:%s" % (reftype,
                                                             reftarget))
        return None

    basename = reftarget.split(".")[0]
    try:
        target_module = getattr(sys.modules['sage.all'], basename).__module__
        debug_inf(app, "++ found %s using sage.all in %s" % (basename, target_module))
    except AttributeError:
        try:
            target_module = getattr(sys.modules[node['py:module']], basename).__module__
            debug_inf(app, "++ found %s in this module" % (basename,))
        except AttributeError:
            debug_inf(app, "-- %s not found in sage.all or this module" % (basename))
            return None
        except KeyError:
            target_module = None
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
    elif len(matches) > 1:
        env.warn(target_module,
                 'more than one target found for cross-reference '
                 '%r: %s' % (newtarget,
                             ', '.join(match[0] for match in matches)),
                 node.line)
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


# lists of basic Python class which are documented as functions
base_class_as_func = [
    'bool', 'complex', 'dict', 'file', 'float',
    'frozenset', 'int', 'list', 'long', 'object',
    'set', 'slice', 'str', 'tuple', 'type', 'unicode', 'xrange']


# nitpicky option configuration: Put here broken links we want to ignore. For
# link to the Python documentation several links where broken because there
# where class listed as functions. Expand the list 'base_class_as_func' above
# instead of marking the link as broken.
nitpick_ignore = [
    ('py:class', 'twisted.web2.resource.Resource'),
    ('py:class', 'twisted.web2.resource.PostableResource')]


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
                app.warn('Pickling of nested class %r is probably broken. '
                         'Please set the metaclass of the parent class to '
                         'sage.misc.nested_class.NestedClassMetaclass.' % (
                        v.__module__ + '.' + name + '.' + nm))


def skip_member(app, what, name, obj, skip, options):
    """
    To suppress Sphinx warnings / errors, we

    - Don't include [aliases of] builtins.

    - Don't include the docstring for any nested class which has been
      inserted into its module by
      :class:`sage.misc.NestedClassMetaclass` only for pickling.  The
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
                    from docutils.nodes import container as Container, label as Label, literal_block as LiteralBlock, Text
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

                    if SAGE_PREPARSED_DOC == 'yes':
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
                        from jupyter_sphinx.ast import JupyterCellNode, CellInputNode
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


# This replaces the setup() in sage.misc.sagedoc_conf
def setup(app):
    app.connect('autodoc-process-docstring', process_docstring_cython)
    app.connect('autodoc-process-docstring', process_directives)
    app.connect('autodoc-process-docstring', process_docstring_module_title)
    app.connect('autodoc-process-docstring', process_dollars)
    app.connect('autodoc-process-docstring', process_inherited)
    if os.environ.get('SAGE_SKIP_TESTS_BLOCKS', False):
        app.connect('autodoc-process-docstring', skip_TESTS_block)
    app.connect('autodoc-skip-member', skip_member)
    app.add_transform(SagemathTransform)
    if SAGE_LIVE_DOC == 'yes' or SAGE_PREPARSED_DOC == 'yes':
        app.add_transform(SagecodeTransform)
    if not JupyterSphinx().is_present():
        app.add_directive("jupyter-execute", Ignore)
        app.add_directive("jupyter-kernel", Ignore)
        app.add_directive("jupyter-input", Ignore)
        app.add_directive("jupyter-output", Ignore)
        app.add_directive("thebe-button", Ignore)

    # When building the standard docs, app.srcdir is set to SAGE_DOC_SRC +
    # 'LANGUAGE/DOCNAME'.
    if app.srcdir.is_relative_to(SAGE_DOC_SRC):
        app.add_config_value('intersphinx_mapping', {}, False)
        app.add_config_value('intersphinx_cache_limit', 5, False)
        app.add_config_value('intersphinx_disabled_reftypes', [], False)
        app.add_config_value('intersphinx_timeout', None, False)
        app.connect('config-inited', set_intersphinx_mappings)
        app.connect('builder-inited', intersphinx.load_mappings)
        # We do *not* fully initialize intersphinx since we call it by hand
        # in find_sage_dangling_links.
        #   app.connect('missing-reference', missing_reference)
        app.connect('missing-reference', find_sage_dangling_links)
        app.connect('html-page-context', add_page_context)


# Conditional content
# https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#tags
# https://www.sphinx-doc.org/en/master/usage/configuration.html#conf-tags
# https://github.com/readthedocs/readthedocs.org/issues/4603#issuecomment-1411594800
def feature_tags():
    for feature in all_features():
        if feature.is_present():
            yield 'feature_' + feature.name.replace('.', '_')
