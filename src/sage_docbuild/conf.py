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

import importlib
import os
import re

import dateutil.parser
from IPython.lib.lexers import IPyLexer, IPythonConsoleLexer
from sphinx import highlighting
from sphinx.util import logging as sphinx_logging

import sage.version
from sage.env import MATHJAX_DIR, SAGE_DOC, SAGE_DOC_SRC
from sage.features.sphinx import JupyterSphinx
from sage.misc.latex_macros import sage_mathjax_macros
from sage.misc.sagedoc import extlinks as extlinks  # noqa: PLC0414
from sage.misc.sagedoc_conf import *  # Load configuration shared with sage.misc.sphinxify

logger = sphinx_logging.getLogger(__name__)

# ---------------------
# General configuration
# ---------------------

SAGE_LIVE_DOC = os.environ.get('SAGE_LIVE_DOC', 'no')

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sage_docbuild.ext.inventory_builder',
    'sage_docbuild.ext.multidocs',
    'sage_docbuild.ext.sage_autodoc',
    'sage_docbuild.ext.crossrefs',
    'sage_docbuild.ext.livedoc',
    'sage_docbuild.ext.members',
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

# By default document is master.
multidocs_is_master = True

# https://sphinx-copybutton.readthedocs.io/en/latest/use.html
copybutton_prompt_text = r"sage: |[.][.][.][.]: |>>> |[.][.][.] |\$ "
copybutton_line_continuation_character = "\\"
copybutton_prompt_is_regexp = True
copybutton_exclude = '.linenos, .c1'  # exclude single comments (in particular, # optional!)
copybutton_only_copy_prompt_lines = True


# https://www.sphinx-doc.org/en/master/usage/extensions/linkcode.html
def linkcode_resolve(domain, info):
    from urllib.parse import quote

    from sage.misc.sageinspect import sage_getsourcelines
    if domain != 'py':
        return None
    if info['module']:
        m = importlib.import_module(info['module'])
        source = getattr(m, '__file__', None)
        if not source:
            # A namespace package has no source file to link to.
            return None
        path = info['module'].replace('.', '/')
        if os.path.basename(source).startswith('__init__.'):
            # The module is a package, and its name points at the directory
            # holding the initializer that defines it.
            path += '/__init__'
        filename = quote(path)
        if source.endswith('py'):
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
            context['theme_source_view_link'] = os.path.join(source_repository, 'blob/develop/src', '{filename}')
            context['theme_source_edit_link'] = os.path.join(source_repository, 'edit/develop/src', '{filename}')


# ---------------------------------------
# Sub-documents of the reference manual
# ---------------------------------------

def reference_subdocument(directory=None):
    r"""
    Return the configuration values proper to one sub-document of the
    reference manual.

    Every sub-document is a Sphinx project of its own, with a configuration
    file that takes the shared configuration from this module and then the
    values returned here, which only depend on where the sub-document lives::

        from sage_docbuild.conf import *
        from sage_docbuild.conf import reference_subdocument

        globals().update(reference_subdocument())

    Anything proper to a single sub-document belongs in its configuration
    file, after that call.

    INPUT:

    - ``directory`` -- the directory of the sub-document; defaults to the
      current directory, which is the one holding the configuration file that
      Sphinx is reading

    As a side effect, the shared ``html_theme_options`` and ``latex_elements``
    are given the entries that depend on the sub-document; a Sphinx run builds
    a single document, so they cannot be shared.

    EXAMPLES::

        sage: import os, tempfile
        sage: from sage_docbuild.conf import reference_subdocument
        sage: def subdocument(name, index):
        ....:     directory = os.path.join(tempfile.mkdtemp(), name)
        ....:     os.mkdir(directory)
        ....:     with open(os.path.join(directory, 'index.rst'), 'w') as f:
        ....:         _ = f.write(index)
        ....:     return reference_subdocument(directory)

        sage: config = subdocument('algebras', 'Algebras\n========\n')
        sage: config['htmlhelp_basename']
        'algebras'
        sage: config['html_title']
        'Algebras'
        sage: config['latex_documents']
        [('index', 'algebras.tex', 'Algebras', 'The Sage Development Team', 'manual')]
        sage: config['multidocs_is_master']
        False

    A title in backticks is math, which the HTML title writes with dollars::

        sage: subdocument('padics', '`p`-adics\n=========\n')['html_title']
        '$p$-adics'

    Without a title, the name of the directory is used::

        sage: subdocument('padics', 'No title here.\n')['html_title']
        'Padics'
    """
    directory = os.path.abspath(directory or '.')
    name = os.path.basename(directory)

    # We use the main document's title, if we can find it.
    title = ''
    with open(os.path.join(directory, 'index.rst'), encoding='utf-8') as rst_file:
        rst_lines = rst_file.read().splitlines()
    for i, line in enumerate(rst_lines):
        if line.startswith('==') and i > 0:
            title = rst_lines[i - 1].strip()
            break
    # Otherwise, we use this directory's name.
    if not title:
        title = name.capitalize()
    title = title.replace('`', '$')

    # We use the directory's name to add small view/edit buttons.
    source = f'src/doc/en/reference/{name}'
    html_theme_options.update({
        'source_view_link': os.path.join(source_repository, 'blob/develop', source, '{filename}'),
        'source_edit_link': os.path.join(source_repository, 'edit/develop', source, '{filename}'),
    })

    latex_elements['hyperref'] = r"""
\usepackage{xr}
\externaldocument[../references/]{../references/references}
% Include hyperref last.
\usepackage{hyperref}
% Fix anchor placement for figures with captions.
\usepackage{hypcap}% it must be loaded after hyperref.
% Set up styles of URL: it should be placed after hyperref.
\urlstyle{same}"""

    return {
        # Paths that contain custom static files (such as style sheets). They
        # are copied after the builtin static files, so a file named
        # "default.css" will overwrite the builtin "default.css".
        'html_static_path': [] + html_common_static_path,
        'project': title,
        'html_title': title,
        'html_short_title': title,
        # Output file base name for HTML help builder.
        'htmlhelp_basename': name,
        # Grouping the document tree into LaTeX files: (source start file,
        # target name, title, author, document class [howto/manual]).
        'latex_documents': [
            ('index', name + '.tex', title, 'The Sage Development Team', 'manual')
        ],
        # Ignore all .rst in the _sage subdirectory
        'exclude_patterns': exclude_patterns + ['_sage'],
        'multidocs_is_master': False,
    }



autodoc_type_aliases = {
    'IntegerMod_abstract': 'sage.rings.finite_rings.integer_mod.IntegerMod_abstract',
    'EllipticCurve_finite_field': 'sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field',
    'EllipticCurvePoint_finite_field': 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field',
}


# nitpicky option configuration: Put here broken links we want to ignore.
# For links to the Python documentation, expand the role-fallback lists above
# instead of marking the link as broken.  For external projects, prefer adding
# a vendored inventory before removing entries from this list.  A link to an
# implementation module resolves through _public_alias(), and a name that no
# module defines usually means an annotation that Sphinx could not evaluate;
# see _type_checking_aliases() in sage_docbuild.ext.sage_autodoc.
nitpick_ignore = []





# This replaces the setup() in sage.misc.sagedoc_conf.  Everything that is not
# tied to a configuration value of this file lives in an extension module of
# sage_docbuild.ext, listed in ``extensions`` above.
def setup(app):
    app.connect('autodoc-process-docstring', process_docstring_cython)
    app.connect('autodoc-process-docstring', process_directives)
    app.connect('autodoc-process-docstring', process_docstring_module_title)
    app.connect('autodoc-process-docstring', process_dollars)
    app.connect('autodoc-process-docstring', process_inherited)
    app.connect('autodoc-process-docstring', process_docstring_aliases)
    if os.environ.get('SAGE_SKIP_TESTS_BLOCKS', False):
        app.connect('autodoc-process-docstring', skip_TESTS_block)
    app.add_transform(SagemathTransform)

    # When building the standard docs, app.srcdir is set to SAGE_DOC_SRC +
    # 'LANGUAGE/DOCNAME'.
    if app.srcdir.is_relative_to(SAGE_DOC_SRC):
        app.connect('html-page-context', add_page_context)
