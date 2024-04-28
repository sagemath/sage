# -*- coding: utf-8 -*-
"""
The Sage ZMQ Kernel

Version of the Jupyter kernel when running Sage inside the Jupyter
notebook or remote Jupyter sessions.
"""

# ***************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import sys
import warnings
with warnings.catch_warnings():
    # When upstream pydevd (as opposed to the bundled version) is used
    # with debugpy, a PEP 420 warning is emitted. Debugpy and/or
    # pydevd will eventually work around this, but as of September
    # 2023, hiding the warning gives us more flexibility in the
    # versions of those packages that we can accept.
    warnings.filterwarnings("ignore",
                            message=r".*pkg_resources\.declare_namespace",
                            category=DeprecationWarning)
    from ipykernel.ipkernel import IPythonKernel

from ipykernel.zmqshell import ZMQInteractiveShell
from traitlets import Type

from sage.env import SAGE_VERSION
from sage.repl.interpreter import SageNotebookInteractiveShell
from sage.repl.ipython_extension import SageJupyterCustomizations


class SageZMQInteractiveShell(SageNotebookInteractiveShell, ZMQInteractiveShell):
    pass


class SageKernel(IPythonKernel):
    implementation = 'sage'
    implementation_version = SAGE_VERSION

    shell_class = Type(SageZMQInteractiveShell)

    def __init__(self, **kwds):
        """
        The Sage Jupyter Kernel

        INPUT:

        See the Jupyter documentation

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.kernel import SageKernel
            sage: SageKernel.__new__(SageKernel)
            <sage.repl.ipython_kernel.kernel.SageKernel object at 0x...>
        """
        super().__init__(**kwds)
        SageJupyterCustomizations(self.shell)

    @property
    def banner(self):
        r"""
        The Sage Banner

        The value of this property is displayed in the Jupyter
        notebook.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.kernel import SageKernel
            sage: sk = SageKernel.__new__(SageKernel)
            sage: print(sk.banner)
            ┌...SageMath version...
        """
        from sage.misc.banner import banner_text
        return banner_text()

    @property
    def help_links(self):
        r"""
        Help in the Jupyter Notebook

        OUTPUT:

        See the Jupyter documentation.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.kernel import SageKernel
            sage: sk = SageKernel.__new__(SageKernel)
            sage: sk.help_links
            [{'text': 'Sage Documentation',
              'url': 'https://doc.sagemath.org/html/en/index.html'},
             ...]
        """
        # DEPRECATED: The URLs in the form 'kernelspecs/...' were used for
        # classical Jupyter notebooks. For instance,
        #
        #  'kernelspecs/sagemath/doc/html/en/index.html'
        #
        # is constructed by kernel_url('doc/html/en/index.html'), but these
        # URLs of local files don't work for JupyterLab. Hence all URLs here
        # have been replaced with URLs of online documents.

        from sage.repl.ipython_kernel.install import SageKernelSpec
        identifier = SageKernelSpec.identifier()

        def kernel_url(x):
            # URLs starting with 'kernelspecs' are prepended by the
            # browser with the appropriate path
            return 'kernelspecs/{0}/{1}'.format(identifier, x)

        return [
            {
                'text': 'Sage Documentation',
                'url': "https://doc.sagemath.org/html/en/index.html",
            },
            {
                'text': 'A Tour of Sage',
                'url': "https://doc.sagemath.org/html/en/a_tour_of_sage/index.html",
            },
            {
                'text': 'Tutorial',
                'url': "https://doc.sagemath.org/html/en/tutorial/index.html",
            },
            {
                'text': 'Thematic Tutorials',
                'url': "https://doc.sagemath.org/html/en/thematic_tutorials/index.html",
            },
            {
                'text': 'PREP Tutorials',
                'url': "https://doc.sagemath.org/html/en/prep/index.html",
            },
            {
                'text': 'Constructions',
                'url': "https://doc.sagemath.org/html/en/constructions/index.html",
            },
            {
                'text': 'FAQ',
                'url': "https://doc.sagemath.org/html/en/faq/index.html",
            },
            {
                'text': 'Reference',
                'url': "https://doc.sagemath.org/html/en/reference/index.html",
            },
            {
                'text': "Developer's Guide",
                'url': "https://doc.sagemath.org/html/en/developer/index.html",
            },
            {
                'text': "Python",
                'url': "http://docs.python.org/%i.%i" % sys.version_info[:2],
            },
            {
                'text': "IPython",
                'url': "http://ipython.org/documentation.html",
            },
            {
                'text': 'Singular',
                'url': 'http://www.singular.uni-kl.de/Manual/latest/index.htm',
            },
            {
                'text': 'GAP',
                'url': 'http://gap-system.org/Manuals/doc/ref/chap0.html',
            },
            {
                'text': "NumPy",
                'url': "http://docs.scipy.org/doc/numpy/reference/",
            },
            {
                'text': "SciPy",
                'url': "http://docs.scipy.org/doc/scipy/reference/",
            },
            {
                'text': "SymPy",
                'url': 'http://docs.sympy.org/latest/index.html',
            },
            {
                'text': "Matplotlib",
                'url': "https://matplotlib.org/contents.html",
            },
            {
                'text': "Markdown",
                'url': "http://help.github.com/articles/github-flavored-markdown",
            },
        ]

    def pre_handler_hook(self):
        """
        Restore the signal handlers to their default values at Sage
        startup, saving the old handler at the ``saved_sigint_handler``
        attribute. This is needed because Jupyter needs to change the
        ``SIGINT`` handler.

        See :issue:`19135`.

        TESTS::

            sage: from sage.repl.ipython_kernel.kernel import SageKernel
            sage: k = SageKernel.__new__(SageKernel)
            sage: k.pre_handler_hook()
            sage: k.saved_sigint_handler
            <cyfunction python_check_interrupt at ...>
        """
        from cysignals import init_cysignals
        self.saved_sigint_handler = init_cysignals()
