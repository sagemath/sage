"""
Installing the Sage IPython Kernel

Kernels have to register themselves with IPython so that they appear
in the IPython notebook's kernel drop-down. This is done by
:class:`SageKernelSpec`.
"""

import os
import errno
import copy

from IPython.kernel.kernelspec import (
    get_kernel_spec, install_kernel_spec, NoSuchKernel)
from IPython.utils.path import get_ipython_dir

from sage.env import (
    SAGE_ROOT, SAGE_DOC, SAGE_LOCAL, SAGE_EXTCODE,
    SAGE_VERSION
)
from sage.misc.temporary_file import tmp_dir


class SageKernelSpec(object):

    def __init__(self):
        """
        Utility to manage Sage kernels

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._display_name    # random output
            'Sage 6.6.beta2'
        """
        self._display_name = 'Sage {0}'.format(SAGE_VERSION)
        self._ipython_dir = get_ipython_dir()
        self._mkdirs()

    def _mkdirs(self):
        """
        Create necessary parent directories

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._mkdirs()
            sage: nbextensions = os.path.join(spec._ipython_dir, 'nbextensions')
            sage: os.path.exists(nbextensions)
            True
        """
        def mkdir_p(*path_components):
            path = os.path.join(*path_components)
            try:
                os.makedirs(path)
            except OSError as err:
                if err.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else:
                    raise
        mkdir_p(self._ipython_dir, 'nbextensions')

    @classmethod
    def identifier(self):
        """
        Internal identifier for the Sage kernel

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: SageKernelSpec.identifier()    # random output
            'sage_6_6_beta3'
            sage: SageKernelSpec.identifier().startswith('sage_')
            True
        """
        return 'Sage {0}'.format(SAGE_VERSION).lower().replace(' ', '_').replace('.', '_')
        
    def symlink(self, src, dst):
        """
        Symlink ``src`` to ``dst``

        This is not an atomic operation.

        Already-existing symlinks will be deleted, already existing
        non-empty directories will be kept.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: path = tmp_dir()
            sage: spec.symlink(os.path.join(path, 'a'), os.path.join(path, 'b'))
            sage: os.listdir(path)
            ['b']
        """
        try:
            os.remove(dst)
        except OSError as err:
            if err.errno == errno.EEXIST:
                return
        os.symlink(src, dst)
        
    def use_local_mathjax(self):
        """
        Symlink Sage's Mathjax Install to the IPython notebook.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: from IPython.utils.path import get_ipython_dir
            sage: spec = SageKernelSpec()
            sage: spec.use_local_mathjax()
            sage: ipython_dir = get_ipython_dir()
            sage: mathjax = os.path.join(ipython_dir, 'nbextensions', 'mathjax')
            sage: os.path.exists(mathjax)
            True
        """
        src = os.path.join(SAGE_LOCAL, 'share', 'mathjax')
        dst = os.path.join(self._ipython_dir, 'nbextensions', 'mathjax')
        self.symlink(src, dst)

    def use_local_jsmol(self):
        src = os.path.join(SAGE_LOCAL, 'share', 'jsmol')
        dst = os.path.join(self._ipython_dir, 'nbextensions', 'jsmol')
        self.symlink(src, dst)

    def _kernel_cmd(self):
        """
        Helper to construct the Sage kernel command.
        
        OUTPUT:

        List of strings. The command to start a new Sage kernel.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._kernel_cmd()
            ['/.../sage',
             '-python',
             '-m',
             'sage.repl.ipython_kernel',
             '-f',
             '{connection_file}']
        """
        return [
            os.path.join(SAGE_ROOT, 'sage'),
            '-python',
            '-m', 'sage.repl.ipython_kernel',
            '-f', '{connection_file}',
        ]
        
    def kernel_spec(self):
        """
        Return the kernel spec as Python dictionary

        OUTPUT:

        A dictionary. See the IPython documentation for details.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec.kernel_spec()
            {'argv': ..., 'display_name': 'Sage ...'}
        """
        return dict(
            argv=self._kernel_cmd(),
            display_name=self._display_name,
        )
    
    def _install_spec(self):
        """
        Install the Sage IPython kernel
        
        It is safe to call this method multiple times, only one Sage
        kernel spec is ever installed for any given Sage
        version. However, it resets the IPython kernel spec directory
        so additional resources symlinked there are lost. See
        :meth:`symlink_resources`.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._install_spec()    # not tested
        """
        import json
        temp = tmp_dir()
        kernel_spec = os.path.join(temp, 'kernel.json')
        with open(kernel_spec, 'w') as f:
            json.dump(self.kernel_spec(), f)
        identifier = self.identifier()
        install_kernel_spec(temp, identifier, user=True, replace=True)
        self._spec = get_kernel_spec(identifier)

    def _symlink_resources(self):
        """
        Symlink miscellaneous resources

        This method symlinks additional resources (like the Sage
        documentation) into the Sage kernel directory. This is
        necessary to make the help links in the notebook work.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._install_spec()         # not tested
            sage: spec._symlink_resources()    # not tested
        """
        assert self._spec, 'call _install_spec() first'
        spec_dir = self._spec.resource_dir
        path = os.path.join(SAGE_EXTCODE, 'notebook-ipython')
        for filename in os.listdir(path):
            self.symlink(
                os.path.join(path, filename),
                os.path.join(spec_dir, filename)
            )
        self.symlink(
            os.path.join(SAGE_DOC, 'output', 'html', 'en'),
            os.path.join(spec_dir, 'doc')
        )
        
    @classmethod
    def update(cls):
        """
        Configure the IPython notebook for the Sage kernel
        
        This method does everything necessary to use the Sage kernel,
        you should never need to call any of the other methods
        directly.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec.update()
        """
        instance = cls()
        instance.use_local_mathjax()
        instance.use_local_jsmol()
        instance._install_spec()
        instance._symlink_resources()

        
def have_prerequisites():
    """
    Check that we have all prerequisites to run the IPython notebook.

    In particular, the IPython notebook requires OpenSSL whether or
    not you are using https. See trac:`17318`.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.install import have_prerequisites
        sage: have_prerequisites() in [True, False]
        True
    """
    try:
        from IPython.html.notebookapp import NotebookApp
        return True
    except ImportError as err:
        return False




