# sage_setup: distribution = sagemath-repl
"""
Installing the SageMath Jupyter Kernel and Extensions

Kernels have to register themselves with Jupyter so that they appear
in the Jupyter notebook's kernel drop-down. This is done by
:class:`SageKernelSpec`.

.. NOTE::

    The doctests in this module run in a temporary directory as the involved
    directories might be different during runs of the tests and actual
    installation and because we might be lacking write permission to places
    such as ``/usr/share``.
"""

import errno
import os
import warnings

from sage.env import (
    SAGE_DOC,
    SAGE_EXTCODE,
    SAGE_VENV,
    SAGE_VERSION,
)


class SageKernelSpec():

    def __init__(self, prefix=None):
        """
        Utility to manage SageMath kernels and extensions.

        INPUT:

        - ``prefix`` -- (default: ``sys.prefix``)
          directory for the installation prefix

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: prefix = tmp_dir()
            sage: spec = SageKernelSpec(prefix=prefix)
            sage: spec._display_name    # random output
            'SageMath 6.9'
            sage: spec.kernel_dir == SageKernelSpec(prefix=prefix).kernel_dir
            True
        """
        self._display_name = 'SageMath {0}'.format(SAGE_VERSION)
        if prefix is None:
            from sys import prefix
        jupyter_dir = os.path.join(prefix, "share", "jupyter")
        self.nbextensions_dir = os.path.join(jupyter_dir, "nbextensions")
        self.kernel_dir = os.path.join(jupyter_dir, "kernels", self.identifier())
        self._mkdirs()

    def _mkdirs(self):
        """
        Create necessary parent directories.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec._mkdirs()
            sage: os.path.isdir(spec.nbextensions_dir)
            True
        """
        def mkdir_p(path):
            try:
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise
        mkdir_p(self.nbextensions_dir)
        mkdir_p(self.kernel_dir)

    @classmethod
    def identifier(cls):
        """
        Internal identifier for the SageMath kernel.

        OUTPUT: the string ``'sagemath'``

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: SageKernelSpec.identifier()
            'sagemath'
        """
        return 'sagemath'

    def symlink(self, src, dst):
        """
        Symlink ``src`` to ``dst``.

        This is not an atomic operation.

        Already-existing symlinks will be deleted, already existing
        non-empty directories will be kept.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
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

    def use_local_threejs(self):
        """
        Symlink threejs to the Jupyter notebook.

        EXAMPLES::

            sage: # needs threejs
            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec.use_local_threejs()
            sage: threejs = os.path.join(spec.nbextensions_dir, 'threejs-sage')
            sage: os.path.isdir(threejs)
            True
        """
        from sage.features.threejs import Threejs
        if not Threejs().is_present():
            return
        src = os.path.dirname(os.path.dirname(Threejs().absolute_filename()))
        dst = os.path.join(self.nbextensions_dir, 'threejs-sage')
        self.symlink(src, dst)

    def _kernel_cmd(self):
        """
        Helper to construct the SageMath kernel command.

        OUTPUT: list of strings; the command to start a new SageMath kernel

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec._kernel_cmd()
            ['/.../sage',
             '--python',
             '-m',
             'sage.repl.ipython_kernel',
             '-f',
             '{connection_file}']
        """
        return [
            os.path.join(SAGE_VENV, 'bin', 'sage'),
            '--python',
            '-m', 'sage.repl.ipython_kernel',
            '-f', '{connection_file}',
        ]

    def kernel_spec(self):
        """
        Return the kernel spec as Python dictionary.

        OUTPUT: a dictionary; see the Jupyter documentation for details

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec.kernel_spec()
            {'argv': ..., 'display_name': 'SageMath ...', 'language': 'sage'}
        """
        return dict(
            argv=self._kernel_cmd(),
            display_name=self._display_name,
            language='sage',
        )

    def _install_spec(self):
        """
        Install the SageMath Jupyter kernel.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec._install_spec()
        """
        jsonfile = os.path.join(self.kernel_dir, "kernel.json")
        import json
        with open(jsonfile, 'w') as f:
            json.dump(self.kernel_spec(), f)

    def _symlink_resources(self):
        """
        Symlink miscellaneous resources.

        This method symlinks additional resources (like the SageMath
        documentation) into the SageMath kernel directory. This is
        necessary to make the help links in the notebook work.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec._install_spec()
            sage: spec._symlink_resources()
        """
        path = os.path.join(SAGE_EXTCODE, 'notebook-ipython')
        for filename in os.listdir(path):
            self.symlink(
                os.path.join(path, filename),
                os.path.join(self.kernel_dir, filename)
            )
        self.symlink(
            SAGE_DOC,
            os.path.join(self.kernel_dir, 'doc')
        )

    @classmethod
    def update(cls, *args, **kwds):
        """
        Configure the Jupyter notebook for the SageMath kernel.

        This method does everything necessary to use the SageMath kernel,
        you should never need to call any of the other methods
        directly.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: SageKernelSpec.update(prefix=tmp_dir())
        """
        instance = cls(*args, **kwds)
        instance.use_local_threejs()
        instance._install_spec()
        instance._symlink_resources()

    @classmethod
    def check(cls):
        """
        Check that the SageMath kernel can be discovered by its name (sagemath).

        This method issues a warning if it cannot -- either because it is not installed,
        or it is shadowed by a different kernel of this name, or Jupyter is
        misconfigured in a different way.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: SageKernelSpec.check()  # random
        """
        from jupyter_client.kernelspec import NoSuchKernel, get_kernel_spec
        ident = cls.identifier()
        try:
            spec = get_kernel_spec(ident)
        except NoSuchKernel:
            warnings.warn(f'no kernel named {ident} is accessible; '
                          'check your Jupyter configuration '
                          '(see https://docs.jupyter.org/en/latest/use/jupyter-directories.html)')
        else:
            from pathlib import Path
            if Path(spec.argv[0]).resolve() != Path(os.path.join(SAGE_VENV, 'bin', 'sage')).resolve():
                warnings.warn(f'the kernel named {ident} does not seem to correspond to this '
                              'installation of SageMath; check your Jupyter configuration '
                              '(see https://docs.jupyter.org/en/latest/use/jupyter-directories.html)')


def have_prerequisites(debug=True):
    """
    Check that we have all prerequisites to run the Jupyter notebook.

    In particular, the Jupyter notebook requires OpenSSL whether or
    not you are using https. See :issue:`17318`.

    INPUT:

    - ``debug`` -- boolean (default: ``True``); whether to print debug
      information in case that prerequisites are missing

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.install import have_prerequisites
        sage: have_prerequisites(debug=False) in [True, False]
        True
    """
    try:
        from notebook.notebookapp import NotebookApp
        return True
    except ImportError:
        if debug:
            import traceback
            traceback.print_exc()
        return False
