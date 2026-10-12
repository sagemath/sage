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

import os
import warnings

from sage.env import (
    SAGE_DOC,
    SAGE_EXTCODE,
    SAGE_LOCAL,
    SAGE_VERSION,
)


class SageKernelSpec:

    def __init__(self, prefix=None, portable=False, jupyter_dir=None):
        """
        Utility to manage SageMath kernels and extensions.

        INPUT:

        - ``prefix`` -- (default: ``sys.prefix``)
          directory for the installation prefix; the Jupyter files are
          installed under ``prefix/share/jupyter``. Ignored if
          ``jupyter_dir`` is given.

        - ``portable`` -- boolean (default: ``False``); if ``True``, generate
          a kernel spec that launches the kernel via the absolute path to this
          Python interpreter and adds Sage's ``bin`` directories to ``PATH``
          (and its ``lib`` directory to the dynamic-loader path), so that it
          works from a Jupyter server running outside Sage's virtual
          environment

        - ``jupyter_dir`` -- (default: ``prefix/share/jupyter``) the Jupyter
          data directory to install into. This is used to install into a
          location such as the per-user Jupyter directory that does not follow
          the ``prefix/share/jupyter`` layout.

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
        self._portable = portable
        if jupyter_dir is None:
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

    def symlink(self, src, dst, replace_dir=False):
        """
        Symlink ``src`` to ``dst``.

        This is not an atomic operation.

        Already-existing files and symlinks at ``dst`` are removed first so
        that they can be replaced by a fresh symlink. An already-existing
        *directory* at ``dst`` is kept (and ``src`` is not linked), unless
        ``replace_dir`` is ``True``, in which case it is removed first.

        INPUT:

        - ``replace_dir`` -- boolean (default: ``False``); if ``True``, replace
          an existing directory at ``dst``. This is used to clean up a real
          directory left behind by older installations that copied resources
          (such as the documentation) here instead of symlinking them; see
          :issue:`40605`.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: path = tmp_dir()
            sage: spec.symlink(os.path.join(path, 'a'), os.path.join(path, 'b'))
            sage: os.listdir(path)
            ['b']

        An existing directory at the destination is kept by default (to avoid
        deleting unrelated data), but is replaced with ``replace_dir=True``::

            sage: d = os.path.join(path, 'c')
            sage: os.mkdir(d)
            sage: with open(os.path.join(d, 'stale'), 'w') as f:
            ....:     _ = f.write('x')
            sage: spec.symlink(os.path.join(path, 'a'), d)
            sage: os.path.isdir(d) and not os.path.islink(d)
            True
            sage: spec.symlink(os.path.join(path, 'a'), d, replace_dir=True)
            sage: os.path.islink(d)
            True
        """
        if os.path.islink(dst) or os.path.isfile(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            if not replace_dir:
                return
            import shutil
            shutil.rmtree(dst)
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
            ['python3',
             '-m',
             'sage.repl.ipython_kernel',
             '-f',
             '{connection_file}']

        A *portable* kernel spec uses the absolute path to this Python
        interpreter, so that it can be started from a Jupyter server running
        outside Sage's virtual environment (see also :meth:`_kernel_env`)::

            sage: spec = SageKernelSpec(prefix=tmp_dir(), portable=True)
            sage: cmd = spec._kernel_cmd()
            sage: import sys
            sage: cmd[0] == os.path.abspath(sys.executable)
            True
            sage: cmd[1:]
            ['-m', 'sage.repl.ipython_kernel', '-f', '{connection_file}']
        """
        if self._portable:
            import sys
            executable = os.path.abspath(sys.executable)
        else:
            executable = 'python3'
        return [
            executable,
            '-m', 'sage.repl.ipython_kernel',
            '-f', '{connection_file}',
        ]

    def _kernel_env(self):
        """
        Environment for a portable kernel spec.

        This puts the ``bin`` directory containing this Python interpreter on
        ``PATH``, along with Sage's installation ``bin`` directory when that is
        separate. This lets helper programs such as ``gap`` be found even when
        the kernel is launched by a Jupyter server running outside Sage's
        virtual environment. When Sage's installation prefix is known, its
        ``lib`` directory is set as ``LD_LIBRARY_PATH`` (and
        ``DYLD_LIBRARY_PATH`` on macOS) so that the dynamic loader can find
        Sage's shared libraries on installations that do not encode an RPATH
        (for example some source builds). The loader variables are not
        appended to, since the launching environment usually does not define
        them; the default loader search path (and any RPATH) still applies.
        Jupyter expands the ``${PATH}`` reference at launch time.

        OUTPUT: a dictionary mapping environment variable names to values

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir(), portable=True)
            sage: import sys
            sage: bindir = os.path.dirname(os.path.abspath(sys.executable))
            sage: spec._kernel_env()['PATH'].endswith(os.pathsep + '${PATH}')
            True

        If Sage's installation prefix is separate from its Python virtual
        environment, the portable kernel also needs the installation ``bin``
        directory on ``PATH`` and the ``lib`` directory on the library path::

            sage: from sage.repl.ipython_kernel import install
            sage: old_sage_local = install.SAGE_LOCAL
            sage: try:
            ....:     install.SAGE_LOCAL = os.path.join(tmp_dir(), 'local')
            ....:     spec = install.SageKernelSpec(prefix=tmp_dir(), portable=True)
            ....:     env = spec._kernel_env()
            ....:     sage_local_bin = os.path.abspath(os.path.join(install.SAGE_LOCAL, 'bin'))
            ....:     sage_local_lib = os.path.abspath(os.path.join(install.SAGE_LOCAL, 'lib'))
            ....:     ok = env['PATH'].split(os.pathsep) == [bindir, sage_local_bin, '${PATH}']
            ....:     ok &= env['LD_LIBRARY_PATH'] == sage_local_lib
            ....:     ok &= env['DYLD_LIBRARY_PATH'] == sage_local_lib
            ....: finally:
            ....:     install.SAGE_LOCAL = old_sage_local
            sage: ok
            True
        """
        import sys
        paths = [os.path.dirname(os.path.abspath(sys.executable))]
        env = {}
        if SAGE_LOCAL:
            sage_local_bin = os.path.abspath(os.path.join(SAGE_LOCAL, 'bin'))
            if sage_local_bin not in paths:
                paths.append(sage_local_bin)
            sage_local_lib = os.path.abspath(os.path.join(SAGE_LOCAL, 'lib'))
            for var in ('LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH'):
                env[var] = sage_local_lib
        paths.append('${PATH}')
        env['PATH'] = os.pathsep.join(paths)
        return env

    def kernel_spec(self):
        """
        Return the kernel spec as Python dictionary.

        OUTPUT: a dictionary; see the Jupyter documentation for details

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec(prefix=tmp_dir())
            sage: spec.kernel_spec()
            {'argv': ..., 'display_name': 'SageMath ...', 'language': 'sage', 'metadata': {'debugger': True}}

        A portable spec additionally carries an ``env`` entry::

            sage: spec = SageKernelSpec(prefix=tmp_dir(), portable=True)
            sage: 'PATH' in spec.kernel_spec()['env']
            True
        """
        spec = dict(
            argv=self._kernel_cmd(),
            display_name=self._display_name,
            language='sage',
            metadata=dict(debugger=True),
        )
        if self._portable:
            spec['env'] = self._kernel_env()
        return spec

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
            sage: 'logo.svg' in os.listdir(spec.kernel_dir)
            True

        The ``kernel.json.in`` template is not part of the kernel spec and is
        not symlinked into the kernel directory::

            sage: 'kernel.json.in' in os.listdir(spec.kernel_dir)
            False
        """
        path = os.path.join(SAGE_EXTCODE, 'notebook-ipython')
        for filename in os.listdir(path):
            if filename.endswith('.in'):
                # Skip templates such as kernel.json.in; the configured
                # kernel.json is written by _install_spec().
                continue
            self.symlink(
                os.path.join(path, filename),
                os.path.join(self.kernel_dir, filename)
            )
        # replace_dir=True cleans up a real doc/ directory left behind by
        # older installations that copied the documentation here (see
        # :issue:`40605`).
        self.symlink(
            SAGE_DOC,
            os.path.join(self.kernel_dir, 'doc'),
            replace_dir=True
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
            warnings.warn(f'No kernel named {ident} is accessible; '
                          'check your Jupyter configuration '
                          '(see https://docs.jupyter.org/en/latest/use/jupyter-directories.html).')
        else:
            import sys
            from pathlib import Path
            from sage.features import Executable
            kernel_executable_feature = Executable(name=spec.argv[0], executable=spec.argv[0])
            if not kernel_executable_feature.is_present():
                warnings.warn(f'The kernel named {ident} does not seem to be runnable; '
                              'check your Jupyter configuration '
                              '(see https://docs.jupyter.org/en/latest/use/jupyter-directories.html).')
                return
            kernel_executable = kernel_executable_feature.absolute_filename()
            if Path(kernel_executable).resolve() != Path(sys.executable).resolve():
                warnings.warn(f'The kernel named {ident} does not seem to correspond to this '
                              'installation of SageMath; check your Jupyter configuration '
                              '(see https://docs.jupyter.org/en/latest/use/jupyter-directories.html).')


def have_prerequisites(debug=True) -> bool:
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


def main(argv=None):
    r"""
    Command-line entry point for ``sage --jupyter-kernel``.

    INPUT:

    - ``argv`` -- (default: ``sys.argv[1:]``) list of command-line arguments

    EXAMPLES:

    Installing a kernel spec into a given prefix writes a ``kernel.json`` that,
    by default (``--portable``), uses the absolute path to this Python
    interpreter and carries an ``env`` entry putting the Sage ``bin`` directory
    on ``PATH``::

        sage: from sage.repl.ipython_kernel.install import main
        sage: import json, os, sys
        sage: prefix = tmp_dir()
        sage: main(['install', '--prefix', prefix])
        sage: kernel_json = os.path.join(prefix, 'share', 'jupyter',
        ....:                            'kernels', 'sagemath', 'kernel.json')
        sage: with open(kernel_json) as f:
        ....:     spec = json.load(f)
        sage: spec['argv'][0] == os.path.abspath(sys.executable)
        True
        sage: spec['argv'][1:]
        ['-m', 'sage.repl.ipython_kernel', '-f', '{connection_file}']
        sage: 'PATH' in spec['env']
        True

    With ``--no-portable`` the in-venv spec (using ``python3``) is written
    instead::

        sage: prefix = tmp_dir()
        sage: main(['install', '--prefix', prefix, '--no-portable'])
        sage: kernel_json = os.path.join(prefix, 'share', 'jupyter',
        ....:                            'kernels', 'sagemath', 'kernel.json')
        sage: with open(kernel_json) as f:
        ....:     json.load(f)['argv']
        ['python3', '-m', 'sage.repl.ipython_kernel', '-f', '{connection_file}']

    The ``--user`` option installs into the per-user Jupyter directory
    (here redirected to a temporary location)::

        sage: import jupyter_core.paths
        sage: datadir = tmp_dir()
        sage: old = jupyter_core.paths.jupyter_data_dir
        sage: try:
        ....:     jupyter_core.paths.jupyter_data_dir = lambda: datadir
        ....:     main(['install', '--user'])
        ....: finally:
        ....:     jupyter_core.paths.jupyter_data_dir = old
        sage: os.path.exists(os.path.join(datadir, 'kernels', 'sagemath', 'kernel.json'))
        True
    """
    import argparse
    parser = argparse.ArgumentParser(
        prog='sage --jupyter-kernel',
        description='Install the SageMath Jupyter kernel spec.')
    subparsers = parser.add_subparsers(dest='command', required=True)
    install = subparsers.add_parser(
        'install', help='install the SageMath Jupyter kernel spec')
    location = install.add_mutually_exclusive_group()
    location.add_argument(
        '--user', action='store_true',
        help='install into the per-user Jupyter directory')
    location.add_argument(
        '--sys-prefix', action='store_true',
        help='install into sys.prefix of the active environment (default)')
    location.add_argument(
        '--prefix', default=None,
        help='install into PREFIX/share/jupyter')
    install.add_argument(
        '--portable', action=argparse.BooleanOptionalAction, default=True,
        help="generate a spec with an absolute interpreter path and PATH / "
             "library-path environment entries so that it works from a "
             "Jupyter server outside Sage's virtual environment (the default). "
             "Use --no-portable to write a spec that uses a bare 'python3', "
             "which only works when the Jupyter server itself runs inside "
             "Sage's environment")
    args = parser.parse_args(argv)

    if args.command == 'install':
        kwds = {'portable': args.portable}
        if args.user:
            from jupyter_core.paths import jupyter_data_dir
            kwds['jupyter_dir'] = jupyter_data_dir()
        elif args.prefix:
            kwds['prefix'] = args.prefix
        else:
            # --sys-prefix (the default)
            import sys
            kwds['prefix'] = sys.prefix
        SageKernelSpec.update(**kwds)


if __name__ == '__main__':
    main()
