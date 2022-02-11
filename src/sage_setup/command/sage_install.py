#########################################################
### Install Jupyter kernel spec and clean stale files
#########################################################

import os
import time

# Import setuptools before importing distutils, so that setuptools
# can replace distutils by its own vendored copy.
import setuptools

from distutils import log
from distutils.command.install import install

class sage_install(install):

    def run(self):
        install.run(self)
        self.install_kernel_spec()

    def install_kernel_spec(self):
        """
        Install the Jupyter kernel spec.

        .. NOTE::

            The files are generated, not copied. Therefore, we cannot
            use ``data_files`` for this.
        """
        from sage.repl.ipython_kernel.install import SageKernelSpec
        # Jupyter packages typically use the data_files option to
        # setup() to install kernels and nbextensions. So we should use
        # the install_data directory for installing our Jupyter files.
        SageKernelSpec.update(prefix=self.install_data)

class sage_install_and_clean(sage_install):

    def run(self):
        sage_install.run(self)
        t = time.time()
        self.clean_stale_files()
        log.info('Finished cleaning, time: %.2f seconds.' % (time.time() - t))

    def clean_stale_files(self):
        """
        Remove stale installed files.

        This removes files which are built/installed but which do not
        exist in the Sage sources (typically because some source file
        has been deleted). Files are removed from the build directory
        ``build/lib-*`` and from the install directory ``site-packages``.
        """
        dist = self.distribution
        cmd_build_py = self.get_finalized_command("build_py")
        cmd_build_cython = self.get_finalized_command("build_cython")

        # Determine all Python modules inside all packages
        py_modules = []
        for package in dist.packages:
            package_dir = cmd_build_py.get_package_dir(package)
            py_modules += cmd_build_py.find_package_modules(package, package_dir)
        # modules is a list of triples (package, module, module_file).
        # Construct the complete module name from this.
        py_modules = ["{0}.{1}".format(*m) for m in py_modules]

        # Determine all files of package data and Cythonized package files
        # example of entries of cmd_build_cython.get_cythonized_package_files():
        #   ('sage/media', ['./sage/media/channels.pyx'])
        data_files = cmd_build_cython.get_cythonized_package_files()
        # examples of entries of build_py.data_files:
        #   ('sage.libs.gap', 'sage/libs/gap', 'build/lib.macosx-10.9-x86_64-3.7/sage/libs/gap', ['sage.gaprc'])
        #   ('sage', 'sage', 'build/lib.macosx-10.9-x86_64-3.7/sage', ['ext_data/nodoctest.py', 'ext_data/kenzo/S4.txt', ...])
        nobase_data_files = [(src_dir, [os.path.join(src_dir, filename) for filename in filenames])
                             for package, src_dir, build_dir, filenames in cmd_build_py.data_files]

        # Clean install directory (usually, purelib and platlib are the same)
        # and build directory.
        output_dirs = [self.install_purelib, self.install_platlib, self.build_lib]
        from sage_setup.clean import clean_install_dir
        for output_dir in set(output_dirs):
            log.info('- cleaning {0}'.format(output_dir))
            clean_install_dir(output_dir,
                    dist.packages,
                    py_modules,
                    dist.ext_modules,
                    data_files,
                    nobase_data_files)
