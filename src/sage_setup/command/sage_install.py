#########################################################
### Install Jupyter kernel spec
#########################################################

import logging

# Import setuptools before importing distutils, so that setuptools
# can replace distutils by its own vendored copy.
import setuptools
from setuptools.command.install import install
from setuptools.command.develop import develop

logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'  # keep as distutils.log the same
)
log = logging.getLogger(__name__)


class install_kernel_spec_mixin:

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


class sage_install(install, install_kernel_spec_mixin):

    def run(self):
        install.run(self)
        self.install_kernel_spec()


class sage_develop(develop, install_kernel_spec_mixin):

    def run(self):
        develop.run(self)
        if not self.uninstall:
            self.install_kernel_spec()
