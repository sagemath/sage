# Import setuptools before importing distutils, so that setuptools
# can replace distutils by its own vendored copy.
import setuptools

from distutils import log
from distutils.command.build import build

class sage_build(build):
    sub_commands = [('build_cython', lambda *args: True)] + build.sub_commands

    def run_autogen(self):
        """
        Generate auto-generated sources.

        This must be done before building the python modules,
        see :trac:`22106`.
        """
        from sage_setup.autogen import autogen_all
        from sage_setup.find import find_python_sources
        from sage.env import SAGE_SRC

        log.info("Generating auto-generated sources")

        pkgs = autogen_all()
        python_packages, python_modules, cython_modules = find_python_sources(
            SAGE_SRC, [ pkg.replace('.', '/') for pkg in pkgs])

        for pkg in python_packages:
            if pkg not in self.distribution.packages:
                self.distribution.packages.append(pkg)

        for cython_module in cython_modules:
            self.distribution.ext_modules.append(cython_module)

    def run(self):
        self.run_autogen()
        build.run(self)
