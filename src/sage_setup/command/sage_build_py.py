from setuptools.command.build_py import build_py


class sage_build_py(build_py):

    def find_package_modules(self, package, package_dir):
        r"""
        Do not expand packages to modules.

        Our :func:`sage_setup.find.find_python_sources` already lists all modules.
        """
        return []
