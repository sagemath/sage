import os
import sys
from distutils.command.build_scripts import \
    build_scripts as distutils_build_scripts
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as setuptools_build_py
from setuptools.errors import SetupError


class build_py(setuptools_build_py):
    def run(self):
        here = Path(__file__).parent
        if self.editable_mode:
            root = here.parent.parent
        else:
            raise SetupError('Not supported')

        conda_prefix = os.environ.get('CONDA_PREFIX', '')
        if not conda_prefix:
            raise SetupError(
                'No conda environment is active. '
                'See https://doc.sagemath.org/html/en/installation/conda.html on how to get started.'
            )

        builddir = here / "builddir"
        cmd = f"cd {root} && meson setup {builddir} --wipe"
        print(f"Running {cmd}")
        sys.stdout.flush()
        if os.system(cmd) != 0:
            raise SetupError("configure failed")

        # Write build info
        with open(builddir / 'build_info.py', 'w', encoding="utf-8") as build_info:
            build_info.write(f'SAGE_ROOT = "{root}"\n')
            build_info.write(f'CONDA_PREFIX = "{conda_prefix}"\n')


class build_scripts(distutils_build_scripts):
    def run(self):
        self.distribution.scripts.append(os.path.join('bin', 'sage-env-config'))
        if not self.distribution.entry_points:
            self.entry_points = self.distribution.entry_points = dict()
        distutils_build_scripts.run(self)


setup(
    cmdclass=dict(
        build_py=build_py, build_scripts=build_scripts
    )
)
