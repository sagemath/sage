import os
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as setuptools_build_py
from setuptools.command.editable_wheel import (
    editable_wheel as setuptools_editable_wheel,
)
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

        builddir = here / "sage_conf" / "builddir"
        cmd = f"cd {root} && meson setup {builddir} --wipe"
        print(f"Running {cmd}")
        sys.stdout.flush()
        if os.system(cmd) != 0:
            raise SetupError("configure failed")

        # Write build info
        with open(builddir / 'build_info.py', 'w', encoding="utf-8") as build_info:
            build_info.write(f'SAGE_ROOT = "{root}"\n')
            build_info.write(f'CONDA_PREFIX = "{conda_prefix}"\n')


class editable_wheel(setuptools_editable_wheel):
    r"""
    Customized so that exceptions raised by our build_py
    do not lead to the "Customization incompatible with editable install" message
    """
    _safely_run = setuptools_editable_wheel.run_command


setup(
    cmdclass=dict(
        build_py=build_py, editable_wheel=editable_wheel
    )
)
