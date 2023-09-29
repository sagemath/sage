import os
import sys
import shutil
import sysconfig
import platform
import fnmatch

from pathlib import Path

from setuptools import setup
from distutils.command.build_scripts import build_scripts as distutils_build_scripts
from setuptools.command.build_py import build_py as setuptools_build_py
from setuptools.command.egg_info import egg_info as setuptools_egg_info
from distutils.errors import (DistutilsSetupError, DistutilsModuleError,
                              DistutilsOptionError)

class build_py(setuptools_build_py):

    def run(self):
        DOT_SAGE = os.environ.get('DOT_SAGE', os.path.join(os.environ.get('HOME'), '.sage'))
        HERE = os.path.dirname(__file__)
        with open(os.path.join(HERE, 'VERSION.txt')) as f:
            sage_version = f.read().strip()
        # After #30534, SAGE_LOCAL no longer contains any Python.  So we key the SAGE_ROOT only to Sage version
        # and architecture.
        system = platform.system()
        machine = platform.machine()
        arch_tag = f'{system}-{machine}'
        # TODO: These two should be user-configurable with options passed to "setup.py install"
        SAGE_ROOT = os.path.join(DOT_SAGE, f'sage-{sage_version}-{arch_tag}-conda')

        def ignore(path, names):
            # exclude all embedded src trees
            if fnmatch.fnmatch(path, f'*/build/pkgs/*'):
                return ['src']
            ### ignore more stuff --- .tox etc.
            return [name for name in names
                    if name in ('.tox',)]

        if os.path.exists(os.path.join(SAGE_ROOT, 'config.status')):
            print(f'Reusing SAGE_ROOT={SAGE_ROOT}')
        else:
            # config.status and other configure output has to be writable.
            # So (until the Sage distribution supports VPATH builds - #21469), we have to make a copy of sage_root.
            try:
                shutil.copytree('sage_root', SAGE_ROOT,
                                ignore=ignore)  # will fail if already exists
            except Exception as e:
                raise DistutilsSetupError(f"the directory SAGE_ROOT={SAGE_ROOT} already exists but it is not configured.  Please remove it and try again. (Exception: {e})")
            cmd = f"cd {SAGE_ROOT} && ./configure --enable-build-as-root --with-system-python3=force --disable-notebook --disable-sagelib --disable-sage_conf --disable-doc"
            cmd += ' --with-python=$CONDA_PREFIX/bin/python --prefix="$CONDA_PREFIX" --enable-system-site-packages'
            cmd += ' $(for pkg in $(build/bin/sage-package list :standard: --exclude rpy2 --has-file spkg-configure.m4 --has-file distros/conda.txt); do echo --with-system-$pkg=force; done)'
            print(f"Running {cmd}")
            sys.stdout.flush()
            if os.system(cmd) != 0:
                print(f"configure failed; this may be caused by missing build prerequisites.")
                sys.stdout.flush()
                PREREQ_SPKG = "_prereq bzip2 xz libffi"  # includes python3 SPKG_DEPCHECK packages
                os.system(f'cd {SAGE_ROOT} && export SYSTEM=$(build/bin/sage-guess-package-system 2>/dev/null) && export PACKAGES="$(build/bin/sage-get-system-packages $SYSTEM {PREREQ_SPKG})" && [ -n "$PACKAGES" ] && echo "You can install the required build prerequisites using the following shell command" && echo "" && build/bin/sage-print-system-package-command $SYSTEM --verbose --sudo install $PACKAGES && echo ""')
                raise DistutilsSetupError("configure failed")

        # In this mode, we never run "make".

        # Install configuration
        shutil.copyfile(os.path.join(SAGE_ROOT, 'pkgs', 'sage-conf', '_sage_conf', '_conf.py'),
                        os.path.join(HERE, '_sage_conf', '_conf.py'))
        shutil.copyfile(os.path.join(SAGE_ROOT, 'src', 'bin', 'sage-env-config'),
                        os.path.join(HERE, 'bin', 'sage-env-config'))
        setuptools_build_py.run(self)

class build_scripts(distutils_build_scripts):

    def run(self):
        self.distribution.scripts.append(os.path.join('bin', 'sage-env-config'))
        if not self.distribution.entry_points:
            self.entry_points = self.distribution.entry_points = dict()
        distutils_build_scripts.run(self)

setup(
    cmdclass=dict(build_py=build_py, build_scripts=build_scripts)
)
