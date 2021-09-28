import os
import sys
import shutil
import sysconfig
import platform

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
        # For convenience, set up the homebrew env automatically. This is a no-op if homebrew is not present.
        SETENV = '(. ./.homebrew-build-env 2> /dev/null || :)'
        # After #30534, SAGE_LOCAL no longer contains any Python.  So we key the SAGE_ROOT only to Sage version
        # and architecture.
        system = platform.system()
        machine = platform.machine()
        arch_tag = f'{system}-{machine}'
        # TODO: These two should be user-configurable with options passed to "setup.py install"
        SAGE_ROOT = os.path.join(DOT_SAGE, f'sage-{sage_version}-{arch_tag}')
        SAGE_LOCAL = os.path.join(SAGE_ROOT, 'local')
        if os.path.exists(os.path.join(SAGE_ROOT, 'config.status')):
            print(f'Reusing SAGE_ROOT={SAGE_ROOT}')
        else:
            # config.status and other configure output has to be writable.
            # So (until the Sage distribution supports VPATH builds - #21469), we have to make a copy of sage_root.
            try:
                shutil.copytree('sage_root', SAGE_ROOT)  # will fail if already exists
            except Exception:
                raise DistutilsSetupError(f"the directory SAGE_ROOT={SAGE_ROOT} already exists but it is not configured.  Please remove it and try again.")
            cmd = f"cd {SAGE_ROOT} && {SETENV} && ./configure --prefix={SAGE_LOCAL} --with-python={sys.executable} --enable-build-as-root --enable-download-from-upstream-url --with-system-python3=force --with-sage-venv --disable-notebook --disable-sagelib"
            print(f"Running {cmd}")
            sys.stdout.flush()
            if os.system(cmd) != 0:
                print(f"configure failed; this may be caused by missing build prerequisites.")
                sys.stdout.flush()
                PREREQ_SPKG = "_prereq bzip2 xz libffi"  # includes python3 SPKG_DEPCHECK packages
                os.system(f'cd {SAGE_ROOT} && export SYSTEM=$(build/bin/sage-guess-package-system 2>/dev/null) && export PACKAGES="$(build/bin/sage-get-system-packages $SYSTEM {PREREQ_SPKG})" && [ -n "$PACKAGES" ] && echo "You can install the required build prerequisites using the following shell command" && echo "" && build/bin/sage-print-system-package-command $SYSTEM --verbose --sudo install $PACKAGES && echo ""')
                raise DistutilsSetupError("configure failed")
        # Here we run "make build" -- which builds everything except for sagelib because we
        # used configure --disable-sagelib
        # Alternative:
        # "make build-local" only builds the non-Python packages of the Sage distribution.
        # It still makes an (empty) venv in SAGE_VENV, which is unused by default;
        # but then a user could use "make build-venv" to build compatible wheels for all Python packages.
        # TODO: A target to only build wheels of tricky packages
        # (that use native libraries shared with other packages).
        SETMAKE = 'if [ -z "$MAKE" ]; then export MAKE="make -j$(PATH=build/bin:$PATH build/bin/sage-build-num-threads | cut -d" " -f 2)"; fi'
        cmd = f'cd {SAGE_ROOT} && {SETENV} && {SETMAKE} && $MAKE V=0 build-local'
        if os.system(cmd) != 0:
            raise DistutilsSetupError("make build-local failed")

        # Install configuration
        shutil.copyfile(os.path.join(SAGE_ROOT, 'pkgs', 'sage-conf', 'sage_conf.py'),
                        os.path.join(HERE, 'sage_conf.py'))
        if not self.distribution.py_modules:
            self.py_modules = self.distribution.py_modules = []
        self.distribution.py_modules.append('sage_conf')
        shutil.copyfile(os.path.join(SAGE_ROOT, 'src', 'bin', 'sage-env-config'),
                        os.path.join(HERE, 'bin', 'sage-env-config'))
        setuptools_build_py.run(self)

class build_scripts(distutils_build_scripts):

    def run(self):
        self.distribution.scripts.append(os.path.join('bin', 'sage-env-config'))
        if not self.distribution.entry_points:
            self.entry_points = self.distribution.entry_points = dict()
        if 'console_scripts' not in self.distribution.entry_points:
            self.distribution.entry_points['console_scripts'] = []
        self.distribution.entry_points['console_scripts'].append('sage-config=sage_conf:_main')
        distutils_build_scripts.run(self)

setup(
    cmdclass=dict(build_py=build_py, build_scripts=build_scripts)
)
