import os
import sys
import glob
import shutil
import sysconfig
from pathlib import Path
import fnmatch

from setuptools import setup
from distutils.command.build_scripts import build_scripts as distutils_build_scripts
from setuptools.command.build_py import build_py as setuptools_build_py
from setuptools.command.egg_info import egg_info as setuptools_egg_info
from distutils.errors import (DistutilsSetupError, DistutilsModuleError,
                              DistutilsOptionError)

class build_py(setuptools_build_py):

    def initialize_options(self):
        setuptools_build_py.initialize_options(self)
        self.plat_name = None

    def finalize_options(self):
        setuptools_build_py.finalize_options(self)
        if self.plat_name is None:
            self.plat_name = self.get_finalized_command('bdist_wheel').plat_name

    def run(self):
        HERE = os.path.dirname(__file__)
        with open(os.path.join(HERE, 'VERSION.txt')) as f:
            sage_version = f.read().strip()
        SETENV = ':'

        # Have sagemath_standard (later, all packages using sage_setup.command.sage_egg_info)
        # use an install_requires of sage_conf @ the GH URL.
        cmd_bdist_wheel = self.get_finalized_command('bdist_wheel')
        # from bdist_wheel.run:
        impl_tag, abi_tag, plat_tag = cmd_bdist_wheel.get_tag()
        archive_basename = "{}-{}-{}-{}".format(cmd_bdist_wheel.wheel_dist_name, impl_tag, abi_tag, plat_tag)

        # Until pynac is repackaged as a pip-installable package (#30534), SAGE_LOCAL still has to be specific to
        # the Python version.  Note that as of pynac-0.7.26.sage-2020-04-03, on Cygwin, pynac is linked through
        # to libpython; whereas on all other platforms, it is not linked through, so we only key it to the
        # wheel tags.
        if sys.platform == 'cygwin':
            libdir_tag = sysconfig.get_config_var('LIBDIR').replace(' ', '-').replace('\\', '-').replace('/', '-')
            ldversion = sysconfig.get_config_var('LDVERSION')
            python_tag = f'{libdir_tag}-{ldversion}'
        else:
            python_tag = "{}-{}-{}".format(impl_tag, abi_tag, plat_tag)

        #
        version = self.distribution.metadata.version
        download_url = f'https://github.com/sagemath/sage-wheels/releases/download/{version}/{archive_basename}.whl'
        SETENV += f' && export SAGE_CONF_WHEEL_URL={download_url}'

        # On macOS, /var -> /private/var; we work around the DESTDIR staging bug #31569.
        STICKY = '/var/tmp'
        STICKY = str(Path(STICKY).resolve())
        # SAGE_ROOT will be a symlink during Sage runtime, but has to be a physical directory during build.
        SAGE_ROOT = os.path.join(STICKY, f'sage-{sage_version}-{python_tag}')
        # After building, we move the directory out of the way to make room for the symlink.
        # We do the wheel packaging from here.
        SAGE_ROOT_BUILD = SAGE_ROOT + '-build'
        # This will resolve via SAGE_ROOT.
        SAGE_LOCAL = os.path.join(SAGE_ROOT, 'local')
        SAGE_LOCAL_BUILD = os.path.join(SAGE_ROOT_BUILD, 'local')
        # The tree containing the wheel-building venv.  Not shipped as part of the
        # The built wheels are to be shipped separately.
        venv_name = f'venv-{python_tag}'
        SAGE_VENV = os.path.join(SAGE_ROOT, venv_name)
        SAGE_VENV_BUILD = os.path.join(SAGE_ROOT_BUILD, venv_name)
        self.SAGE_SPKG_WHEELS_BUILD = os.path.join(SAGE_VENV_BUILD,
                                                   'var', 'lib', 'sage', 'wheels')
        # Also logs
        SAGE_LOGS = os.path.join(SAGE_ROOT, 'logs')
        SAGE_LOGS_BUILD = os.path.join(SAGE_ROOT_BUILD, 'logs')

        if Path(SAGE_ROOT).is_symlink():
            # Remove symlink created by the sage_conf runtime
            os.remove(SAGE_ROOT)

        # config.status and other configure output has to be writable.
        # So (until the Sage distribution supports VPATH builds - #21469), we have to make a copy of sage_root_source.
        #
        # The file exclusions here duplicate what is done in MANIFEST.in

        spkg_with_src = ['sagelib', 'sage_conf', 'sage_docbuild', 'sage_setup', 'sage_sws2rst',
                         'sagemath_objects', 'sagemath_categories', 'sagemath_environment', 'sagemath_repl']

        def ignore(path, names):
            # exclude all embedded src trees
            if any(fnmatch.fnmatch(path, f'*/build/pkgs/{spkg}')
                   for spkg in spkg_with_src):
                return ['src']
            ### ignore more stuff --- .tox etc.
            return [name for name in names
                    if name in ('.tox',)]
        try:
            shutil.copytree(os.path.join(HERE, 'sage_root_source'), SAGE_ROOT,
                            ignore=ignore)  # will fail if already exists
        except Exception:
            print(f"hint:\n\n    To reuse a build tree left in {SAGE_ROOT} by a failed or interrupted build, first use\n\n    $ mv {SAGE_ROOT} {SAGE_ROOT_BUILD}\n")
            raise
        try:
            print(f"Created SAGE_ROOT={SAGE_ROOT}", flush=True)
            # Because of the above 'copytree',
            # within this try...finally block, SAGE_ROOT is a physical directory.

            # Re-create the "src" symlinks ignored above
            for spkg in spkg_with_src:
                shutil.copy(os.path.join(HERE, 'sage_root_source', 'build', 'pkgs', spkg, 'src'),
                            os.path.join(SAGE_ROOT, 'build', 'pkgs', spkg),
                            follow_symlinks=False)

            # Use our copy of the sage_conf template, which contains the relocation logic
            shutil.copyfile(os.path.join(HERE, '_sage_conf', '_conf.py.in'),
                            os.path.join(SAGE_ROOT, 'pkgs', 'sage-conf', '_sage_conf', '_conf.py.in'))

            if os.path.exists(SAGE_LOCAL_BUILD):
                # Previously built, start from there
                print(f"### Reusing {SAGE_LOCAL_BUILD}", flush=True)
                shutil.move(SAGE_LOCAL_BUILD, SAGE_LOCAL)

            if os.path.exists(SAGE_VENV_BUILD):
                print(f"### Reusing {SAGE_VENV_BUILD}", flush=True)
                shutil.move(SAGE_VENV_BUILD, SAGE_VENV)

            if os.path.exists(SAGE_LOGS_BUILD):
                print(f"### Reusing {SAGE_LOGS_BUILD}", flush=True)
                shutil.move(SAGE_LOGS_BUILD, SAGE_LOGS)

            # Delete old SAGE_ROOT_BUILD (if any), all the useful things have been moved from there.
            shutil.rmtree(SAGE_ROOT_BUILD, ignore_errors=True)

            CONFIGURE_ARGS = f"--prefix={SAGE_LOCAL} --with-sage-venv={SAGE_VENV} --with-python={sys.executable} --enable-build-as-root --with-system-python3=force --without-system-mpfr --without-system-readline --without-system-boost_cropped --without-system-zeromq --enable-download-from-upstream-url --enable-fat-binary --disable-notebook --disable-r --disable-doc --disable-sagelib"
            # These may be set by tox.ini
            if 'CONFIGURED_CC' in os.environ:
                CONFIGURE_ARGS += ' CC="$CONFIGURED_CC"'
            if 'CONFIGURED_CXX' in os.environ:
                CONFIGURE_ARGS += ' CXX="$CONFIGURED_CXX"'
            with open(os.path.join(SAGE_ROOT, "run-configure.sh"), "w") as f:
                f.write(f"./configure {CONFIGURE_ARGS}")
            # We run it through sage-logger... just to prefix all outputs
            cmd = fr'cd {SAGE_ROOT} && {SETENV} && V=1 build/bin/sage-logger -p "sh run-configure.sh" logs/sage-conf_relocatable.log'
            print(f"Running {cmd}", flush=True)
            if os.system(cmd) != 0:
                raise DistutilsSetupError("configure failed")

            shutil.copyfile(os.path.join(SAGE_ROOT, 'src', 'bin', 'sage-env-config'),
                            os.path.join(SAGE_ROOT, 'build', 'pkgs', 'sage_conf', 'src', 'bin', 'sage-env-config'))

            SETMAKE = 'if [ -z "$MAKE" ]; then export MAKE="make -j$(PATH=build/bin:$PATH build/bin/sage-build-num-threads | cut -d" " -f 2)"; fi'
            TARGETS = 'build'
            # We run it through sage-logger... just to prefix all outputs
            cmd = f'cd {SAGE_ROOT} && {SETENV} && {SETMAKE} && V=1 build/bin/sage-logger -p "$MAKE V=0 {TARGETS}" logs/sage-conf_relocatable.log'
            print(f"Running {cmd}", flush=True)
            if os.system(cmd) != 0:
                raise DistutilsSetupError(f"make {TARGETS} failed")

        except Exception as e:
            print(e)
            print(f"Left SAGE_ROOT={SAGE_ROOT} in place. To reuse its contents for the next build, use 'mv {SAGE_ROOT} {SAGE_ROOT_BUILD}' first", flush=True)
            raise

        else:
            shutil.move(SAGE_ROOT, SAGE_ROOT_BUILD)

        # Install configuration
        shutil.copyfile(os.path.join(SAGE_ROOT_BUILD, 'pkgs', 'sage-conf', '_sage_conf', '_conf.py'),
                        os.path.join(HERE, '_sage_conf', '_conf.py'))
        shutil.copyfile(os.path.join(SAGE_ROOT_BUILD, 'src', 'bin', 'sage-env-config'),
                        os.path.join(HERE, 'bin', 'sage-env-config'))
        # Install built SAGE_ROOT as package data
        if not self.distribution.package_data:
            self.package_data = self.distribution.package_data = {}

        # symlink into build dir
        HERE_SAGE_ROOT = os.path.join(HERE, '_sage_conf', 'sage_root')
        if os.path.islink(HERE_SAGE_ROOT):
            os.remove(HERE_SAGE_ROOT)
        os.symlink(SAGE_ROOT_BUILD, HERE_SAGE_ROOT)

        # We do not include lib64 (a symlink) because all symlinks are followed,
        # causing another copy to be installed.
        old_cwd = os.getcwd()
        os.chdir('_sage_conf')
        self.distribution.package_data['_sage_conf'] = (
            glob.glob('sage_root/*')
            + glob.glob('sage_root/config/*')
            + glob.glob('sage_root/m4/*')
            + glob.glob('sage_root/build/**', recursive=True)
            + glob.glob('sage_root/local/*')
            + glob.glob('sage_root/local/bin/**', recursive=True)
            + glob.glob('sage_root/local/include/**', recursive=True)
            + glob.glob('sage_root/local/lib/**', recursive=True)
            + glob.glob('sage_root/local/libexec/**', recursive=True)
            + glob.glob('sage_root/local/share/**', recursive=True)
            + glob.glob('sage_root/local/var/lib/**', recursive=True)  # omit /var/tmp
            )
        os.chdir(old_cwd)
        #
        setuptools_build_py.run(self)

class build_scripts(distutils_build_scripts):

    def run(self):
        self.distribution.scripts.append(os.path.join('bin', 'sage-env-config'))
        if not self.distribution.entry_points:
            self.entry_points = self.distribution.entry_points = dict()
        distutils_build_scripts.run(self)

class egg_info(setuptools_egg_info):

    def finalize_options(self):
        cmd_build_py = self.get_finalized_command('build_py')
        ## FIXME: Tried to make sure that egg_info is run _after_ build_py
        ## cmd_build_py.run()    # <-- runs it a second time, not what we want
        self.distribution.install_requires = []
        version = self.distribution.metadata.version
        download_url = f'https://github.com/sagemath/sage-wheels/releases/download/{version}'
        try:
            SAGE_SPKG_WHEELS_BUILD = cmd_build_py.SAGE_SPKG_WHEELS_BUILD
        except AttributeError:
            # egg_info invoked separately, not as part of bdist_wheel
            pass
        else:
            for f in Path(SAGE_SPKG_WHEELS_BUILD).glob("*.whl"):
                parts = str(f.stem).split('-')
                if parts[-1] != 'any':
                    # Binary wheel, include as an install_requires
                    distribution = parts[0]
                    if distribution not in ('Cython',              # has extension modules but binary compatibility is not needed
                                            'sagemath_standard',   # we are a dependency of that
                                            'tornado',             # does not need to run in the same process
                                            ):
                        self.distribution.install_requires.append(
                            f'{distribution} @ {download_url}/{f.name}')

        setuptools_egg_info.finalize_options(self)

setup(
    cmdclass=dict(build_py=build_py, build_scripts=build_scripts, egg_info=egg_info),
    # Do not mark the wheel as pure
    has_ext_modules=lambda: True
)
