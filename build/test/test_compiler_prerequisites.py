"""Exercise the compiler checks from configure.ac with real compilers."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(
    all(shutil.which(tool) for tool in ('autoconf', 'cc', 'c++', 'gfortran')),
    'requires Autoconf, C/C++, and gfortran compilers',
)
class CompilerPrerequisitesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='sage-compiler-checks-')
        cls.directory = Path(cls.temporary.name)
        source = (ROOT / 'configure.ac').read_text()
        checks = source.split('# Require suitable system compilers;', 1)[1]
        checks = '# Require suitable system compilers;' + checks.split(
            '\n\n###############################################################################', 1
        )[0]
        macros = (
            'ax_compiler_vendor', 'ax_cxx_compile_stdcxx',
            'ax_cxx_compile_stdcxx_11', 'ax_check_compile_flag', 'ax_openmp',
        )
        includes = '\n'.join(
            f'm4_include([{ROOT / "m4" / (name + ".m4")}])' for name in macros
        )
        (cls.directory / 'configure.ac').write_text(f'''AC_INIT([compiler-checks], [1])
m4_define([AX_REQUIRE_DEFINED], [m4_ifndef([$1], [m4_fatal([missing macro $1])])])
{includes}
m4_define([SAGE_PREREQ_URL], [See the compiler prerequisites documentation])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
host=compiler-test-linux
{checks}
AC_CONFIG_FILES([settings])
AC_OUTPUT
''')
        (cls.directory / 'settings.in').write_text(
            'FCFLAGS=@FCFLAGS@\nCXX=@CXX@\n'
            'CXX_WITHOUT_STD=@SAGE_CXX_WITHOUT_STD@\n'
            'CFLAGS_MARCH=@CFLAGS_MARCH@\n'
        )
        result = subprocess.run(['autoconf'], cwd=cls.directory, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def configure(self, versions=None, **variables):
        with tempfile.TemporaryDirectory(dir=self.directory) as temporary:
            directory = Path(temporary)
            env = os.environ.copy()
            for variable in ('CC', 'CXX', 'FC', 'CFLAGS', 'CXXFLAGS', 'FCFLAGS', 'AS', 'LD'):
                env.pop(variable, None)
            env.update(CC=shutil.which('cc'), CXX=shutil.which('c++'),
                       FC=shutil.which('gfortran'), sage_use_march_native='no')
            if versions:
                env['ax_cv_c_compiler_vendor'] = 'gnu'
                for variable, name, version in zip(
                    ('CC', 'CXX', 'FC'), ('gcc-test', 'g++-test', 'gfortran-test'), versions
                ):
                    wrapper = directory / name
                    wrapper.write_text(
                        '#!/bin/sh\nfor argument do\ncase "$argument" in\n'
                        f'  -dumpfullversion|-dumpversion) echo {shlex.quote(version)}; exit 0;;\n'
                        f'esac\ndone\nexec {shlex.quote(env[variable])} "$@"\n'
                    )
                    wrapper.chmod(0o755)
                    env[variable] = str(wrapper)
            env.update(variables)
            result = subprocess.run(
                [str(self.directory / 'configure'), f'--srcdir={self.directory}'],
                cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            settings = directory / 'settings'
            return result, settings.read_text() if settings.exists() else ''

    def test_system_compilers_and_flags(self):
        result, settings = self.configure(FCFLAGS='-O1')
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('FCFLAGS=-O1\n', settings)
        self.assertIn('CFLAGS_MARCH=\n', settings)
        self.assertIn('CXX_WITHOUT_STD=', settings)

    def test_gcc_minimum_version(self):
        result, _ = self.configure(('10.3.0', '10.3.0', '10.3.0'))
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_old_gcc(self):
        result, _ = self.configure(('10.2.1', '10.2.1', '15.2.0'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('GCC >= 10.3 is required', result.stdout)

    def test_recent_gcc(self):
        result, _ = self.configure(('17.0.0', '17.0.0', '15.2.0'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('too recent', result.stdout)

    def test_mismatched_gcc(self):
        result, _ = self.configure(('15.2.0', '14.3.0', '15.2.0'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not the same version', result.stdout)

    def test_old_gfortran(self):
        result, _ = self.configure(('15.2.0', '15.2.0', '10.2.1'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('gfortran >= 10.3 is required', result.stdout)

    def test_recent_gfortran(self):
        result, _ = self.configure(('15.2.0', '15.2.0', '17.0.0'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('too recent', result.stdout)

    def test_mismatched_assembler(self):
        result, _ = self.configure(AS=shutil.which('false'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('unset AS or set it to match', result.stdout)

    def test_broken_fortran(self):
        result, _ = self.configure(FC=shutil.which('false'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not accept free-form source code', result.stdout)


if __name__ == '__main__':
    unittest.main()
