"""Test the terminal-library configure checks with isolated headers/libraries."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(
    all(shutil.which(tool) for tool in ('aclocal', 'autoconf', 'cc', 'ar', 'pkg-config')),
    'requires Autotools, a C toolchain, and pkg-config',
)
class TerminalPrerequisitesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='sage-terminal-checks-')
        cls.directory = Path(cls.temporary.name)
        source = (ROOT / 'configure.ac').read_text()
        checks = source.split('sage_terminal_save_CPPFLAGS=', 1)[1]
        checks = 'sage_terminal_save_CPPFLAGS=' + checks.split(
            '\n# Check that we are not building in a directory containing spaces', 1
        )[0]
        (cls.directory / 'configure.ac').write_text(f'''AC_INIT([terminal-checks], [1])
AC_PROG_CC
PKG_PROG_PKG_CONFIG
m4_define([SAGE_PREREQ_URL], [See terminal prerequisites])
{checks}
AC_CONFIG_FILES([settings])
AC_OUTPUT
''')
        (cls.directory / 'settings.in').write_text(
            'CPPFLAGS=@CPPFLAGS@\nLIBS=@LIBS@\n'
            'READLINE_LIBS=@READLINE_LIBS@\nNCURSES_LIBS=@NCURSES_LIBS@\n'
            'PREFIX=@SAGE_READLINE_PREFIX@\n'
        )
        for command in (['aclocal'], ['autoconf']):
            result = subprocess.run(command, cwd=cls.directory, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def configure(self, wide=False, pkg_config=True, readline_api=True,
                  readline_header=True, **variables):
        with tempfile.TemporaryDirectory(dir=self.directory) as temporary:
            directory = Path(temporary)
            include = directory / 'include'
            library = directory / 'lib'
            pc = library / 'pkgconfig'
            (include / 'readline').mkdir(parents=True)
            pc.mkdir(parents=True)
            (include / 'ncurses.h').write_text('int wresize(void *, int, int);\n')
            (include / 'readline/readline.h').write_text(
                'int rl_bind_keyseq(const char *, void *);\n' if readline_header else '#error missing header\n'
            )
            curses = 'ncursesw' if wide else 'ncurses'
            for name, contents in (
                (curses, 'int wresize(void *p, int h, int w) { return 0; }\n'),
                ('readline', 'extern int wresize(void *, int, int);\n'
                 'int rl_bind_keyseq(const char *s, void *p) { return wresize(p, 0, 0); }\n'
                 if readline_api else 'int readline_compatibility_only(void) { return 0; }\n'),
            ):
                source = directory / (name + '.c')
                source.write_text(contents)
                obj = directory / (name + '.o')
                subprocess.run(['cc', '-c', str(source), '-o', str(obj)], check=True, capture_output=True)
                subprocess.run(['ar', 'rcs', str(library / ('lib' + name + '.a')), str(obj)],
                               check=True, capture_output=True)
                if pkg_config:
                    (pc / (name + '.pc')).write_text(f'''prefix={directory}
Name: {name}
Description: Isolated test library
Version: 8.0
Libs: -L{library} -l{name}
Cflags: -I{include}
''')
            env = os.environ.copy()
            for variable in ('CC', 'CFLAGS', 'CPPFLAGS', 'LDFLAGS', 'LIBS', 'PKG_CONFIG',
                             'READLINE_CFLAGS', 'READLINE_LIBS', 'NCURSES_CFLAGS', 'NCURSES_LIBS',
                             'NCURSESW_CFLAGS', 'NCURSESW_LIBS', 'SAGE_READLINE_PREFIX'):
                env.pop(variable, None)
            env.update(CPPFLAGS=f'-I{include} -DPRESERVED_FLAG', LDFLAGS=f'-L{library}',
                       LIBS='-lm', PKG_CONFIG_PATH='', PKG_CONFIG_LIBDIR=str(pc))
            env.update(variables)
            result = subprocess.run(
                [str(self.directory / 'configure'), f'--srcdir={self.directory}'],
                cwd=directory, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            settings = directory / 'settings'
            return result, settings.read_text() if settings.exists() else ''

    def test_pkg_config_and_flag_restoration(self):
        result, settings = self.configure()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('LIBS=-lm\n', settings)
        self.assertIn('-DPRESERVED_FLAG\n', settings)
        self.assertIn('-lreadline', settings)
        self.assertIn('-lncurses', settings)
        self.assertNotIn('PREFIX=\n', settings)

    def test_wide_ncurses(self):
        result, settings = self.configure(wide=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('-lncursesw', settings)

    def test_without_pkg_config_metadata(self):
        result, settings = self.configure(pkg_config=False)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('PREFIX=\n', settings)
        self.assertIn('READLINE_LIBS=-lreadline -lncurses', settings)

    def test_unusable_ncurses_metadata(self):
        result, _ = self.configure(NCURSES_LIBS='-lsage_missing_ncurses')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('could not link against ncurses', result.stdout)

    def test_readline_compatibility_library(self):
        result, _ = self.configure(readline_api=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('GNU readline with rl_bind_keyseq support', result.stdout)

    def test_missing_readline_header(self):
        result, _ = self.configure(readline_header=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('could not find readline/readline.h', result.stdout)


if __name__ == '__main__':
    unittest.main()
