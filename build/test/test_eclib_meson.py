"""Exercise eclib selection without configuring all of Sage's dependencies.

The system pkg-config metadata and fallback are fixtures; the dependency call
and wrap metadata are taken from the Sage build. Building and testing the real
fallback additionally requires FLINT, NTL and PARI.
"""

import configparser
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


@unittest.skipUnless(shutil.which('meson') and shutil.which('pkg-config'),
                     'meson and pkg-config are required')
class EclibMesonTestCase(unittest.TestCase):
    def check_selection(self, system_version, expected, *, option='auto',
                        arguments=(), windows=False, succeeds=True):
        source = Path(__file__).resolve().parents[2]
        dependency = re.search(
            r'^ec = dependency\(\n.*?^\)',
            (source / 'src/meson.build').read_text(), re.M | re.S,
        ).group()
        wrap = configparser.ConfigParser()
        wrap.read(source / 'subprojects/eclib.wrap')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pc = root / 'pkgconfig'
            pc.mkdir()
            if system_version is not None:
                (pc / 'eclib.pc').write_text(
                    'Name: eclib\nDescription: eclib selection fixture\n'
                    f'Version: {system_version}\nLibs:\nCflags:\n')
            fallback = root / 'subprojects' / wrap['wrap-file']['directory']
            fallback.mkdir(parents=True)
            shutil.copy2(source / 'subprojects/eclib.wrap', root / 'subprojects')
            (fallback / 'meson.build').write_text("""
project('eclib', version: '20250627')
eclib_dep = declare_dependency(version: meson.project_version())
meson.override_dependency('eclib', eclib_dep)
""")
            (root / 'meson.options').write_text(
                "option('eclib', type: 'feature', value: 'auto')\n")
            (root / 'meson.build').write_text(
                "project('eclib-selection', 'cpp', meson_version: '>=1.11.2')\n"
                f"is_windows = {'true' if windows else 'false'}\n"
                + dependency + """
if ec.found()
  message('SELECTED: ' + ec.type_name() + ' ' + ec.version())
else
  message('SELECTED: disabled')
endif
""")
            env = os.environ.copy()
            env['PKG_CONFIG_LIBDIR'] = str(pc)
            env['PKG_CONFIG_PATH'] = ''
            result = subprocess.run(
                ['meson', 'setup', str(root / 'build'), str(root),
                 f'-Declib={option}', '--wrap-mode=nodownload', *arguments],
                env=env, capture_output=True, text=True, timeout=60, check=False,
            )
            output = result.stdout + result.stderr
            self.assertEqual(result.returncode == 0, succeeds, output)
            if succeeds:
                self.assertIn('SELECTED: ' + expected, output)

    def test_compatible_system(self):
        # Include the old version reported by some compatible conda packages.
        for version in ('20231211', '20250122', '20250627'):
            with self.subTest(version=version):
                self.check_selection(version, 'pkgconfig ' + version)

    def test_incompatible_or_missing_system(self):
        for version in ('20231210', '20260422', '20260707', None):
            with self.subTest(version=version):
                self.check_selection(version, 'internal 20250627')

    def test_disabled(self):
        for version in ('20250627', '20260707', None):
            with self.subTest(version=version):
                self.check_selection(version, 'disabled', option='disabled')

    def test_enabled(self):
        self.check_selection('20260707', 'internal 20250627', option='enabled')

    def test_forced_fallback(self):
        self.check_selection('20250627', 'internal 20250627',
                             arguments=('--force-fallback-for=eclib',))

    def test_no_fallback(self):
        self.check_selection('20260707', 'disabled',
                             arguments=('--wrap-mode=nofallback',))
        self.check_selection('20260707', None, option='enabled', succeeds=False,
                             arguments=('--wrap-mode=nofallback',))

    def test_no_automatic_fallback_on_windows(self):
        self.check_selection(None, 'disabled', windows=True)


if __name__ == '__main__':
    unittest.main()
