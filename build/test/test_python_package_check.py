"""Integration tests for the system Python package configure check."""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import packaging
except ImportError:
    packaging = None


@unittest.skipUnless(shutil.which("autoconf"), "autoconf is required")
@unittest.skipUnless(packaging is not None, "packaging is required")
class PythonPackageCheckTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.root = Path(cls.workspace.name)
        source = Path(__file__).resolve().parents[2]
        (cls.root / "build" / "bin").mkdir(parents=True)
        shutil.copy2(source / "build/bin/sage-venv", cls.root / "build/bin/sage-venv")
        shutil.copy2(source / "m4/sage_python_package_check.m4", cls.root / "package-check.m4")
        (cls.root / "configure.ac").write_text("""AC_INIT([package-check-test], [1])
AC_CONFIG_SRCDIR([configure.ac])
m4_include([package-check.m4])
m4_define([SPKG_INSTALL_REQUIRES_probe], [$TEST_REQUIREMENTS])
PYTHON_FOR_VENV="$TEST_PYTHON"
enable_system_site_packages="$TEST_SYSTEM_SITE_PACKAGES"
sage_spkg_install_probe=no
sage_use_system_probe=yes
SAGE_PYTHON_PACKAGE_CHECK([probe])
AS_ECHO(["$sage_spkg_install_probe $sage_use_system_probe"]) > check-result
AC_OUTPUT
""")
        subprocess.run(["autoconf"], cwd=cls.root, check=True, capture_output=True)

    def setUp(self):
        self.fixtures = tempfile.TemporaryDirectory()
        self.addCleanup(self.fixtures.cleanup)
        self.site_packages = Path(self.fixtures.name)
        # A successful check must work even when importing pkg_resources fails.
        (self.site_packages / "pkg_resources.py").write_text(
            "raise ImportError('pkg_resources must not be used')\n")

    def install_metadata(self, name="sage-configure-probe", version="1.5",
                         requires=(), extras=()):
        directory = self.site_packages / f"{name.replace('-', '_')}-{version}.dist-info"
        directory.mkdir()
        metadata = f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n"
        metadata += "".join(f"Requires-Dist: {requirement}\n" for requirement in requires)
        metadata += "".join(f"Provides-Extra: {extra}\n" for extra in extras)
        (directory / "METADATA").write_text(metadata)

    def check(self, *requirements, install="no", system_site_packages="yes"):
        env = os.environ.copy()
        env.update(
            TEST_PYTHON=sys.executable,
            TEST_SYSTEM_SITE_PACKAGES=system_site_packages,
            TEST_REQUIREMENTS="".join(f"{requirement!r}," for requirement in requirements),
            PYTHONPATH=os.pathsep.join((str(self.site_packages),
                                       str(Path(packaging.__file__).resolve().parents[1]))),
        )
        result = subprocess.run(["sh", "configure"], cwd=self.root, env=env,
                                capture_output=True, text=True, timeout=30, check=False)
        log = (self.root / "config.log").read_text()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr + log)
        use_system = "yes" if system_site_packages == "yes" else "no"
        self.assertEqual((self.root / "check-result").read_text().strip(),
                         f"{install} {use_system}", log)
        self.assertFalse((self.root / "config.venv").exists())

    def test_version_constraints(self):
        self.install_metadata()
        for requirement, install in [
            ("sage-configure-probe", "no"),
            ("sage-configure-probe >=1.2,<2", "no"),
            ("sage-configure-probe ==1.5", "no"),
            ("sage-configure-probe >=2", "yes"),
            ("sage-configure-probe !=1.5", "yes"),
            ("sage-configure-missing >=1", "yes"),
        ]:
            with self.subTest(requirement=requirement):
                self.check(requirement, install=install)

    def test_installed_prerelease(self):
        self.install_metadata(version="1.5rc1")
        self.check("sage-configure-probe >=1,<2")

    def test_environment_markers(self):
        self.install_metadata()
        self.check("sage-configure-missing; python_version < '0'")
        self.check("sage-configure-probe >=2; python_version >= '3'", install="yes")

    def test_multiple_requirements(self):
        self.install_metadata()
        self.install_metadata(name="sage-configure-dependency", version="2")
        self.check("sage-configure-probe >=1", "sage-configure-dependency ==2")

    def test_missing_transitive_dependency(self):
        self.install_metadata(requires=("sage-configure-missing >=1",))
        self.check("sage-configure-probe", install="yes")

    def test_conflicting_transitive_dependency(self):
        self.install_metadata(requires=("sage-configure-dependency >=2",))
        self.install_metadata(name="sage-configure-dependency", version="1")
        self.check("sage-configure-probe", install="yes")

    def test_conflicting_requirement_after_dependency_was_visited(self):
        self.install_metadata()
        self.check("sage-configure-probe >=2", "sage-configure-probe >=1", install="yes")

    def test_dependency_cycle(self):
        self.install_metadata(requires=("sage-configure-dependency",))
        self.install_metadata(name="sage-configure-dependency",
                              requires=("sage-configure-probe",))
        self.check("sage-configure-probe")

    def test_transitive_markers(self):
        self.install_metadata(requires=("sage-configure-missing; python_version < '0'",))
        self.check("sage-configure-probe")

    def test_optional_dependency(self):
        self.install_metadata(requires=("sage-configure-missing; extra == 'test-extra'",),
                              extras=("test-extra",))
        self.check("sage-configure-probe")
        self.check("sage-configure-probe[test_extra]", install="yes")

    def test_satisfied_extra(self):
        self.install_metadata(requires=("sage-configure-dependency; extra == 'test-extra'",),
                              extras=("test-extra",))
        self.install_metadata(name="sage-configure-dependency")
        self.check("sage-configure-probe[test_extra]")

    def test_unknown_extra(self):
        self.install_metadata()
        self.check("sage-configure-probe[unknown]", install="yes")

    def test_base_dependency_with_extra(self):
        self.install_metadata(requires=("sage-configure-missing; extra != 'test-extra'",),
                              extras=("test-extra",))
        self.check("sage-configure-probe[test-extra]", install="yes")

    def test_missing_packaging(self):
        self.install_metadata()
        directory = self.site_packages / "packaging"
        directory.mkdir()
        (directory / "__init__.py").write_text("raise ModuleNotFoundError('packaging')\n")
        self.check("sage-configure-probe", install="yes")

    def test_system_site_packages_disabled(self):
        self.install_metadata()
        self.check("sage-configure-probe", install="yes", system_site_packages="no")
