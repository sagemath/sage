"""Regression tests for build-time discovery of shared bibliographies."""

import runpy
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


COPY = runpy.run_path(
    str(Path(__file__).resolve().parents[2] / "tools/copy-bibliographies.py")
)["copy_bibliographies"]


class CopyBibliographiesTest(unittest.TestCase):
    @unittest.skipUnless(shutil.which("meson"), "Meson is not installed")
    def test_meson_discovers_changes_without_reconfiguration(self):
        repo = Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory(prefix="meson bibliography ") as directory:
            root = Path(directory)
            source, build = root / "source", root / "build"
            refs = source / "src/doc/references"
            refs.mkdir(parents=True)
            (source / "tools").mkdir()
            shutil.copyfile(
                repo / "tools/copy-bibliographies.py",
                source / "tools/copy-bibliographies.py",
            )
            shutil.copyfile(repo / "src/doc/references/meson.build", refs / "meson.build")
            (source / "meson.build").write_text('''project('bibliography-discovery')
py = import('python').find_installation()
doc_src = []
subdir('src/doc/references')
custom_target('document', output: 'document-count',
  command: [py, '-c', 'from pathlib import Path; import sys; p = Path(sys.argv[1]); p.write_text(str(int(p.read_text()) + 1) if p.exists() else "1")', '@OUTPUT@'],
  depends: doc_bibliographies,
  build_by_default: true)
''')

            def run(*args):
                result = subprocess.run(args, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            def compile_and_check(expected):
                run("meson", "compile", "-C", str(build))
                self.assertEqual((build / "document-count").read_text(), str(expected))
                self.assertEqual(
                    {p.name: p.read_bytes() for p in (build / "src/doc/references").glob("*.bib")},
                    {p.name: p.read_bytes() for p in refs.glob("*.bib")},
                )

            first = refs / "reference-database.bib"
            first.write_text("@misc{First, title={First}}\n")
            run("meson", "setup", str(build), str(source))
            compile_and_check(1)
            compile_and_check(1)
            added = refs / "PR123.bib"
            added.write_text("@misc{Added, title={Added}}\n")
            compile_and_check(2)
            first.write_text("@misc{First, title={Edited}}\n")
            compile_and_check(3)
            added.rename(refs / "PR124.bib")
            compile_and_check(4)
            (refs / "PR124.bib").unlink()
            compile_and_check(5)
            compile_and_check(5)
            (build / "src/doc/references/reference-database.bib").unlink()
            compile_and_check(6)

    def test_add_edit_remove_and_no_op(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source with spaces"
            destination = root / "build with spaces"
            source.mkdir()
            stamp = destination / "bibliographies.stamp"
            first = source / "reference-database.bib"
            first.write_text("@misc{First, title={First}}\n")
            COPY(source, destination, stamp)
            copied = destination / first.name
            self.assertEqual(copied.read_bytes(), first.read_bytes())
            original_times = (copied.stat().st_mtime_ns, stamp.stat().st_mtime_ns)
            COPY(source, destination, stamp)
            self.assertEqual(
                (copied.stat().st_mtime_ns, stamp.stat().st_mtime_ns), original_times
            )
            added = source / "PR123.bib"
            added.write_text("@misc{Added, title={Added}}\n")
            COPY(source, destination, stamp)
            self.assertEqual((destination / added.name).read_bytes(), added.read_bytes())
            first.write_text("@misc{First, title={Edited}}\n")
            COPY(source, destination, stamp)
            self.assertEqual(copied.read_bytes(), first.read_bytes())
            unrelated = destination / "keep.txt"
            unrelated.write_text("keep")
            added.unlink()
            COPY(source, destination, stamp)
            self.assertFalse((destination / added.name).exists())
            self.assertEqual(unrelated.read_text(), "keep")

    def test_reject_source_as_destination(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(ValueError):
                COPY(root, root, root / "stamp")

    def test_missing_source_preserves_existing_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old = root / "PR123.bib"
            old.write_text("keep")
            with self.assertRaises(FileNotFoundError):
                COPY(root / "missing", root, root / "stamp")
            self.assertEqual(old.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
