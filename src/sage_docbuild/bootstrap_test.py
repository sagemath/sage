"""Regression tests for the standalone documentation-source generator."""

import errno
import hashlib
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SAGE_ROOT = Path(__file__).resolve().parents[2]


def load_bootstrap_docs():
    """Load ``tools/bootstrap-docs.py`` despite the hyphen in its name."""

    name = "_sage_bootstrap_docs_test"
    path = SAGE_ROOT / "tools" / "bootstrap-docs.py"
    if not path.is_file():
        raise unittest.SkipTest(f"standalone generator is unavailable: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    old_path = sys.path[:]
    old_dont_write_bytecode = sys.dont_write_bytecode
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    finally:
        sys.path[:] = old_path
        sys.dont_write_bytecode = old_dont_write_bytecode
    return module


bootstrap_docs = load_bootstrap_docs()


class BodyError(OSError):
    """An error from inside a lock body, distinct from a locking error."""


@unittest.skipIf(bootstrap_docs.fcntl is None, "fcntl is unavailable")
class GenerationLockTest(unittest.TestCase):
    def assert_body_survives_lock_failure(self, lock):
        ran = False
        with lock():
            ran = True
        self.assertTrue(ran)

        with self.assertRaises(BodyError):
            with lock():
                raise BodyError("body failed")

    def test_open_failure_runs_unlocked_without_hiding_body_error(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)

            def lock():
                return bootstrap_docs._generation_lock(target, exclusive=True)

            error = OSError(errno.EIO, "lock file unavailable")
            with mock.patch.object(bootstrap_docs.os, "open", side_effect=error):
                self.assert_body_survives_lock_failure(lock)

    def test_flock_failure_runs_unlocked_without_hiding_body_error(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)

            def lock():
                return bootstrap_docs._generation_lock(target, exclusive=True)

            error = OSError(errno.ENOLCK, "no locks available")
            with mock.patch.object(
                bootstrap_docs.fcntl, "flock", side_effect=error
            ):
                self.assert_body_survives_lock_failure(lock)

    def test_unlock_failure_does_not_replace_body_result_or_error(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)

            def flock(file, operation):
                if operation == bootstrap_docs.fcntl.LOCK_UN:
                    raise OSError(errno.ENOLCK, "no locks available")

            def lock():
                return bootstrap_docs._generation_lock(target, exclusive=True)

            with mock.patch.object(bootstrap_docs.fcntl, "flock", side_effect=flock):
                self.assert_body_survives_lock_failure(lock)


class GeneratedOutputCheckTest(unittest.TestCase):
    def complaints_for(self, target, output, *, reader=None):
        content = output.read_bytes()
        manifest = {
            "version": bootstrap_docs.MANIFEST_VERSION,
            "inputs": "inputs",
            "outputs": {
                output.relative_to(target).as_posix(): hashlib.sha256(
                    content
                ).hexdigest()
            },
        }
        patches = [
            mock.patch.object(bootstrap_docs, "input_digest", return_value="inputs"),
            mock.patch.object(
                bootstrap_docs, "expected_targets", return_value=[output]
            ),
            mock.patch.object(bootstrap_docs, "_read_manifest", return_value=manifest),
            mock.patch.object(
                bootstrap_docs, "_generated_extras", return_value=([], [])
            ),
        ]
        if reader is not None:
            patches.append(
                mock.patch.object(
                    bootstrap_docs,
                    "_read_regular_file_stably",
                    side_effect=reader,
                )
            )
        with patches[0], patches[1], patches[2], patches[3]:
            if reader is None:
                return bootstrap_docs.complaints_about(target)
            with patches[4]:
                return bootstrap_docs.complaints_about(target)

    @unittest.skipIf(os.name == "nt", "POSIX permission bits are required")
    def test_group_readable_output_is_accepted_independently_of_umask(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            output = target / "generated.rst"
            output.write_text("generated\n", encoding="utf-8")
            output.chmod(0o640)

            # Under the old check, this simulated umask made the readable 0640
            # file fail because it did not also grant read access to "other".
            with mock.patch.object(bootstrap_docs, "_file_mode", return_value=0o644):
                self.assertEqual(self.complaints_for(target, output), [])

    def test_stable_read_failure_is_reported_as_unreadable(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            output = target / "generated.rst"
            output.write_text("generated\n", encoding="utf-8")

            complaints = self.complaints_for(
                target,
                output,
                reader=OSError("file moved while it was read"),
            )
            self.assertEqual(
                complaints,
                [f"unreadable: {output}: file moved while it was read"],
            )


if __name__ == "__main__":
    unittest.main()
