# -*- coding: utf-8 -*-
"""
Tests for entrypoint shebang rewriting.
"""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from sage_bootstrap.shebang import rewrite_distribution_entrypoint_shebangs


class ShebangRewriteTestCase(unittest.TestCase):
    def test_rewrite_entrypoint_shebang_to_env_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            scripts_dir = Path(tmp)
            script = scripts_dir / "jupyter"
            script.write_text(
                "#!/Users/someone/build/venv/bin/python3\nprint('ok')\n",
                encoding="utf-8",
            )
            untouched = scripts_dir / "nonpython"
            untouched.write_text("#!/bin/sh\necho ok\n", encoding="utf-8")

            fake_dist = SimpleNamespace(
                entry_points=[
                    SimpleNamespace(group="console_scripts", name="jupyter"),
                    SimpleNamespace(group="console_scripts", name="missing"),
                    SimpleNamespace(group="console_scripts", name="nonpython"),
                ]
            )
            with patch("sage_bootstrap.shebang.distribution", return_value=fake_dist):
                rewritten = rewrite_distribution_entrypoint_shebangs(
                    "jupyter-core", scripts_dir=scripts_dir
                )

            self.assertEqual(rewritten, 1)
            lines = script.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0], "#!/usr/bin/env python3")
            self.assertNotIn("/Users/", lines[0])
            self.assertEqual(
                untouched.read_text(encoding="utf-8").splitlines()[0], "#!/bin/sh"
            )
