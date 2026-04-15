#!/usr/bin/env python3
"""
Helpers for normalizing Python entrypoint shebangs.
"""

from __future__ import annotations

import argparse
import os
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path
import sysconfig


ENV_PYTHON3_SHEBANG = "#!/usr/bin/env python3"


def _entrypoint_script_names(dist_name: str) -> set[str]:
    try:
        dist = distribution(dist_name)
    except PackageNotFoundError:
        return set()
    return {
        ep.name
        for ep in dist.entry_points
        if ep.group in ("console_scripts", "gui_scripts")
    }


def rewrite_distribution_entrypoint_shebangs(
    dist_name: str,
    scripts_dir: str | os.PathLike[str] | None = None,
    shebang: str = ENV_PYTHON3_SHEBANG,
) -> int:
    """
    Rewrite shebangs of a distribution's generated entrypoint scripts.

    Only scripts with a Python shebang are modified.
    """
    scripts_path = Path(scripts_dir) if scripts_dir else Path(sysconfig.get_path("scripts"))
    rewritten = 0
    desired = shebang.encode("utf-8")
    for script_name in _entrypoint_script_names(dist_name):
        script_path = scripts_path / script_name
        if not script_path.is_file():
            continue
        content = script_path.read_bytes()
        if not content:
            continue
        lines = content.splitlines(keepends=True)
        first = lines[0].rstrip(b"\r\n")
        if not first.startswith(b"#!") or b"python" not in first.lower():
            continue
        if first == desired:
            continue
        newline = b"\n"
        if lines[0].endswith(b"\r\n"):
            newline = b"\r\n"
        lines[0] = desired + newline
        script_path.write_bytes(b"".join(lines))
        rewritten += 1
    return rewritten


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rewrite Python entrypoint shebangs for one distribution.",
    )
    parser.add_argument("--distribution", required=True, help="Installed distribution name")
    parser.add_argument(
        "--scripts-dir",
        default=None,
        help="Directory containing generated scripts (default: interpreter scripts dir)",
    )
    parser.add_argument(
        "--shebang",
        default=ENV_PYTHON3_SHEBANG,
        help=f"Replacement shebang (default: {ENV_PYTHON3_SHEBANG})",
    )
    args = parser.parse_args()
    rewrite_distribution_entrypoint_shebangs(
        dist_name=args.distribution,
        scripts_dir=args.scripts_dir,
        shebang=args.shebang,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
