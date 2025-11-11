#!/usr/bin/env python3

"""Generate some of the documentation sources."""

from __future__ import annotations

import inspect
import re
import shlex
import sys
import textwrap
from pathlib import Path
from typing import Sequence

# So that Python can find the sage_bootstrap package
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "build"))

from sage_bootstrap.expand_class import PackageClass
from sage_bootstrap.package import Package


def log_install(target_pattern: str) -> None:
    frame = inspect.currentframe()
    caller = frame.f_back if frame else None
    lineno = caller.f_lineno if caller else 0
    script = Path(__file__).name
    print(f"{script}:{lineno}: installing {target_pattern}")


def write_text(path: Path, content: str) -> None:
    """Write ``content`` to ``path``, always terminating with a newline."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if not content.endswith("\n"):
        content += "\n"
    path.write_text(content, encoding="utf-8")


def has_python_package_check(pkg) -> bool:
    """Detect whether ``pkg`` relies on ``SAGE_PYTHON_PACKAGE_CHECK``."""

    if not pkg.has_file("spkg-configure.m4"):
        return False
    spkg_configure = Path(pkg.path) / "spkg-configure.m4"
    return "SAGE_PYTHON_PACKAGE_CHECK" in spkg_configure.read_text()


SYSTEM_COMMANDS: dict[str, Sequence[str]] = {
    "arch": ("sudo", "pacman", "-S"),
    "debian": ("sudo", "apt-get", "install"),
    "fedora": ("sudo", "dnf", "install"),
    "homebrew": ("brew", "install"),
    "opensuse": ("sudo", "zypper", "install"),
    "void": ("sudo", "xbps-install"),
    "alpine": ("apk", "add"),
    "conda": ("conda", "install"),
    "freebsd": ("sudo", "pkg", "install"),
    "gentoo": ("sudo", "emerge"),
    "macports": ("sudo", "port", "install"),
    "nix": ("nix-env", "--install"),
    "openbsd": ("sudo", "pkg_add"),
    "slackware": ("sudo", "slackpkg", "install"),
}


def format_shell_command(tokens: Sequence[str], wrap: int | None = 78) -> str:
    """Format a shell command."""

    if not tokens:
        return ""
    quoted = " ".join([shlex.quote(token) for token in tokens])
    if wrap is None:
        return "$ " + quoted + " "
    lines = textwrap.wrap(
        quoted,
        width=wrap,
        initial_indent="       $ ",
        subsequent_indent="             ",
        break_long_words=False,
        break_on_hyphens=False,
    )
    lines = [line + " " for line in lines]
    return " \\\n".join(lines)


def collect_installation_packages(
    system: str,
    recommended: set[str],
    develop: set[str],
) -> dict[str, list[str]]:
    """Collect distro packages grouped by documentation category."""

    categories = {
        "standard": set(),
        "optional": set(),
        "recommended": set(),
        "develop": set(),
    }

    selector = PackageClass(":all:", has_files=[f"distros/{system}.txt"])
    for pkg_name in selector.names:
        pkg = Package(pkg_name)
        system_packages = pkg.read_system_packages(system)
        if not system_packages:
            continue
        if pkg_name == "_develop":
            categories["develop"].update(system_packages)
            continue
        if pkg_name == "_recommended":
            categories["recommended"].update(system_packages)
            continue
        if pkg.is_python_package():
            continue
        has_configure = pkg.has_file("spkg-configure.m4")
        if has_configure:
            if pkg.type == "standard":
                categories["standard"].update(system_packages)
            elif pkg_name in develop:
                categories["develop"].update(system_packages)
            elif pkg_name in recommended:
                categories["recommended"].update(system_packages)
            else:
                categories["optional"].update(system_packages)
        elif pkg_name in develop:
            categories["develop"].update(system_packages)

    return {
        key: sorted(pkg for pkg in value if pkg) for key, value in categories.items()
    }


def generate_installation_docs(base_dir: Path) -> None:
    """Generate ``installation`` command snippets."""

    base_dir.mkdir(parents=True, exist_ok=True)
    recommended = Package("_recommended").dependencies
    develop = Package("_develop").dependencies

    for system, install_command in SYSTEM_COMMANDS.items():
        log_install(f"{base_dir}/{system}*.txt")
        categories = collect_installation_packages(
            system,
            recommended,
            develop,
        )
        for suffix, key in (
            ("", "standard"),
            ("-optional", "optional"),
            ("-recommended", "recommended"),
            ("-develop", "develop"),
        ):
            path = base_dir / f"{system}{suffix}.txt"
            packages = categories.get(key, [])
            if packages:
                tokens = list(install_command) + packages
                write_text(path, format_shell_command(tokens, wrap=None))
            else:
                write_text(path, "")


def packages_for_index(
    selector: str,
    *,
    has_files: list[str] = [],
    no_files: list[str] = [],
    exclude_prefix: str | None = None,
) -> list[str]:
    """Collect package names suited for the index sections."""

    names = PackageClass(selector, has_files=has_files, no_files=no_files).names
    if exclude_prefix is not None:
        names = [name for name in names if not name.startswith(exclude_prefix)]
    return names


def write_index_sections(
    path: Path, sections: Sequence[tuple[str, Sequence[str]]]
) -> None:
    """Write an index file composed of titled sections."""

    lines: list[str] = []
    for title, packages in sections:
        lines.append(title)
        lines.append("~" * len(title))
        lines.append("")
        for name in packages:
            lines.append(f"* :ref:`spkg_{name}`")
        lines.append("")
    write_text(path, "\n".join(lines))


def write_bullet_list(path: Path, package_names: Sequence[str]) -> None:
    """Write a simple bullet list of package references."""

    lines = [f"* :ref:`spkg_{name}`" for name in package_names]
    write_text(path, "\n".join(lines))


def write_alphabetical_index(path: Path, package_names: Sequence[str]) -> None:
    """Write the alphabetical ``index_alph.rst`` file."""

    lines = [
        "",
        "Details of external packages",
        "============================",
        "",
        "Packages are in alphabetical order.",
        "",
        ".. default-role:: code",
        "",
        ".. toctree::",
        "   :maxdepth: 1",
        "",
    ]
    lines.extend(f"   {name}" for name in package_names)
    lines.extend(["", ".. default-role::", ""])
    write_text(path, "\n".join(lines))


def generate_spkg_indexes(base_dir: Path) -> None:
    """Generate the ``reference/spkg`` index files."""

    base_dir.mkdir(parents=True, exist_ok=True)
    log_install(f"{base_dir}/*.rst")

    write_index_sections(
        base_dir / "index_standard.rst",
        (
            (
                "Mathematics",
                packages_for_index(
                    ":standard:",
                    has_files=["SPKG.rst", "math"],
                    exclude_prefix="sagemath_",
                ),
            ),
            (
                "Front-end, graphics, document preparation",
                packages_for_index(
                    ":standard:",
                    has_files=["SPKG.rst", "front-end"],
                    no_files=["math"],
                    exclude_prefix="sagemath_",
                ),
            ),
            (
                "Other dependencies",
                packages_for_index(
                    ":standard:",
                    has_files=["SPKG.rst"],
                    no_files=["math", "front-end"],
                    exclude_prefix="sagemath_",
                ),
            ),
        ),
    )

    write_index_sections(
        base_dir / "index_optional.rst",
        (
            (
                "Mathematics",
                packages_for_index(
                    ":optional:",
                    has_files=["SPKG.rst", "math"],
                    exclude_prefix="sagemath_",
                ),
            ),
            (
                "Front-end, graphics, document preparation",
                packages_for_index(
                    ":optional:",
                    has_files=["SPKG.rst", "front-end"],
                    no_files=["math"],
                    exclude_prefix="sagemath_",
                ),
            ),
            (
                "Other dependencies",
                packages_for_index(
                    ":optional:",
                    has_files=["SPKG.rst"],
                    no_files=["math", "front-end"],
                    exclude_prefix="sagemath_",
                ),
            ),
        ),
    )

    sagemath_packages = [
        name
        for name in PackageClass(":all:", has_files=["SPKG.rst"]).names
        if name.startswith("sagemath_")
    ]
    write_bullet_list(base_dir / "index_sagemath.rst", sagemath_packages)

    write_index_sections(
        base_dir / "index_experimental.rst",
        (
            (
                "Mathematics",
                packages_for_index(
                    ":experimental:",
                    has_files=["SPKG.rst", "math"],
                    exclude_prefix="sagemath_",
                ),
            ),
            (
                "Other dependencies",
                packages_for_index(
                    ":experimental:",
                    has_files=["SPKG.rst"],
                    no_files=["math"],
                    exclude_prefix="sagemath_",
                ),
            ),
        ),
    )

    write_alphabetical_index(
        base_dir / "index_alph.rst",
        PackageClass(":all:", has_files=["SPKG.rst"]).names,
    )


ISSUE_RE = re.compile(r"https://github.com/sagemath/sage/issues/([0-9]+)")
ARXIV_RE = re.compile(r"https://arxiv.org/abs/cs/([0-9]+)")


def transform_spkg_rst(package_name: str, content: str) -> str:
    """Apply post-processing."""

    lines = content.splitlines()
    for idx in range(min(3, len(lines))):
        lines[idx] = re.sub(
            r"^ *Sage: Open Source Mathematics Software:",
            f"{package_name}:",
            lines[idx],
        )
    text = "\n".join(lines)
    text = ISSUE_RE.sub(r":issue:`\\1`", text)
    text = ARXIV_RE.sub(r":arxiv:`cs/\\1`", text)
    return text.rstrip() + "\n"


def _read_distro_packages(path: Path) -> list[str]:
    """Read a distro .txt file, stripping comments and blank entries."""
    if not path.is_file():
        return []
    tokens: list[str] = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            tokens.extend(line.split())
    return tokens


def build_additional_sections(pkg: Package) -> str:
    """Generate additional sections formerly produced by the shell pipeline.

    Sections:
      Type (contents of the 'type' file if present)
      Dependencies (package references)
      Version Information (best-effort from a 'version' file)
      Equivalent System Packages (commands to install system packages)
      Configuration notes (use / non-use of system packages)
    """

    lines: list[str] = []

    # Type section
    if pkg.type:
        lines.extend(
            [
                "",
                "Type",
                "----",
                "",
                pkg.type,
                "",
            ]
        )

    # Dependencies section
    lines.extend(
        [
            "",
            "Dependencies",
            "------------",
            "",
        ]
    )
    dependencies = pkg.dependencies + pkg.dependencies_order_only
    if dependencies:
        for dep in sorted(dependencies):
            if dep:
                if dep.startswith("$") or dep.startswith("sage"):
                    lines.append(f"- {dep}")
                else:
                    lines.append(f"- :ref:`spkg_{dep}`")
        lines.append("")
    else:
        lines.append("")

    # Version Information section (heuristic)
    lines.extend(["Version Information", "-------------------", ""])
    for candidate in (
        "package-version.txt",
        "requirements.txt",
        "pyproject.toml",
        "version_requirements.txt",
    ):
        path = Path(pkg.path) / candidate
        if path.is_file():
            version_text = path.read_text().strip()
            if version_text:
                lines.extend(
                    [
                        f"{candidate}::",
                        "",
                        *(
                            "    " + line
                            for line in version_text.splitlines()
                            if not line.startswith("#")
                        ),
                        "",
                    ]
                )
    lines.append("")

    # Equivalent System Packages
    lines.extend(
        [
            "Equivalent System Packages",
            "--------------------------",
            "",
        ]
    )
    distros_dir = Path(pkg.path) / "distros"
    system_files: list[Path] = []
    if distros_dir.is_dir():
        system_files = sorted(p for p in distros_dir.glob("*.txt"))
    have_any = False
    have_repology = any(p.stem == "repology" for p in system_files)
    for p in system_files:
        system = p.stem
        if system == "repology":
            continue  # defer
        packages = _read_distro_packages(p)
        # Heading for system (simulate tab title)
        pretty = {
            "alpine": "Alpine",
            "arch": "Arch Linux",
            "conda": "conda-forge",
            "debian": "Debian/Ubuntu",
            "fedora": "Fedora/Redhat/CentOS",
            "freebsd": "FreeBSD",
            "gentoo": "Gentoo Linux",
            "homebrew": "Homebrew",
            "macports": "MacPorts",
            "nix": "Nixpkgs",
            "openbsd": "OpenBSD",
            "opensuse": "openSUSE",
            "slackware": "Slackware",
            "void": "Void Linux",
        }.get(system, system)
        if packages:
            if system == "pyodide":
                lines.extend(
                    [
                        f".. tab:: {pretty}",
                        "",
                        f"   install the following packages: {', '.join(packages)}",
                        "",
                    ]
                )
            else:
                lines.extend([f".. tab:: {pretty}", "", "   .. CODE-BLOCK:: bash"])
                have_any = True
                cmd_tokens = list(SYSTEM_COMMANDS.get(system, ())) + packages
                if cmd_tokens:
                    lines.append("")
                    lines.append(format_shell_command(cmd_tokens))
                    lines.append("")
                    lines.append("")
                else:
                    lines.append("(packages: " + " ".join(packages) + ")")
                    lines.append("")
        else:
            lines.append(f".. tab:: {pretty}")
            lines.append("")
            lines.append("   No package needed.")
            lines.append("")

    # Repology shown after others
    if have_repology:
        repology_file = distros_dir / "repology.txt"
        packages = _read_distro_packages(repology_file)
        if packages:
            lines.append("")
            urls = ", ".join(
                f"https://repology.org/project/{p}/versions" for p in packages
            )
            lines.append(f"See {urls}")
            lines.append("")
            have_any = True
    if not have_any:
        lines.append("(none known)")
        lines.append("")

    # Configuration notes
    if pkg.has_file("spkg-configure.m4"):
        if has_python_package_check(pkg):
            lines.extend(
                [
                    "If the system package is installed and if the (experimental) option",
                    "``--enable-system-site-packages`` is passed to ``./configure``, then ``./configure``",
                    "will check if the system package can be used.",
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    "If the system package is installed, ``./configure`` will check if it can be used.",
                    "",
                ]
            )
    elif not pkg.name.startswith("_"):
        lines.extend(
            [
                "However, these system packages will not be used for building Sage",
                "because ``spkg-configure.m4`` has not been written for this package;",
                "see :issue:`27330` for more information.",
                "",
            ]
        )

    return "\n".join(lines).rstrip() + "\n\n"


def generate_spkg_details(base_dir: Path) -> None:
    """Generate per-package ``spkg`` documentation."""

    base_dir.mkdir(parents=True, exist_ok=True)
    log_install(f"{base_dir}/*.rst")
    names = PackageClass(":all:", has_files=["SPKG.rst"]).names

    for name in names:
        pkg = Package(name)
        source = Path(pkg.path) / "SPKG.rst"
        if not source.is_file():
            continue
        transformed = transform_spkg_rst(name, source.read_text())
        extra = build_additional_sections(pkg)
        content = f".. _spkg_{name}:\n\n{transformed}{extra}"
        write_text(base_dir / f"{name}.rst", content)


def main() -> int:
    target_dir = Path(sys.argv[1])

    generate_installation_docs(target_dir / "en" / "installation")
    generate_spkg_indexes(target_dir / "en" / "reference" / "spkg")
    generate_spkg_details(target_dir / "en" / "reference" / "spkg")

    return 0


if __name__ == "__main__":
    sys.exit(main())
