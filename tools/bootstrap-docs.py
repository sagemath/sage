#!/usr/bin/env python3

"""Generate some of the documentation sources."""

from __future__ import annotations

import functools
import hashlib
import inspect
import json
import os
import re
import shlex
import stat
import sys
import tempfile
import textwrap
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

try:
    import fcntl
except ModuleNotFoundError:
    # Windows has no fcntl, and Sage does not build there.  This generator
    # runs as a program of its own, with only <root>/build on its path, so it
    # cannot share sage_docbuild's copy of this and states it again; see
    # _generation_lock().
    fcntl = None

if TYPE_CHECKING:
    from collections.abc import Sequence

# Bytecode caches are deliberately excluded from the generator inputs.  This
# is a dedicated process, so prevent its imports from changing the metadata of
# an otherwise stable input directory between the before/after snapshots.
sys.dont_write_bytecode = True

# So that Python can find the sage_bootstrap package
sys.path.insert(0, str(Path(__file__).parent.parent / "build"))

from sage_bootstrap.env import SAGE_ROOT
from sage_bootstrap.expand_class import PackageClass
from sage_bootstrap.package import Package

MANIFEST_VERSION = 1
MANIFEST_NAME = "sage-docbuild-manifest.json"
PACKAGE_INPUT_NAMES = (
    "SPKG.rst",
    "dependencies",
    "dependencies_order_only",
    "front-end",
    "math",
    "package-version.txt",
    "pyproject.toml",
    "requirements.txt",
    "spkg-configure.m4",
    "type",
    "version_requirements.txt",
)


def _stat_identity(info) -> tuple[int, ...]:
    """Return metadata that changes when a filesystem object is replaced.

    The content digest stored in the manifest deliberately ignores metadata,
    but a second, short-lived token is needed while generating: contents can
    change from A to B and back to A while a writer observes B.  ``ctime`` is
    included because an ordinary process cannot restore it after such a
    change.
    """

    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_uid,
        info.st_gid,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _add_identity(digest, label: str, kind: str, info) -> None:
    """Add a labelled filesystem identity to a transient race digest."""

    digest.update(label.encode("utf-8", errors="surrogateescape"))
    digest.update(b"\0")
    digest.update(kind.encode("ascii"))
    digest.update(b"\0")
    digest.update(repr(_stat_identity(info)).encode("ascii"))
    digest.update(b"\0")


def _walk_tree(root: Path):
    """Yield a tree without suppressing directory traversal errors.

    ``Path.rglob`` follows the modern ``Path.is_*`` convention of treating
    every ``OSError`` as a negative result.  Package metadata is an input, so
    an unreadable directory must abort generation rather than look absent.
    Directory symlinks are yielded but, like ``Path.rglob``'s default, are not
    followed.
    """

    for path in sorted(root.iterdir()):
        info = path.lstat()
        if stat.S_ISDIR(info.st_mode) and path.name == "__pycache__":
            continue
        yield path
        if stat.S_ISDIR(info.st_mode):
            yield from _walk_tree(path)


@contextmanager
def _generation_lock(target_dir: Path, *, exclusive: bool):
    """Serialize writers and share their lock when a check can do so read-only.

    A check must remain usable on an installed or shared read-only source tree.
    It therefore opens an existing lock read-only, and runs unlocked when no
    lock exists; a false stale result during that narrow race merely makes the
    caller wait for the writer and regenerate unchanged files atomically.

    Where the platform has no :mod:`fcntl`, or where the filesystem refuses
    opening or locking the lock file, the lock is absent rather than emulated.
    Generated files are still published atomically and obsolete-file ownership
    is rechecked immediately before removal, so an unavailable advisory lock
    must not make an otherwise usable documentation tree fail to build.
    """

    if fcntl is None:
        yield
        return
    lock_path = manifest_path(target_dir).with_suffix(".lock")
    try:
        if exclusive:
            lock_path.parent.mkdir(parents=True, exist_ok=True)
            descriptor = os.open(lock_path, os.O_RDONLY | os.O_CREAT, 0o666)
        else:
            descriptor = os.open(lock_path, os.O_RDONLY)
    except OSError:
        # The lock is an optimization, not a reason to reject an otherwise
        # usable tree.  In particular, network and read-only filesystems can
        # refuse either creating or opening the persistent lock file.
        yield
        return
    with os.fdopen(descriptor, "rb") as lock:
        operation = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        try:
            fcntl.flock(lock, operation)
        except OSError:
            # Some network filesystems provide fcntl but not advisory locks.
            # The writes themselves are atomic, so retain the pre-lock behavior
            # instead of making documentation builds on those filesystems fail.
            yield
            return
        try:
            yield
        finally:
            # Closing the descriptor releases the lock in any case.  Do not
            # let a refused explicit unlock hide an exception from the body or
            # turn a completed generation into a failure.
            try:
                fcntl.flock(lock, fcntl.LOCK_UN)
            except OSError:
                pass


@functools.cache
def _file_mode() -> int:
    """Return the ordinary creation mode selected by the process umask.

    Probing costs a temporary directory and half a dozen system calls, and
    every one of the hundreds of generated files asks.  Nothing that runs
    while this program generates changes the umask, so the first answer
    stands for the whole run.
    """

    directory = Path(tempfile.mkdtemp(prefix="sage-docbuild-mode-"))
    descriptor = None
    probe = directory / "probe"
    try:
        descriptor = os.open(
            probe, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666
        )
        return stat.S_IMODE(os.fstat(descriptor).st_mode)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        probe.unlink(missing_ok=True)
        directory.rmdir()


def _input_status(path: Path):
    """Return the followed status of an input, or ``None`` only if absent."""

    try:
        path.lstat()
    except FileNotFoundError:
        return None
    try:
        return path.stat()
    except FileNotFoundError as error:
        raise OSError(f"broken or moving input path: {path}") from error


def _input_is_file(path: Path) -> bool:
    """Test an input file without hiding permission or I/O errors."""

    info = _input_status(path)
    return info is not None and stat.S_ISREG(info.st_mode)


def _input_is_dir(path: Path) -> bool:
    """Test an input directory without hiding permission or I/O errors."""

    info = _input_status(path)
    return info is not None and stat.S_ISDIR(info.st_mode)


def _package_directories() -> list[Path]:
    """Return package directories, validating every visible ``type`` input."""

    packages = Path(SAGE_ROOT) / "build" / "pkgs"
    if not _input_is_dir(packages):
        raise FileNotFoundError(
            f"package metadata directory does not exist: {packages}"
        )
    answer = []
    for package in sorted(packages.iterdir()):
        try:
            package_info = package.stat()
        except FileNotFoundError as error:
            raise OSError(
                f"package metadata changed while it was read: {package}"
            ) from error
        if not stat.S_ISDIR(package_info.st_mode):
            continue
        type_path = package / "type"
        try:
            type_link = type_path.lstat()
        except FileNotFoundError:
            continue  # a directory in build/pkgs, but not a Sage package
        type_info = type_path.stat()
        if not stat.S_ISREG(type_info.st_mode):
            raise OSError(f"package type metadata is not a regular file: {type_path}")
        # Keep the lstat read meaningful: a type symlink that changes between
        # lstat and stat is a moving input and the later digest check will fail.
        _ = type_link
        answer.append(package)
    return answer


def _documented_package_names() -> list[str]:
    """Return packages with a readable regular ``SPKG.rst`` input."""

    names = []
    for package in _package_directories():
        source = package / "SPKG.rst"
        try:
            source.lstat()
        except FileNotFoundError:
            continue
        if not _input_is_file(source):
            raise OSError(f"package documentation is not a regular file: {source}")
        names.append(package.name)
    return names


def log_install(target_pattern: str) -> None:
    frame = inspect.currentframe()
    caller = frame.f_back if frame else None
    lineno = caller.f_lineno if caller else 0
    script = Path(__file__).name
    print(f"{script}:{lineno}: installing {target_pattern}")


def write_text(path: Path, content: str) -> None:
    """Write ``content`` to ``path``, always terminating with a newline.

    The file appears whole or not at all: an interrupted run would otherwise
    leave a truncated page behind, which every later run takes for a file it
    has already written.
    """

    if not content.endswith("\n"):
        content += "\n"
    encoded = content.encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)

    # Path.write_text(), which this replaces, creates a file according to the
    # process umask.  mkstemp deliberately creates mode 0600, so set the mode
    # explicitly before publishing its file.  An unchanged page is left in
    # place: Sphinx uses its mtime for incremental builds.
    mode = _file_mode()
    try:
        info = path.lstat()
    except FileNotFoundError:
        info = None
    if (info is not None
            and not (stat.S_ISREG(info.st_mode)
                     or stat.S_ISLNK(info.st_mode))):
        raise OSError(f"refusing to replace non-file generated path: {path}")
    if info is not None and stat.S_ISREG(info.st_mode):
        unchanged = path.read_bytes() == encoded
        if unchanged:
            if stat.S_IMODE(info.st_mode) != mode:
                path.chmod(mode)
            return

    handle, temporary = tempfile.mkstemp(dir=path.parent, prefix=path.name,
                                         suffix=".tmp")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as file:
            file.write(content)
            os.fchmod(file.fileno(), mode)
        # os.replace() atomically replaces a regular file or a symlink.  A
        # directory or special file cannot be a generated artifact and is
        # preserved rather than recursively removed.
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def remove_path(path: Path) -> None:
    """Remove a generated file or link, refusing directories and special files."""

    try:
        info = path.lstat()
    except FileNotFoundError:
        return
    if not (stat.S_ISREG(info.st_mode) or stat.S_ISLNK(info.st_mode)):
        raise OSError(f"refusing to remove non-file generated path: {path}")
    path.unlink(missing_ok=True)


def has_python_package_check(pkg) -> bool:
    """Detect whether ``pkg`` relies on ``SAGE_PYTHON_PACKAGE_CHECK``."""

    if not pkg.has_file("spkg-configure.m4"):
        return False
    spkg_configure = Path(pkg.path) / "spkg-configure.m4"
    return "SAGE_PYTHON_PACKAGE_CHECK" in spkg_configure.read_text()


# The suffix that each category of packages is written under, next to the name
# of the system.  Both the generator and :func:`expected_targets` read this, so
# that the list of files written and the list of files looked for cannot drift
# apart.
INSTALLATION_CATEGORIES: Sequence[tuple[str, str]] = (
    ("", "standard"),
    ("-optional", "optional"),
    ("-recommended", "recommended"),
    ("-develop", "develop"),
)

# The index files that :func:`generate_spkg_indexes` writes.
SPKG_INDEXES: Sequence[str] = (
    "standard",
    "optional",
    "sagemath",
    "experimental",
    "alph",
)

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

    # PackageClass uses Path predicates internally.  Validate the package tree
    # explicitly first so a permission or I/O error cannot turn a package into
    # an apparently absent one on Python versions that suppress such errors.
    _package_directories()
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
        for suffix, key in INSTALLATION_CATEGORIES:
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


def spkg_index_path(base_dir: Path, name: str) -> Path:
    """Return the file that the ``name`` index of the packages is written to.

    Both :func:`generate_spkg_indexes` and :func:`expected_targets` name their
    files through here, so that neither can come to write or to look for a file
    that the other does not know about.
    """

    if name not in SPKG_INDEXES:
        raise KeyError(f"{name!r} is not one of the package indexes")
    return base_dir / f"index_{name}.rst"


def generate_spkg_indexes(base_dir: Path) -> None:
    """Generate the ``reference/spkg`` index files."""

    documented_names = _documented_package_names()
    base_dir.mkdir(parents=True, exist_ok=True)
    log_install(f"{base_dir}/*.rst")
    written: set[str] = set()

    def write_index(name: str, writer, *args) -> None:
        writer(spkg_index_path(base_dir, name), *args)
        written.add(name)

    write_index(
        "standard",
        write_index_sections,
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

    write_index(
        "optional",
        write_index_sections,
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
        for name in documented_names
        if name.startswith("sagemath_")
    ]
    write_index("sagemath", write_bullet_list, sagemath_packages)

    write_index(
        "experimental",
        write_index_sections,
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

    write_index(
        "alph",
        write_alphabetical_index,
        documented_names,
    )

    unwritten = [name for name in SPKG_INDEXES
                 if name not in written]
    if unwritten:
        raise AssertionError(
            f"SPKG_INDEXES names indexes that nothing writes: {unwritten}")


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
    info = _input_status(path)
    if info is None:
        return []
    if not stat.S_ISREG(info.st_mode):
        raise OSError(f"distro metadata is not a regular file: {path}")
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
                if dep == "FORCE":
                    # Suppress FORCE
                    continue
                elif dep.startswith("$") or dep.startswith("sage"):
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
        info = _input_status(path)
        if info is not None:
            if not stat.S_ISREG(info.st_mode):
                raise OSError(f"package metadata is not a regular file: {path}")
            version_text = path.read_text().strip()
            version_text = [
                line for line in version_text.splitlines() if not line.startswith("#")
            ]
            if version_text:
                lines.extend(
                    [
                        f"{candidate}::",
                        "",
                        *("    " + line for line in version_text),
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
    distros_info = _input_status(distros_dir)
    if distros_info is not None:
        if not stat.S_ISDIR(distros_info.st_mode):
            raise OSError(
                f"package distro metadata is not a directory: {distros_dir}"
            )
        system_files = []
        for path in sorted(distros_dir.iterdir()):
            if path.suffix != ".txt":
                continue
            info = path.stat()
            if not stat.S_ISREG(info.st_mode):
                raise OSError(f"distro metadata is not a regular file: {path}")
            system_files.append(path)
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
    names = _documented_package_names()

    for name in names:
        pkg = Package(name)
        source = Path(pkg.path) / "SPKG.rst"
        if not _input_is_file(source):
            raise OSError(f"package documentation is not a regular file: {source}")
        transformed = transform_spkg_rst(name, source.read_text())
        extra = build_additional_sections(pkg)
        content = f".. _spkg_{name}:\n\n{transformed}{extra}"
        write_text(base_dir / f"{name}.rst", content)


def expected_targets(target_dir: Path) -> list[Path]:
    """Return every file that :func:`main` writes under ``target_dir``.

    Which files those are depends on the package metadata, so whoever wants to
    know whether a documentation tree still has to be generated has to ask
    here: a list kept anywhere else goes stale as soon as a package or a
    supported system is added.
    """

    installation = target_dir / "en" / "installation"
    spkg = target_dir / "en" / "reference" / "spkg"
    targets = [
        installation / f"{system}{suffix}.txt"
        for system in SYSTEM_COMMANDS
        for suffix, _ in INSTALLATION_CATEGORIES
    ]
    targets += [spkg_index_path(spkg, name) for name in SPKG_INDEXES]
    targets += [
        spkg / f"{name}.rst"
        for name in _documented_package_names()
    ]
    return targets


def _hash_file(path: Path) -> str:
    """Return the SHA-256 digest of ``path``."""

    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_regular_file_stably(path: Path) -> tuple[bytes, tuple[int, ...]]:
    """Read one regular path and prove it was not replaced during the read."""

    before = path.lstat()
    if not stat.S_ISREG(before.st_mode):
        raise OSError(f"not a regular file: {path}")
    with path.open("rb") as file:
        opened_before = os.fstat(file.fileno())
        if (_stat_identity(opened_before) != _stat_identity(before)):
            raise OSError(f"file changed before it could be read: {path}")
        content = file.read()
        opened_after = os.fstat(file.fileno())
    after = path.lstat()
    identity = _stat_identity(before)
    if (_stat_identity(opened_after) != identity
            or _stat_identity(after) != identity):
        raise OSError(f"file changed while it was read: {path}")
    return content, identity


def _add_input(content_digest, race_digest, label: str, path: Path) -> None:
    """Add an input's contents and stable identity to two digests."""

    before = path.lstat()
    _add_identity(race_digest, label, "link", before)
    content_digest.update(label.encode("utf-8", errors="surrogateescape"))
    content_digest.update(b"\0")
    target_before = None
    if stat.S_ISLNK(before.st_mode):
        content_digest.update(b"link\0")
        content_digest.update(
            os.readlink(path).encode("utf-8", errors="surrogateescape")
        )
        content_digest.update(b"\0")
        # Package metadata in build/pkgs commonly links to pkgs/*; changing
        # that target is an input change even though the link itself is not.
        try:
            target_before = path.stat()
        except FileNotFoundError:
            target_before = None
        if target_before is None:
            content_digest.update(b"broken")
        elif stat.S_ISREG(target_before.st_mode):
            _add_identity(race_digest, label, "target", target_before)
            content_digest.update(b"file\0")
            content_digest.update(_hash_file(path).encode("ascii"))
        elif stat.S_ISDIR(target_before.st_mode):
            _add_identity(race_digest, label, "target", target_before)
            content_digest.update(b"directory")
        else:
            _add_identity(race_digest, label, "target", target_before)
            content_digest.update(
                f"mode:{stat.S_IFMT(target_before.st_mode)}".encode("ascii")
            )
    elif stat.S_ISREG(before.st_mode):
        content_digest.update(b"file\0")
        content_digest.update(_hash_file(path).encode("ascii"))
    elif stat.S_ISDIR(before.st_mode):
        content_digest.update(b"directory")
    else:
        content_digest.update(
            f"mode:{stat.S_IFMT(before.st_mode)}".encode("ascii")
        )
    content_digest.update(b"\0")

    after = path.lstat()
    if _stat_identity(after) != _stat_identity(before):
        raise RuntimeError(f"input changed while it was read: {path}")
    if stat.S_ISLNK(before.st_mode):
        try:
            target_after = path.stat()
        except FileNotFoundError:
            target_after = None
        if ((target_before is None) != (target_after is None)
                or (target_before is not None
                    and _stat_identity(target_after)
                    != _stat_identity(target_before))):
            raise RuntimeError(f"input link target changed while it was read: {path}")


def _input_state() -> tuple[str, str]:
    """Return content and transient race digests of every generator input.

    Modification times are deliberately absent.  They can be in the future,
    have coarse resolution, or change while the contents do not.  File names,
    file contents, symlink destinations and the contents of symlinked files
    all matter to the first digest instead.  The second digest does include
    identity and timestamps, and is compared only across one generation to
    detect an A-to-B-to-A race without making the persistent manifest stale.
    """

    root = Path(SAGE_ROOT)
    packages = root / "build" / "pkgs"
    helpers = root / "build" / "sage_bootstrap"
    if not _input_is_dir(packages):
        raise FileNotFoundError(
            f"package metadata directory does not exist: {packages}"
        )
    if not _input_is_dir(helpers):
        raise FileNotFoundError(f"bootstrap package does not exist: {helpers}")

    content = hashlib.sha256()
    race = hashlib.sha256()
    _add_input(content, race, "tools/bootstrap-docs.py", Path(__file__))
    _add_input(content, race, "build/sage_bootstrap", helpers)
    for path in _walk_tree(helpers):
        if "__pycache__" in path.parts:
            continue
        relative = path.relative_to(helpers)
        info = path.lstat()
        if stat.S_ISDIR(info.st_mode) or path.suffix == ".py":
            label = f"build/sage_bootstrap/{relative}"
            _add_input(content, race, label, path)
    content.update(b"PYTHON_MINOR\0")
    content.update(os.environ.get("PYTHON_MINOR", "").encode("utf-8"))
    content.update(b"\0")
    race.update(b"PYTHON_MINOR\0")
    race.update(os.environ.get("PYTHON_MINOR", "").encode("utf-8"))
    race.update(b"\0")
    _add_input(content, race, "build/pkgs", packages)
    for package in _package_directories():
        label = f"build/pkgs/{package.name}"
        _add_input(content, race, label, package)
        for name in PACKAGE_INPUT_NAMES:
            path = package / name
            try:
                path.lstat()
            except FileNotFoundError:
                continue
            _add_input(content, race, f"{label}/{name}", path)
        distros = package / "distros"
        distros_info = _input_status(distros)
        if distros_info is not None:
            if not stat.S_ISDIR(distros_info.st_mode):
                raise OSError(
                    f"package distro metadata is not a directory: {distros}"
                )
            _add_input(content, race, f"{label}/distros", distros)
            for path in sorted(distros.iterdir()):
                if path.suffix != ".txt":
                    continue
                _add_input(content, race, f"{label}/distros/{path.name}", path)
    return content.hexdigest(), race.hexdigest()


def input_digest() -> str:
    """Return the stable content digest stored in generated manifests."""

    return _input_state()[0]


def manifest_path(target_dir: Path) -> Path:
    """Return the ignored metadata file describing generated sources."""

    return (
        target_dir / "en" / "installation" / "__pycache__" / MANIFEST_NAME
    )


def _read_manifest(target_dir: Path) -> dict:
    path = manifest_path(target_dir)
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"not a regular file: {path}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or manifest.get("version") != MANIFEST_VERSION:
        raise ValueError(f"unsupported manifest: {path}")
    if not isinstance(manifest.get("inputs"), str):
        raise ValueError(f"invalid input digest in {path}")
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict) or not all(
        isinstance(name, str) and isinstance(value, str)
        for name, value in outputs.items()
    ):
        raise ValueError(f"invalid output digests in {path}")
    return manifest


def write_manifest(target_dir: Path, inputs: str,
                   targets: Sequence[Path] | None = None) -> None:
    """Record exact digests of the inputs and generated regular files."""

    outputs = {}
    if targets is None:
        targets = expected_targets(target_dir)
    for path in targets:
        try:
            content, _ = _read_regular_file_stably(path)
        except OSError as error:
            raise OSError(
                f"generated file is missing, nonregular, or moving: {path}"
            ) from error
        outputs[path.relative_to(target_dir).as_posix()] = hashlib.sha256(
            content
        ).hexdigest()
    manifest = {
        "version": MANIFEST_VERSION,
        "inputs": inputs,
        "outputs": outputs,
    }
    write_text(
        manifest_path(target_dir),
        json.dumps(manifest, indent=2, sort_keys=True),
    )


def complaints_about(target_dir: Path,
                     preserved: list[Path] | None = None) -> list[str]:
    """Return what keeps ``target_dir`` from being a tree Sphinx can read.

    A file that is missing, that an interrupted run left truncated, that the
    package metadata has outgrown, or that a package that no longer exists left
    behind: each of them makes the documentation build fail, and the point of
    asking is to write the sources again before it does.
    """

    complaints = []
    try:
        inputs = input_digest()
    except OSError as error:
        return [f"unreadable inputs: {error}"]
    try:
        targets = expected_targets(target_dir)
    except (OSError, ValueError) as error:
        return [f"unreadable inputs: {error}"]

    manifest = None
    try:
        manifest = _read_manifest(target_dir)
    except (OSError, ValueError) as error:
        complaints.append(
            f"missing or invalid manifest: {manifest_path(target_dir)}: {error}"
        )
    else:
        if manifest["inputs"] != inputs:
            complaints.append("stale: generator inputs changed")

    recorded = manifest["outputs"] if manifest is not None else None
    for path in targets:
        try:
            info = path.lstat()
        except OSError:
            complaints.append(f"missing: {path}")
            continue
        if not stat.S_ISREG(info.st_mode):
            complaints.append(f"not a regular file: {path}")
        elif not info.st_size:
            complaints.append(f"empty: {path}")
        elif recorded is not None:
            relative = path.relative_to(target_dir).as_posix()
            expected = recorded.get(relative)
            if expected is None:
                complaints.append(f"unrecorded: {path}")
            else:
                try:
                    content, _ = _read_regular_file_stably(path)
                except OSError as error:
                    complaints.append(f"unreadable: {path}: {error}")
                else:
                    actual = hashlib.sha256(content).hexdigest()
                    if actual != expected:
                        complaints.append(f"changed or truncated: {path}")

    obsolete, unowned = _generated_extras(target_dir)
    complaints += [f"obsolete: {proof.path}" for proof in obsolete]
    if preserved is not None:
        preserved.extend(unowned)
    return complaints


def _warn_preserved(paths: Sequence[Path]) -> None:
    """Report entries that this generator deliberately leaves untouched."""

    for path in paths:
        print(
            f"warning: unrecognized entry in generated directory "
            f"(preserved): {path}",
            file=sys.stderr,
        )


@dataclass(frozen=True)
class _RemovalProof:
    """The exact obsolete file whose ownership was established."""

    path: Path
    identity: tuple[int, ...]
    digest: str


def _generated_candidates(directory: Path, suffix: str,
                          *, exclude: str | None = None) -> list[Path]:
    """List generated-directory candidates without hiding access errors."""

    info = _input_status(directory)
    if info is None:
        return []
    if not stat.S_ISDIR(info.st_mode):
        raise OSError(f"generated path is not a directory: {directory}")
    return [path for path in sorted(directory.iterdir())
            if path.suffix == suffix and path.name != exclude]


def _generated_extras(
    target_dir: Path, expected_targets_: Sequence[Path] | None = None
) -> tuple[list[_RemovalProof], list[Path]]:
    """Classify extra entries in directories owned by this generator.

    An old manifest proves ownership only while the recorded contents remain
    unchanged.  Before manifests existed, a package page identified itself by
    its first-line ``spkg`` label.  Everything else is preserved and reported:
    a filename pattern alone is not permission to delete a hand-written file.
    """

    if expected_targets_ is None:
        expected_targets_ = expected_targets(target_dir)
    expected = set(expected_targets_)
    installation = target_dir / "en" / "installation"
    spkg = target_dir / "en" / "reference" / "spkg"
    candidates = _generated_candidates(installation, ".txt")
    candidates += _generated_candidates(spkg, ".rst", exclude="index.rst")

    legacy_proof_allowed = False
    try:
        recorded = _read_manifest(target_dir)["outputs"]
    except FileNotFoundError:
        recorded = {}
        legacy_proof_allowed = True
    except (OSError, ValueError):
        # A present but unreadable/corrupt manifest cannot establish that an
        # extra page predates manifests.  Preserve it rather than weakening a
        # failed digest proof to the legacy first-line heuristic.
        recorded = {}

    obsolete = []
    unowned = []
    for path in candidates:
        if path in expected:
            continue
        relative = path.relative_to(target_dir).as_posix()
        recorded_digest = recorded.get(relative)
        try:
            content, identity = _read_regular_file_stably(path)
        except OSError:
            unowned.append(path)
            continue
        digest = hashlib.sha256(content).hexdigest()
        first_line = content.split(b"\n", 1)[0].decode("utf-8", "replace")
        unchanged_recorded_output = recorded_digest == digest
        # The label is only a migration proof for a page predating manifests.
        # Once a manifest names the path, a digest mismatch instead proves
        # that somebody changed it and it must be preserved.
        legacy_page = (legacy_proof_allowed and recorded_digest is None
                       and path.parent == spkg
                       and first_line == f".. _spkg_{path.stem}:")
        if unchanged_recorded_output or legacy_page:
            obsolete.append(_RemovalProof(path, identity, digest))
        else:
            unowned.append(path)
    return obsolete, unowned


def _remove_proven_path(proof: _RemovalProof) -> None:
    """Remove an obsolete output only if the ownership proof is still exact."""

    content, identity = _read_regular_file_stably(proof.path)
    digest = hashlib.sha256(content).hexdigest()
    if identity != proof.identity or digest != proof.digest:
        raise OSError(
            f"refusing to remove generated file changed after inspection: "
            f"{proof.path}"
        )
    # Keep the check adjacent to unlink.  Holding the generator lock prevents
    # another generator from replacing the path between the proof and here.
    remove_path(proof.path)


def main() -> int:
    args = sys.argv[1:]
    check = "--check" in args
    if check:
        args.remove("--check")
    target_dir = Path(args[0])

    if check:
        with _generation_lock(target_dir, exclusive=False):
            preserved: list[Path] = []
            complaints = complaints_about(target_dir, preserved)
            for complaint in complaints:
                print(complaint)
            _warn_preserved(preserved)
            return 1 if complaints else 0

    with _generation_lock(target_dir, exclusive=True):
        inputs, race = _input_state()
        targets = expected_targets(target_dir)
        generate_installation_docs(target_dir / "en" / "installation")
        generate_spkg_indexes(target_dir / "en" / "reference" / "spkg")
        generate_spkg_details(target_dir / "en" / "reference" / "spkg")

        obsolete, unowned = _generated_extras(target_dir, targets)
        _warn_preserved(unowned)

        # Never certify output produced from a moving input tree.  Comparing
        # both content and metadata catches A-to-B-to-A changes as well as a
        # persistent edit.  This check precedes the only destructive action.
        after, after_race = _input_state()
        if (after, after_race) != (inputs, race):
            raise RuntimeError(
                "package metadata changed while documentation sources were generated"
            )
        for proof in obsolete:
            log_install(f"removing {proof.path}")
            _remove_proven_path(proof)

        write_manifest(target_dir, after, targets)
        final, final_race = _input_state()
        if (final, final_race) != (inputs, race):
            raise RuntimeError(
                "package metadata changed while documentation sources were generated"
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
