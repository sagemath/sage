# pyright: strict
"""Configuration and fixtures for pytest.

This file configures pytest and provides some global fixtures.
See https://docs.pytest.org/en/latest/index.html for more details.
"""

from __future__ import annotations

import doctest
import inspect
import sys
import warnings
from typing import TYPE_CHECKING, Any, Optional

import pytest
from _pytest.doctest import (
    DoctestItem,
    DoctestModule,
    _get_continue_on_failure,
    _get_runner,
    _is_mocked,
    _patch_unwrap_mock_aware,
    get_optionflags,
)
from _pytest.pathlib import ImportMode, import_path

from sage.doctest.forker import (
    init_sage,
    showwarning_with_traceback,
)
from sage.doctest.parsing import SageDocTestParser, SageOutputChecker

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path


# Stash key holding the per-session Sage random seed (resolved in
# ``pytest_configure`` and applied by the ``set_random_seed`` fixture).
_random_seed_key = pytest.StashKey[int]()


def _resolve_lazy_members(obj: object) -> None:
    """
    Eagerly resolve lazy members that cache themselves into ``obj``'s namespace.

    The stdlib doctest finder iterates over ``obj.__dict__`` for both modules
    and classes. Touching a lazy member during that walk resolves it and writes
    the result back into the namespace, raising ``RuntimeError: dictionary
    changed size during iteration``. Resolving them up front avoids that, and
    does no extra work overall: the finder would resolve the same members while
    walking ``obj``, and the resolved values are cached globally.

    Two kinds of lazy member cause this:

    - :class:`~sage.misc.lazy_import.LazyImport` objects stored directly in a
      module (or class) namespace; and

    - :class:`~sage.misc.lazy_attribute.lazy_class_attribute` descriptors, which
      are looked up along the MRO but cache their value into the *subclass* that
      they are accessed on (e.g. ``_axiom`` on a category with axiom).
    """
    from sage.misc.lazy_attribute import lazy_class_attribute
    from sage.misc.lazy_import import LazyImport

    if isinstance(obj, type):
        names = {
            name
            for klass in obj.__mro__
            for name, value in list(vars(klass).items())
            if isinstance(value, (LazyImport, lazy_class_attribute))
        }
        for name in names:
            try:
                getattr(obj, name)
            except Exception:
                # Leave unresolvable members in place; the finder's usual
                # missing-feature/module handling applies when they are reached.
                pass
    elif inspect.ismodule(obj):
        for value in list(vars(obj).values()):
            if isinstance(value, LazyImport):
                try:
                    value._get_object()
                except Exception:
                    pass


def is_subpath(path: Path, parent: Path) -> bool:
    # Check if the path is in a subdirectory, or a subsubdirectory, ... of the parent
    path = path.resolve()
    parent = parent.resolve()
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


class SageDoctestModule(DoctestModule):
    """
    This is essentially a copy of `DoctestModule` from
    https://github.com/pytest-dev/pytest/blob/main/src/_pytest/doctest.py.
    The only change is that we use `SageDocTestParser` to extract the doctests
    and `SageOutputChecker` to verify the output.
    """

    def collect(self) -> Iterable[DoctestItem]:
        import doctest

        class MockAwareDocTestFinder(doctest.DocTestFinder):
            """A hackish doctest finder that overrides stdlib internals to fix a stdlib bug.
            https://github.com/pytest-dev/pytest/issues/3456
            https://bugs.python.org/issue25532
            """

            def __init__(self) -> None:
                super().__init__(parser=SageDocTestParser(set(["sage"])))

            def _find_lineno(self, obj, source_lines):
                """Doctest code does not take into account `@property`, this
                is a hackish way to fix it. https://bugs.python.org/issue17446
                Wrapped Doctests will need to be unwrapped so the correct
                line number is returned. This will be reported upstream. #8796
                """
                if isinstance(obj, property):
                    obj = getattr(obj, "fget", obj)

                if hasattr(obj, "__wrapped__"):
                    # Get the main obj in case of it being wrapped
                    obj = inspect.unwrap(obj)

                # Type ignored because this is a private function.
                return super()._find_lineno(  # type:ignore[misc]
                    obj,
                    source_lines,
                )

            def _find(
                self, tests, obj, name, module, source_lines, globs, seen
            ) -> None:
                if _is_mocked(obj):
                    return
                # Resolve lazy members of obj before super()._find iterates
                # obj.__dict__, to avoid "dictionary changed size during
                # iteration" when the walk triggers a lazy import or
                # lazy_class_attribute.
                _resolve_lazy_members(obj)
                with _patch_unwrap_mock_aware():
                    # Type ignored because this is a private function.
                    super()._find(  # type:ignore[misc]
                        tests, obj, name, module, source_lines, globs, seen
                    )

        if self.path.name == "conftest.py":
            module = self.config.pluginmanager._importconftest(
                self.path,
                self.config.getoption("importmode"),
                rootpath=self.config.rootpath,
                consider_namespace_packages=True,
            )
        else:
            try:
                module = import_path(
                    self.path,
                    mode=ImportMode.importlib,
                    root=self.config.rootpath,
                    consider_namespace_packages=True,
                )
            except ImportError:
                if self.config.getvalue("doctest_ignore_import_errors"):
                    pytest.skip("unable to import module %r" % self.path)
                else:
                    raise
        # Uses internal doctest module parsing mechanism.
        finder = MockAwareDocTestFinder()
        optionflags = get_optionflags(self.config)
        from sage.features import FeatureNotPresentError

        runner = _get_runner(
            verbose=False,
            optionflags=optionflags,
            checker=SageOutputChecker(),
            continue_on_failure=_get_continue_on_failure(self.config),
        )
        try:
            for test in finder.find(module, module.__name__):
                if test.examples:  # skip empty doctests
                    yield DoctestItem.from_parent(
                        self, name=test.name, runner=runner, dtest=test
                    )
        except FeatureNotPresentError as exception:
            pytest.skip(
                f"unable to import module {self.path} due to missing feature {exception.feature.name}"
            )
        except ModuleNotFoundError as exception:
            # TODO: Remove this once all optional things are using Features
            pytest.skip(
                f"unable to import module {self.path} due to missing module {exception.name}"
            )


class IgnoreCollector(pytest.Collector):
    """
    Ignore a file.
    """

    def __init__(self, parent: pytest.Collector) -> None:
        super().__init__("ignore", parent)

    def collect(self) -> Iterable[pytest.Item | pytest.Collector]:
        return []


def pytest_collect_file(
    file_path: Path, parent: pytest.Collector
) -> pytest.Collector | None:
    """
    This hook is called when collecting test files, and can be used to
    modify the file or test selection logic by returning a list of
    ``pytest.Item`` objects which the ``pytest`` command will directly
    add to the list of test items.

    See `pytest documentation <https://docs.pytest.org/en/latest/reference/reference.html#std-hook-pytest_collect_file>`_.
    """
    if (
        file_path.parent.name == "combinat"
        or file_path.parent.parent.name == "combinat"
    ):
        # Crashes CI for some reason
        return IgnoreCollector.from_parent(parent)
    if file_path.suffix == ".pyx":
        # We don't allow pytests to be defined in Cython files.
        # Normally, Cython files are filtered out already by pytest and we only
        # hit this here if someone explicitly runs `pytest some_file.pyx`.
        return IgnoreCollector.from_parent(parent)
    if file_path.suffix == ".py":
        if parent.config.option.doctest:
            if file_path.name == "__main__.py" or file_path.name == "setup.py":
                # We don't allow tests to be defined in __main__.py/setup.py files (because their import will fail).
                return IgnoreCollector.from_parent(parent)
            if (
                file_path.name == "postprocess.py"
                and file_path.parent.name == "nbconvert"
            ):
                # This is an executable file.
                return IgnoreCollector.from_parent(parent)

            if (
                file_path.name in ("forker.py", "reporting.py")
            ) and file_path.parent.name == "doctest":
                # Fails with many errors due to different testing framework
                return IgnoreCollector.from_parent(parent)

            if (
                file_path.name == "finitely_presented.py"
                and file_path.parent.name == "groups"
            ):
                # Passes under `sage -t` but not here: pytest extracts doctests
                # from __doc__ (where Python collapses backslash-continuations
                # and processes escapes in non-raw docstrings) while sage -t
                # reads the raw source, plus some order-dependent GAP state.
                return IgnoreCollector.from_parent(parent)

            return SageDoctestModule.from_parent(parent, path=file_path)


def pytest_ignore_collect(
    collection_path: Path, config: pytest.Config
) -> bool | None:
    """
    This hook is called when collecting test files, and can be used to
    prevent considering this path for collection by returning ``True``.

    See `pytest documentation <https://docs.pytest.org/en/latest/reference/reference.html#pytest.hookspec.pytest_ignore_collect>`_.
    """
    root = config.rootpath
    if (
        is_subpath(collection_path, root / "src" / "sage_docbuild")
        or is_subpath(collection_path, root / "src" / "sage_setup")
        or collection_path == root / "src" / "build-docs.py"
    ):
        # Fails to import with Meson
        return True
    if collection_path.name == "all.py":
        # all.py do not contain tests and may fail when imported twice / in the wrong order
        return True


def pytest_addoption(parser):
    # Add a command line option to run doctests
    # (we don't use the built-in --doctest-modules option because then doctests are collected twice)
    group = parser.getgroup("collect")
    group.addoption(
        "--doctest",
        action="store_true",
        default=False,
        help="Run doctests in all .py modules",
        dest="doctest",
    )
    # Mirror `sage -t`: long-running tests are skipped unless explicitly
    # requested. Tests are tagged with the ``long`` / ``longlong`` markers
    # (declared in ``pyproject.toml``).
    group.addoption(
        "--long",
        action="store_true",
        default=False,
        help="Also run tests marked as long (skipped by default)",
        dest="run_long",
    )
    group.addoption(
        "--longlong",
        action="store_true",
        default=False,
        help="Also run tests marked as long or longlong (skipped by default)",
        dest="run_longlong",
    )
    group.addoption(
        "--random-seed",
        type=int,
        default=None,
        metavar="SEED",
        help=(
            "Seed for Sage's random number generator, set before each test "
            "(default: a fresh random seed, reported in the test header). "
            "Can also be set via the SAGE_PYTEST_RANDOM_SEED environment variable."
        ),
        dest="random_seed",
    )


def pytest_configure(config: pytest.Config) -> None:
    """
    Resolve the Sage random seed once per session.

    The seed is taken from ``--random-seed``, falling back to the
    ``SAGE_PYTEST_RANDOM_SEED`` environment variable, and finally to a fresh
    random seed. It is stashed on the config and applied before each test by
    the :func:`set_random_seed` fixture, mirroring ``sage -t``.
    """
    import os

    from sage.misc import randstate

    seed = config.getoption("random_seed")
    if seed is None:
        env_seed = os.environ.get("SAGE_PYTEST_RANDOM_SEED")
        seed = int(env_seed) if env_seed else None
    if seed is None:
        # Let Sage pick a fresh seed and record it so the run is reproducible.
        randstate.set_random_seed()
        seed = randstate.initial_seed()
    config.stash[_random_seed_key] = seed


def pytest_report_header(config: pytest.Config) -> str:
    """Report the random seed so a failing run can be reproduced."""
    seed = config.stash[_random_seed_key]
    return f"Sage random seed: {seed} (re-run with --random-seed={seed})"


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]):
    """
    Skip tests marked ``long`` / ``longlong`` unless the corresponding
    command-line option is given.

    This mirrors the behaviour of ``sage -t``, where long-running tests are
    only executed when ``--long`` is passed. ``--longlong`` implies ``--long``.

    See `pytest documentation <https://docs.pytest.org/en/stable/reference/reference.html#std-hook-pytest_collection_modifyitems>`_.
    """
    run_longlong = config.getoption("run_longlong")
    run_long = config.getoption("run_long") or run_longlong

    skip_long = pytest.mark.skip(reason="need --long option to run")
    skip_longlong = pytest.mark.skip(reason="need --longlong option to run")

    for item in items:
        if not run_longlong and "longlong" in item.keywords:
            item.add_marker(skip_longlong)
        elif not run_long and "long" in item.keywords:
            item.add_marker(skip_long)

        _skip_if_features_missing(item)


def _skip_if_features_missing(item: pytest.Item) -> None:
    """
    Honour the ``optional`` marker: skip ``item`` unless every named Sage
    feature is available, using the same detection as the doctest framework.

    Feature names are those accepted by the doctest ``# optional - ...`` and
    ``# needs ...`` tags, e.g. ``sage.symbolic``, ``sage.plot``, ``latex``,
    ``pynormaliz``. Usage::

        @pytest.mark.optional("sage.plot", "latex")
        def test_something():
            ...
    """
    features = [
        name for marker in item.iter_markers(name="optional") for name in marker.args
    ]
    if not features:
        return

    from sage.doctest.external import available_software

    missing = [name for name in features if name not in available_software]
    if missing:
        item.add_marker(
            pytest.mark.skip(reason="missing Sage feature(s): " + ", ".join(missing))
        )


# Monkey patch exception printing to replace the full qualified name of the exception by its short name
# TODO: Remove this hack once migration to pytest is complete
import traceback

old_format_exception_only = traceback.format_exception_only


def format_exception_only(etype: type, value: BaseException) -> list[str]:
    formatted_exception = old_format_exception_only(etype, value)
    exception_name = etype.__name__
    if etype.__module__:
        exception_full_name = etype.__module__ + "." + etype.__qualname__
    else:
        exception_full_name = etype.__qualname__

    for i, line in enumerate(formatted_exception):
        if line.startswith(exception_full_name):
            formatted_exception[i] = line.replace(
                exception_full_name, exception_name, 1
            )
    return formatted_exception


# Initialize Sage-specific doctest stuff
init_sage()

# Monkey patch doctest to use our custom printer etc
old_run = doctest.DocTestRunner.run


def doctest_run(
    self: doctest.DocTestRunner,
    test: doctest.DocTest,
    compileflags: Optional[int] = None,
    out: Any = None,
    clear_globs: bool = True,
) -> doctest.TestResults:
    from sage.repl.rich_output import get_display_manager
    from sage.repl.user_globals import set_globals

    traceback.format_exception_only = format_exception_only

    # Display warnings in doctests
    warnings.showwarning = showwarning_with_traceback
    setattr(sys, "__displayhook__", get_display_manager().displayhook)

    # Ensure that injecting globals works as expected in doctests
    set_globals(test.globs)
    return old_run(self, test, compileflags, out, clear_globs)


doctest.DocTestRunner.run = doctest_run


@pytest.fixture(autouse=True)
def set_random_seed(request: pytest.FixtureRequest):
    """
    Seed Sage's random number generator before each test.

    This mirrors ``sage -t``: the same seed (resolved once per session in
    :func:`pytest_configure`) is applied before every test, so a run is
    reproducible with ``--random-seed=<seed>`` (the seed is printed in the
    test header).
    """
    from sage.misc.randstate import set_random_seed as sage_set_random_seed

    sage_set_random_seed(request.config.stash[_random_seed_key])


@pytest.fixture(autouse=True, scope="session")
def add_imports(doctest_namespace: dict[str, Any]):
    """
    Add global imports for doctests.

    See `pytest documentation <https://docs.pytest.org/en/stable/doctest.html#doctest-namespace-fixture>`.
    """
    # Inject sage.all into each doctest
    import sage.repl.ipython_kernel.all_jupyter

    dict_all = sage.repl.ipython_kernel.all_jupyter.__dict__

    # Remove '__package__' item from the globals since it is not
    # always in the globals in an actual Sage session.
    dict_all.pop("__package__", None)

    sage_namespace = dict(dict_all)
    sage_namespace["__name__"] = "__main__"

    doctest_namespace.update(**sage_namespace)


@pytest.fixture
def tmpfile():
    r"""
    Temporary file fixture that can be reopened/closed and still
    clean itself up afterwards.

    Similar to the built-in ``tmpdir`` fixture, but safer for now:

    * https://github.com/pytest-dev/pytest/issues/13669

    """
    from os import unlink
    from tempfile import NamedTemporaryFile
    t = NamedTemporaryFile(delete=False)
    yield t
    unlink(t.name)


@pytest.fixture
def assert_close():
    r"""
    Assert that two (possibly symbolic or exact) numbers are numerically close.

    This is the Sage-aware counterpart of :func:`pytest.approx` for the cases
    it does not handle directly, in particular symbolic and exact values:
    both arguments are evaluated numerically (via :func:`complex`) before
    being compared. Real, complex, rational, and symbolic inputs are all
    accepted.

    The tolerance follows :func:`math.isclose`: the values are close when
    ``abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)``. This mirrors
    the ``# rel tol`` / ``# abs tol`` doctest flags.

    EXAMPLES (used as a pytest fixture)::

        def test_sqrt(assert_close):
            from sage.all import sqrt
            assert_close(sqrt(2), 1.4142135623730951)
            assert_close(sqrt(2), 1.41421, rel_tol=1e-5)
    """
    import cmath

    def _assert_close(actual, expected, *, rel_tol=1e-9, abs_tol=0.0):
        a = complex(actual)
        e = complex(expected)
        if not cmath.isclose(a, e, rel_tol=rel_tol, abs_tol=abs_tol):
            raise AssertionError(
                f"{actual!r} is not close to {expected!r} "
                f"(|diff| = {abs(a - e):.3e}, rel_tol={rel_tol}, abs_tol={abs_tol})"
            )

    return _assert_close
