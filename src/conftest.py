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
from pathlib import Path
from typing import Any, Iterable, Optional

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
            except ImportError as exception:
                if self.config.getvalue("doctest_ignore_import_errors"):
                    pytest.skip("unable to import module %r" % self.path)
                else:
                    if isinstance(exception, ModuleNotFoundError):
                        # Ignore some missing features/modules for now
                        # TODO: Remove this once all optional things are using Features
                        if exception.name in (
                            "valgrind",
                            "rpy2",
                            "sage.libs.coxeter3.coxeter",
                        ):
                            pytest.skip(
                                f"unable to import module { self.path } due to missing feature { exception.name }"
                            )
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
                f"unable to import module { self.path } due to missing feature { exception.feature.name }"
            )
        except ModuleNotFoundError as exception:
            # TODO: Remove this once all optional things are using Features
            pytest.skip(
                f"unable to import module { self.path } due to missing feature { exception.name }"
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
    elif file_path.suffix == ".py":
        if parent.config.option.doctest:
            if file_path.name == "__main__.py" or file_path.name == "setup.py":
                # We don't allow tests to be defined in __main__.py/setup.py files (because their import will fail).
                return IgnoreCollector.from_parent(parent)
            if (
                (
                    file_path.name == "postprocess.py"
                    and file_path.parent.name == "nbconvert"
                )
                or (
                    file_path.name == "giacpy-mkkeywords.py"
                    and file_path.parent.name == "autogen"
                )
                or (
                    file_path.name == "flint_autogen.py"
                    and file_path.parent.name == "autogen"
                )
            ):
                # This is an executable file.
                return IgnoreCollector.from_parent(parent)

            if file_path.name == "conftest_inputtest.py":
                # This is an input file for testing the doctest machinery (and contains broken doctests).
                return IgnoreCollector.from_parent(parent)

            if (
                (
                    file_path.name == "finite_dimensional_lie_algebras_with_basis.py"
                    and file_path.parent.name == "categories"
                )
                or (
                    file_path.name == "__init__.py"
                    and file_path.parent.name == "crypto"
                )
                or (file_path.name == "__init__.py" and file_path.parent.name == "mq")
            ):
                # TODO: Fix these (import fails with "RuntimeError: dictionary changed size during iteration")
                return IgnoreCollector.from_parent(parent)
            return SageDoctestModule.from_parent(parent, path=file_path)


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
