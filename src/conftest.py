# pyright: strict

"""Configuration and fixtures for pytest.

This file configures pytest and provides some global fixtures.
See https://docs.pytest.org/en/latest/index.html for more details.
"""

from __future__ import annotations

import inspect
from pathlib import Path
from typing import Any, Iterable

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
from _pytest.pathlib import import_path, ImportMode

# Import sage.all is necessary to:
# - avoid cyclic import errors, see Trac #33580
# - inject it into globals namespace for doctests
import sage.all
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
            )
        else:
            try:
                module = import_path(
                    self.path,
                    mode=ImportMode.importlib,
                    root=self.config.rootpath,
                )
            except ImportError:
                if self.config.getvalue("doctest_ignore_import_errors"):
                    pytest.skip("unable to import module %r" % self.path)
                else:
                    raise
        # Uses internal doctest module parsing mechanism.
        finder = MockAwareDocTestFinder()
        optionflags = get_optionflags(self)
        runner = _get_runner(
            verbose=False,
            optionflags=optionflags,
            checker=SageOutputChecker(),
            continue_on_failure=_get_continue_on_failure(self.config),
        )

        for test in finder.find(module, module.__name__):
            if test.examples:  # skip empty doctests
                yield DoctestItem.from_parent(
                    self, name=test.name, runner=runner, dtest=test
                )


def pytest_collect_file(
    file_path: Path, parent: pytest.File
) -> pytest.Collector | None:
    """
    This hook is called when collecting test files, and can be used to
    modify the file or test selection logic by returning a list of
    ``pytest.Item`` objects which the ``pytest`` command will directly
    add to the list of test items.

    See `pytest documentation <https://docs.pytest.org/en/latest/reference/reference.html#std-hook-pytest_collect_file>`_.
    """
    if file_path.suffix == ".pyx":
        # We don't allow pytests to be defined in Cython files.
        # Normally, Cython files are filtered out already by pytest and we only
        # hit this here if someone explicitly runs `pytest some_file.pyx`.
        return pytest.skip("Skipping Cython file")
    elif file_path.suffix == ".py":
        if parent.config.option.doctestmodules:
            return SageDoctestModule.from_parent(parent, path=file_path)


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace: dict[str, Any]):
    """
    Add global imports for doctests.

    See `pytest documentation <https://docs.pytest.org/en/stable/doctest.html#doctest-namespace-fixture>`.
    """
    # Inject sage.all into each doctest
    dict_all = sage.all.__dict__

    # Remove '__package__' item from the globals since it is not
    # always in the globals in an actual Sage session.
    dict_all.pop("__package__", None)

    sage_namespace = dict(dict_all)
    sage_namespace["__name__"] = "__main__"

    doctest_namespace.update(**sage_namespace)
