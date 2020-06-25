# -*- coding: utf-8 -*-
r"""
Testing for features of the environment at runtime

A computation can require a certain package to be installed in the runtime
environment. Abstractly such a package describes a :class`Feature` which can
be tested for at runtime. It can be of various kinds, most prominently an
:class:`Executable` in the PATH or an additional package for some installed
system such as a :class:`GapPackage`.

AUTHORS:

- Julian Rüth (2016-04-07): Initial version

- Jeroen Demeyer (2018-02-12): Refactoring and clean up

EXAMPLES:

Some generic features are available for common cases. For example, to
test for the existence of a binary, one can use an :class:`Executable`
feature::

    sage: from sage.features import Executable
    sage: Executable(name="sh", executable="sh").is_present()
    FeatureTestResult('sh', True)

Here we test whether the grape GAP package is available::

    sage: from sage.features.gap import GapPackage
    sage: GapPackage("grape", spkg="gap_packages").is_present()  # optional: gap_packages
    FeatureTestResult('GAP package grape', True)

Note that a :class:`FeatureTestResult` acts like a bool in most contexts::

    sage: if Executable(name="sh", executable="sh").is_present(): "present."
    'present.'

When one wants to raise an error if the feature is not available, one
can use the ``require`` method::

    sage: Executable(name="sh", executable="sh").require()

    sage: Executable(name="random", executable="randomOochoz6x", spkg="random", url="http://rand.om").require()
    Traceback (most recent call last):
    ...
    FeatureNotPresentError: random is not available.
    Executable 'randomOochoz6x' not found on PATH.
    To install random you can try to run 'sage -i random'.
    Further installation instructions might be available at http://rand.om.

As can be seen above, features try to produce helpful error messages.
"""

import os
from distutils.errors import CCompilerError
from distutils.spawn import find_executable

from sage.env import SAGE_SHARE


class TrivialClasscallMetaClass(type):
    """
    A trivial version of :class:`ClasscallMetaclass` without Cython dependencies.
    """
    def __call__(cls, *args, **kwds):
        r"""
        This method implements ``cls(<some arguments>)``.
        """
        if hasattr(cls, '__classcall__'):
            return cls.__classcall__(cls, *args, **kwds)
        else:
            return type.__call__(cls, *args, **kwds)

_trivial_unique_representation_cache = dict()

class TrivialUniqueRepresentation(metaclass=TrivialClasscallMetaClass):
    """
    A trivial version of :class:`UniqueRepresentation` without Cython dependencies.
    """

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        Construct a new object of this class or reuse an existing one.
        """
        key = (cls, tuple(args), frozenset(options.items()))
        cached = _trivial_unique_representation_cache.get(key, None)
        if cached is None:
            cached = _trivial_unique_representation_cache[key] = type.__call__(cls, *args, **options)
        return cached

class Feature(TrivialUniqueRepresentation):
    r"""
    A feature of the runtime environment

    Overwrite :meth:`_is_present` to add feature checks.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")  # indirect doctest
        Feature('GAP package grape')

    For efficiency, features are unique::

        sage: GapPackage("grape") is GapPackage("grape")
        True
    """
    def __init__(self, name, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.features import Feature
            sage: from sage.features.gap import GapPackage
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), Feature)  # indirect doctest
            True
        """
        self.name = name
        self.spkg = spkg
        self.url = url
        self._cache_is_present = None

    def is_present(self):
        r"""
        Return whether the feature is present.

        OUTPUT:

        A :class:`FeatureTestResult` which can be used as a boolean and
        contains additional information about the feature test.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present()  # optional: gap_packages
            FeatureTestResult('GAP package grape', True)
            sage: GapPackage("NOT_A_PACKAGE", spkg="gap_packages").is_present()
            FeatureTestResult('GAP package NOT_A_PACKAGE', False)

        The result is cached::

            sage: from sage.features import Feature
            sage: class TestFeature(Feature):
            ....:     def _is_present(self):
            ....:         print("checking presence")
            ....:         return True
            sage: TestFeature("test").is_present()
            checking presence
            FeatureTestResult('test', True)
            sage: TestFeature("test").is_present()
            FeatureTestResult('test', True)
            sage: TestFeature("other").is_present()
            checking presence
            FeatureTestResult('other', True)
            sage: TestFeature("other").is_present()
            FeatureTestResult('other', True)
        """
        # We do not use @cached_method here because we wish to use
        # Feature early in the build system of sagelib.
        if self._cache_is_present is None:
            res = self._is_present()
            if not isinstance(res, FeatureTestResult):
                res = FeatureTestResult(self, res)
            self._cache_is_present = res
        return self._cache_is_present

    def _is_present(self):
        """
        Override this in a derived class to implement the feature check.

        This should return either an instance of
        :class:`FeatureTestResult` or a boolean.
        """
        raise NotImplementedError("_is_present not implemented for feature {!r}".format(self.name))

    def require(self):
        r"""
        Raise a :class:`FeatureNotPresentError` if the feature is not present.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("ve1EeThu").require()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package ve1EeThu is not available.
            `TestPackageAvailability("ve1EeThu")` evaluated to `fail` in GAP.
        """
        presence = self.is_present()
        if not presence:
            raise FeatureNotPresentError(self, presence.reason, presence.resolution)

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape")  # indirect doctest
            Feature('GAP package grape')
        """
        return 'Feature({name!r})'.format(name=self.name)

    def resolution(self):
        r"""
        Return a suggestion on how to make :meth:`is_present` pass if it did not
        pass.

        EXAMPLES::

            sage: from sage.features import Executable
            sage: Executable(name="CSDP", spkg="csdp", executable="theta", url="http://github.org/dimpase/csdp").resolution()
            "To install CSDP you can try to run 'sage -i csdp'.\nFurther installation instructions might be available at http://github.org/dimpase/csdp."
        """
        lines = []
        if self.spkg:
            lines.append("To install {feature} you can try to run 'sage -i {spkg}'.".format(feature=self.name, spkg=self.spkg))
        if self.url:
            lines.append("Further installation instructions might be available at {url}.".format(url=self.url))
        return "\n".join(lines) or None


class FeatureNotPresentError(RuntimeError):
    r"""
    A missing feature error.

    EXAMPLES::

        sage: from sage.features import Feature, FeatureTestResult
        sage: class Missing(Feature):
        ....:     def _is_present(self):
        ....:         return False

        sage: Missing(name="missing").require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: missing is not available.
    """
    def __init__(self, feature, reason=None, resolution=None):
        self.feature = feature
        self.reason = reason
        self.resolution = resolution

    def __str__(self):
        r"""
        Return the error message.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("gapZuHoh8Uu").require()  # indirect doctest
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package gapZuHoh8Uu is not available.
            `TestPackageAvailability("gapZuHoh8Uu")` evaluated to `fail` in GAP.
        """
        lines = ["{feature} is not available.".format(feature=self.feature.name)]
        if self.reason:
            lines.append(self.reason)
        if self.resolution:
            lines.append(self.resolution)
        return "\n".join(lines)


class FeatureTestResult(object):
    r"""
    The result of a :meth:`Feature.is_present` call.

    Behaves like a boolean with some extra data which may explain why a feature
    is not present and how this may be resolved.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: presence = GapPackage("NOT_A_PACKAGE").is_present(); presence  # indirect doctest
        FeatureTestResult('GAP package NOT_A_PACKAGE', False)
        sage: bool(presence)
        False

    Explanatory messages might be available as ``reason`` and
    ``resolution``::

        sage: presence.reason
        '`TestPackageAvailability("NOT_A_PACKAGE")` evaluated to `fail` in GAP.'
        sage: print(presence.resolution)
        None

    If a feature is not present, ``resolution`` defaults to
    ``feature.resolution()`` if this is defined. If you do not want to use this
    default you need explicitly set ``resolution`` to a string::

        sage: from sage.features import FeatureTestResult
        sage: package = GapPackage("NOT_A_PACKAGE", spkg="no_package")
        sage: FeatureTestResult(package, True).resolution
        "To install GAP package NOT_A_PACKAGE you can try to run 'sage -i no_package'."
        sage: FeatureTestResult(package, False).resolution
        "To install GAP package NOT_A_PACKAGE you can try to run 'sage -i no_package'."
        sage: FeatureTestResult(package, False, resolution="rtm").resolution
        'rtm'
    """
    def __init__(self, feature, is_present, reason=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.features import Executable, FeatureTestResult
            sage: isinstance(Executable(name="sh", executable="sh").is_present(), FeatureTestResult)
            True
        """
        self.feature = feature
        self.is_present = is_present
        self.reason = reason
        self.resolution = resolution or feature.resolution()

    def __bool__(self):
        r"""
        Whether the tested :class:`Feature` is present.

        TESTS::

            sage: from sage.features import Feature, FeatureTestResult
            sage: bool(FeatureTestResult(Feature("SomePresentFeature"), True))  # indirect doctest
            True
            sage: bool(FeatureTestResult(Feature("SomeMissingFeature"), False))
            False
        """
        return bool(self.is_present)

    __nonzero__ = __bool__

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.features import Feature, FeatureTestResult
            sage: FeatureTestResult(Feature("SomePresentFeature"), True)  # indirect doctest
            FeatureTestResult('SomePresentFeature', True)
        """
        return "FeatureTestResult({feature!r}, {is_present!r})".format(feature=self.feature.name, is_present=self.is_present)


class Executable(Feature):
    r"""
    A feature describing an executable in the PATH.

    .. NOTE::

        Overwrite :meth:`is_functional` if you also want to check whether
        the executable shows proper behaviour.

        Calls to :meth:`is_present` are cached. You might want to cache the
        :class:`Executable` object to prevent unnecessary calls to the
        executable.

    EXAMPLES::

        sage: from sage.features import Executable
        sage: Executable(name="sh", executable="sh").is_present()
        FeatureTestResult('sh', True)
    """
    def __init__(self, name, executable, **kwds):
        r"""
        TESTS::

            sage: from sage.features import Executable
            sage: isinstance(Executable(name="sh", executable="sh"), Executable)
            True
        """
        Feature.__init__(self, name, **kwds)
        self.executable = executable

    def _is_present(self):
        r"""
        Test whether the executable is on the current PATH and functional.

        .. SEEALSO:: :meth:`is_functional`

        EXAMPLES::

            sage: from sage.features import Executable
            sage: Executable(name="sh", executable="sh").is_present()
            FeatureTestResult('sh', True)
        """
        if find_executable(self.executable) is None:
            return FeatureTestResult(self, False, "Executable {executable!r} not found on PATH.".format(executable=self.executable))
        return self.is_functional()

    def is_functional(self):
        r"""
        Return whether an executable in the path is functional.

        EXAMPLES:

        Returns ``True`` unless explicitly overwritten::

            sage: from sage.features import Executable
            sage: Executable(name="sh", executable="sh").is_functional()
            FeatureTestResult('sh', True)
        """
        return FeatureTestResult(self, True)


class StaticFile(Feature):
    r"""
    A :class:`Feature` which describes the presence of a certain file such as a
    database.

    EXAMPLES::

        sage: from sage.features import StaticFile
        sage: StaticFile(name="no_such_file", filename="KaT1aihu", search_path=("/",), spkg="some_spkg", url="http://rand.om").require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: no_such_file is not available.
        'KaT1aihu' not found in any of ['/']
        To install no_such_file you can try to run 'sage -i some_spkg'.
        Further installation instructions might be available at http://rand.om.
    """
    def __init__(self, name, filename, search_path=None, **kwds):
        r"""
        TESTS::

            sage: from sage.features import StaticFile
            sage: StaticFile(name="null", filename="null", search_path=("/dev",))
            Feature('null')
        """
        Feature.__init__(self, name, **kwds)
        self.filename = filename
        if search_path is None:
            self.search_path = [SAGE_SHARE]
        else:
            self.search_path = list(search_path)

    def _is_present(self):
        r"""
        Whether the static file is present.

           sage: from sage.features import StaticFile
           sage: StaticFile(name="no_such_file", filename="KaT1aihu", spkg="some_spkg", url="http://rand.om").is_present()
           FeatureTestResult('no_such_file', False)
        """
        try:
            abspath = self.absolute_path()
            return FeatureTestResult(self, True, reason="Found at `{abspath}`.".format(abspath=abspath))
        except FeatureNotPresentError as e:
            return FeatureTestResult(self, False, reason=e.reason, resolution=e.resolution)

    def absolute_path(self):
        r"""
        The absolute path of the file.

        EXAMPLES::

            sage: from sage.features.databases import DatabaseCremona
            sage: DatabaseCremona().absolute_path()  # optional: database_cremona_ellcurve
            '.../local/share/cremona/cremona.db'

        A ``FeatureNotPresentError`` is raised if the file can not be found::

            sage: from sage.features import StaticFile
            sage: StaticFile(name="no_such_file", filename="KaT1aihu", search_path=(), spkg="some_spkg", url="http://rand.om").absolute_path()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: no_such_file is not available.
            'KaT1aihu' not found in any of []
            To install no_such_file you can try to run 'sage -i some_spkg'.
            Further installation instructions might be available at http://rand.om.
        """
        for directory in self.search_path:
            path = os.path.join(directory, self.filename)
            if os.path.isfile(path):
                return os.path.abspath(path)
        raise FeatureNotPresentError(self,
            reason="{filename!r} not found in any of {search_path}".format(filename=self.filename, search_path=self.search_path),
            resolution=self.resolution())


class CythonFeature(Feature):
    r"""
    A :class:`Feature` which describes the ability to compile and import
    a particular piece of Cython code.

    To test the presence of ``name``, the cython compiler is run on
    ``test_code`` and the resulting module is imported.

    EXAMPLES::

        sage: from sage.features import CythonFeature
        sage: fabs_test_code = '''
        ....: cdef extern from "<math.h>":
        ....:     double fabs(double x)
        ....:
        ....: assert fabs(-1) == 1
        ....: '''
        sage: fabs = CythonFeature("fabs", test_code=fabs_test_code, spkg="gcc", url="https://gnu.org")
        sage: fabs.is_present()
        FeatureTestResult('fabs', True)

    Test various failures::

        sage: broken_code = '''this is not a valid Cython program!'''
        sage: broken = CythonFeature("broken", test_code=broken_code)
        sage: broken.is_present()
        FeatureTestResult('broken', False)

    ::

        sage: broken_code = '''cdef extern from "no_such_header_file": pass'''
        sage: broken = CythonFeature("broken", test_code=broken_code)
        sage: broken.is_present()
        FeatureTestResult('broken', False)

    ::

        sage: broken_code = '''import no_such_python_module'''
        sage: broken = CythonFeature("broken", test_code=broken_code)
        sage: broken.is_present()
        FeatureTestResult('broken', False)

    ::

        sage: broken_code = '''raise AssertionError("sorry!")'''
        sage: broken = CythonFeature("broken", test_code=broken_code)
        sage: broken.is_present()
        FeatureTestResult('broken', False)
    """
    def __init__(self, name, test_code, **kwds):
        r"""
        TESTS::

            sage: from sage.features import CythonFeature
            sage: from sage.features.fes import LibFESLibrary
            sage: isinstance(LibFESLibrary(), CythonFeature)  # indirect doctest
            True
        """
        Feature.__init__(self, name, **kwds)
        self.test_code = test_code

    def _is_present(self):
        r"""
        Run test code to determine whether the shared library is present.

        EXAMPLES::

            sage: from sage.features import CythonFeature
            sage: empty = CythonFeature("empty", test_code="")
            sage: empty.is_present()
            FeatureTestResult('empty', True)
        """
        from sage.misc.temporary_file import tmp_filename
        with open(tmp_filename(ext=".pyx"), 'w') as pyx:
            pyx.write(self.test_code)
        from sage.misc.cython import cython_import
        try:
            cython_import(pyx.name, verbose=-1)
        except CCompilerError:
            return FeatureTestResult(self, False, reason="Failed to compile test code.")
        except ImportError:
            return FeatureTestResult(self, False, reason="Failed to import test code.")
        except Exception:
            return FeatureTestResult(self, False, reason="Failed to run test code.")
        return FeatureTestResult(self, True, reason="Test code compiled and imported.")


class PythonModule(Feature):
    r"""
    A :class:`Feature` which describes whether a python module can be imported.

    EXAMPLES:

    Not all builds of python include the ``ssl`` module, so you could check
    whether it is available::

        sage: from sage.features import PythonModule
        sage: PythonModule("ssl").require()  # not tested - output depends on the python build
    """
    def __init__(self, name, **kwds):
        r"""
        TESTS::

            sage: from sage.features import PythonModule
            sage: from sage.features.fes import LibFES
            sage: isinstance(LibFES(), PythonModule)  # indirect doctest
            True
        """
        Feature.__init__(self, name, **kwds)

    def _is_present(self):
        r"""
        Return whether the module can be imported. This is determined by
        actually importing it.

        EXAMPLES::

            sage: from sage.features import PythonModule
            sage: PythonModule("sys").is_present()
            FeatureTestResult('sys', True)
            sage: PythonModule("_no_such_module_").is_present()
            FeatureTestResult('_no_such_module_', False)
        """
        import importlib
        try:
            importlib.import_module(self.name)
        except ImportError:
            return FeatureTestResult(self, False, reason="Failed to import `{name}`.".format(name=self.name))
        return FeatureTestResult(self, True, reason="Successfully imported `{name}`.".format(name=self.name))
