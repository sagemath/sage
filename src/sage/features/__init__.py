# sage_setup: distribution = sagemath-environment
r"""
Testing for features of the environment at runtime

A computation can require a certain package to be installed in the runtime
environment. Abstractly such a package describes a :class:`Feature` which can
be tested for at runtime. It can be of various kinds, most prominently an
:class:`Executable` in the ``PATH``, a :class:`PythonModule`, or an additional
package for some installed
system such as a :class:`~sage.features.gap.GapPackage`.

AUTHORS:

- Julian Rüth (2016-04-07): Initial version

- Jeroen Demeyer (2018-02-12): Refactoring and clean up

EXAMPLES:

Some generic features are available for common cases. For example, to
test for the existence of a binary, one can use an :class:`Executable`
feature::

    sage: from sage.features import Executable
    sage: Executable(name='sh', executable='sh').is_present()
    FeatureTestResult('sh', True)

Here we test whether the grape GAP package is available::

    sage: from sage.features.gap import GapPackage
    sage: GapPackage("grape", spkg='gap_packages').is_present()  # optional - gap_package_grape
    FeatureTestResult('gap_package_grape', True)

Note that a :class:`FeatureTestResult` acts like a bool in most contexts::

    sage: if Executable(name='sh', executable='sh').is_present(): "present."
    'present.'

When one wants to raise an error if the feature is not available, one
can use the ``require`` method::

    sage: Executable(name='sh', executable='sh').require()

    sage: Executable(name='random', executable='randomOochoz6x', spkg='random', url='http://rand.om').require() # optional - sage_spkg
    Traceback (most recent call last):
    ...
    FeatureNotPresentError: random is not available.
    Executable 'randomOochoz6x' not found on PATH.
    ...try to run...sage -i random...
    Further installation instructions might be available at http://rand.om.

As can be seen above, features try to produce helpful error messages.
"""

# *****************************************************************************
#       Copyright (C) 2016      Julian Rüth
#                     2018      Jeroen Demeyer
#                     2018      Timo Kaufmann
#                     2019-2022 Matthias Koeppe
#                     2021      Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import annotations

import os
import shutil
from pathlib import Path

from sage.env import SAGE_SHARE, SAGE_LOCAL, SAGE_VENV


class TrivialClasscallMetaClass(type):
    """
    A trivial version of :class:`sage.misc.classcall_metaclass.ClasscallMetaclass` without Cython dependencies.
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
    r"""
    A trivial version of :class:`UniqueRepresentation` without Cython dependencies.
    """

    @staticmethod
    def __classcall__(cls, *args, **options):
        r"""
        Construct a new object of this class or reuse an existing one.
        """
        key = (cls, tuple(args), frozenset(options.items()))
        cached = _trivial_unique_representation_cache.get(key, None)
        if cached is None:
            cached = _trivial_unique_representation_cache[key] = type.__call__(cls, *args, **options)
        return cached


class Feature(TrivialUniqueRepresentation):
    r"""
    A feature of the runtime environment.

    INPUT:

    - ``name`` -- string; name of the feature. This should be suitable as an optional tag
      for the Sage doctester, i.e., lowercase alphanumeric with underscores (``_``) allowed;
      features that correspond to Python modules/packages may use periods (``.``)

    - ``spkg`` -- string; name of the SPKG providing the feature

    - ``description`` -- string (optional); plain English description of the feature

    - ``url`` -- a URL for the upstream package providing the feature

    - ``type`` -- string; one of ``'standard'``, ``'optional'`` (default), ``'experimental'``

    Overwrite :meth:`_is_present` to add feature checks.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: GapPackage("grape", spkg='gap_packages')  # indirect doctest
        Feature('gap_package_grape')

    For efficiency, features are unique::

        sage: GapPackage("grape") is GapPackage("grape")
        True
    """
    def __init__(self, name, spkg=None, url=None, description=None, type='optional'):
        r"""
        TESTS::

            sage: from sage.features import Feature
            sage: from sage.features.gap import GapPackage
            sage: isinstance(GapPackage("grape", spkg='gap_packages'), Feature)  # indirect doctest
            True
        """
        self.name = name
        self.spkg = spkg
        self.url = url
        self.description = description

        self._cache_is_present = None
        self._cache_resolution = None
        self._hidden = False
        self._type = type

        try:
            from sage.misc.package import spkg_type
        except ImportError:  # may have been surgically removed in a downstream distribution
            pass
        else:
            if spkg and (t := spkg_type(spkg)) not in (type, None):
                from warnings import warn
                warn(f'Feature {name} is declared {type}, '
                     f'but it is provided by {spkg}, '
                     f'which is declared {t} in SAGE_ROOT/build/pkgs',
                     stacklevel=3)

    def is_present(self):
        r"""
        Return whether the feature is present.

        OUTPUT:

        A :class:`FeatureTestResult` which can be used as a boolean and
        contains additional information about the feature test.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape", spkg='gap_packages').is_present()  # optional - gap_package_grape
            FeatureTestResult('gap_package_grape', True)
            sage: GapPackage("NOT_A_PACKAGE", spkg='gap_packages').is_present()
            FeatureTestResult('gap_package_NOT_A_PACKAGE', False)

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

        if self._hidden:
            return FeatureTestResult(self, False, reason="Feature `{name}` is hidden.".format(name=self.name))

        return self._cache_is_present

    def _is_present(self):
        r"""
        Override this in a derived class to implement the feature check.

        This should return either an instance of
        :class:`FeatureTestResult` or a boolean.
        """
        raise NotImplementedError("_is_present not implemented for feature {!r}".format(self.name))

    def require(self):
        r"""
        Raise a :exc:`FeatureNotPresentError` if the feature is not present.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("ve1EeThu").require()                                      # needs sage.libs.gap
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: gap_package_ve1EeThu is not available.
            `LoadPackage("ve1EeThu")` evaluated to `fail` in GAP.
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
            Feature('gap_package_grape')
        """
        description = f'{self.name!r}: {self.description}' if self.description else f'{self.name!r}'
        return f'Feature({description})'

    def _spkg_type(self):
        r"""
        Return the type of this feature.

        For features provided by an SPKG in the Sage distribution,
        this should match the SPKG type, or a warning will be issued.

        EXAMPLES::

            sage: from sage.features.databases import DatabaseCremona
            sage: DatabaseCremona()._spkg_type()
            'optional'

        OUTPUT:

        The type as a string in ``('base', 'standard', 'optional', 'experimental')``.
        """
        return self._type

    def resolution(self):
        r"""
        Return a suggestion on how to make :meth:`is_present` pass if it did not
        pass.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.features import Executable
            sage: Executable(name='CSDP', spkg='csdp', executable='theta', url='https://github.com/dimpase/csdp').resolution()  # optional - sage_spkg
            '...To install CSDP...you can try to run...sage -i csdp...Further installation instructions might be available at https://github.com/dimpase/csdp.'
        """
        if self._hidden:
            return "Use method `unhide` to make it available again."
        if self._cache_resolution is not None:
            return self._cache_resolution
        lines = []
        if self.spkg:
            for ps in package_systems():
                lines.append(ps.spkg_installation_hint(self.spkg, feature=self.name))
        if self.url:
            lines.append("Further installation instructions might be available at {url}.".format(url=self.url))
        self._cache_resolution = "\n".join(lines)
        return self._cache_resolution

    def joined_features(self):
        r"""
        Return a list of features that ``self`` is the join of.

        OUTPUT:

        A (possibly empty) list of instances of :class:`Feature`.

        EXAMPLES::

            sage: from sage.features.graphviz import Graphviz
            sage: Graphviz().joined_features()
            [Feature('dot'), Feature('neato'), Feature('twopi')]
            sage: from sage.features.sagemath import sage__rings__function_field
            sage: sage__rings__function_field().joined_features()
            [Feature('sage.rings.function_field.function_field_polymod'),
            Feature('sage.libs.singular'),
            Feature('sage.libs.singular.singular'),
            Feature('sage.interfaces.singular')]
            sage: from sage.features.interfaces import Mathematica
            sage: Mathematica().joined_features()
            []
        """
        from sage.features.join_feature import JoinFeature
        res = []
        if isinstance(self, JoinFeature):
            for f in self._features:
                res += [f] + f.joined_features()
        return res

    def is_standard(self):
        r"""
        Return whether this feature corresponds to a standard SPKG.

        EXAMPLES::

            sage: from sage.features.databases import DatabaseCremona
            sage: DatabaseCremona().is_standard()
            False
        """
        if self.name.startswith('sage.'):
            return True
        return self._spkg_type() == 'standard'

    def is_optional(self):
        r"""
        Return whether this feature corresponds to an optional SPKG.

        EXAMPLES::

            sage: from sage.features.databases import DatabaseCremona
            sage: DatabaseCremona().is_optional()
            True
        """
        return self._spkg_type() == 'optional'

    def hide(self):
        r"""
        Hide this feature. For example this is used when the doctest option
        ``--hide`` is set. Setting an installed feature as hidden pretends
        that it is not available. To revert this use :meth:`unhide`.

        EXAMPLES:

        Benzene is an optional SPKG. The following test fails if it is hidden or
        not installed. Thus, in the second invocation the optional tag is needed::

            sage: from sage.features.graph_generators import Benzene
            sage: Benzene().hide()
            sage: len(list(graphs.fusenes(2)))                                          # needs sage.graphs
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: benzene is not available.
            Feature `benzene` is hidden.
            Use method `unhide` to make it available again.

            sage: Benzene().unhide()            # optional - benzene, needs sage.graphs
            sage: len(list(graphs.fusenes(2)))  # optional - benzene, needs sage.graphs
            1
        """
        self._hidden = True

    def unhide(self):
        r"""
        Revert what :meth:`hide` did.

        EXAMPLES:

            sage: from sage.features.sagemath import sage__plot
            sage: sage__plot().hide()
            sage: sage__plot().is_present()
            FeatureTestResult('sage.plot', False)
            sage: sage__plot().unhide()                                                 # needs sage.plot
            sage: sage__plot().is_present()                                             # needs sage.plot
            FeatureTestResult('sage.plot', True)
        """
        self._hidden = False

    def is_hidden(self):
        r"""
        Return whether ``self`` is present but currently hidden.

        EXAMPLES:

            sage: from sage.features.sagemath import sage__plot
            sage: sage__plot().hide()
            sage: sage__plot().is_hidden()                                              # needs sage.plot
            True
            sage: sage__plot().unhide()
            sage: sage__plot().is_hidden()
            False
        """
        if self._hidden and self._is_present():
            return True
        return False

class FeatureNotPresentError(RuntimeError):
    r"""
    A missing feature error.

    EXAMPLES::

        sage: from sage.features import Feature, FeatureTestResult
        sage: class Missing(Feature):
        ....:     def _is_present(self):
        ....:         return False

        sage: Missing(name='missing').require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: missing is not available.
    """
    def __init__(self, feature, reason=None, resolution=None):
        self.feature = feature
        self.reason = reason
        self._resolution = resolution

    @property
    def resolution(self):
        if self._resolution:
            return self._resolution
        return self.feature.resolution()

    def __str__(self):
        r"""
        Return the error message.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("gapZuHoh8Uu").require()  # indirect doctest               # needs sage.libs.gap
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: gap_package_gapZuHoh8Uu is not available.
            `LoadPackage("gapZuHoh8Uu")` evaluated to `fail` in GAP.
        """
        lines = ["{feature} is not available.".format(feature=self.feature.name)]
        if self.reason:
            lines.append(self.reason)
        resolution = self.resolution
        if resolution:
            lines.append(str(resolution))
        return "\n".join(lines)


class FeatureTestResult:
    r"""
    The result of a :meth:`Feature.is_present` call.

    Behaves like a boolean with some extra data which may explain why a feature
    is not present and how this may be resolved.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: presence = GapPackage("NOT_A_PACKAGE").is_present(); presence  # indirect doctest
        FeatureTestResult('gap_package_NOT_A_PACKAGE', False)
        sage: bool(presence)
        False

    Explanatory messages might be available as ``reason`` and
    ``resolution``::

        sage: presence.reason                                                           # needs sage.libs.gap
        '`LoadPackage("NOT_A_PACKAGE")` evaluated to `fail` in GAP.'
        sage: bool(presence.resolution)
        False

    If a feature is not present, ``resolution`` defaults to
    ``feature.resolution()`` if this is defined. If you do not want to use this
    default you need explicitly set ``resolution`` to a string::

        sage: from sage.features import FeatureTestResult
        sage: package = GapPackage("NOT_A_PACKAGE", spkg='no_package')
        sage: str(FeatureTestResult(package, True).resolution)  # optional - sage_spkg
        '...To install gap_package_NOT_A_PACKAGE...you can try to run...sage -i no_package...'
        sage: str(FeatureTestResult(package, False).resolution) # optional - sage_spkg
        '...To install gap_package_NOT_A_PACKAGE...you can try to run...sage -i no_package...'
        sage: FeatureTestResult(package, False, resolution='rtm').resolution
        'rtm'
    """
    def __init__(self, feature, is_present, reason=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.features import Executable, FeatureTestResult
            sage: isinstance(Executable(name='sh', executable='sh').is_present(), FeatureTestResult)
            True
        """
        self.feature = feature
        self.is_present = is_present
        self.reason = reason
        self._resolution = resolution

    @property
    def resolution(self):
        if self._resolution:
            return self._resolution
        return self.feature.resolution()

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

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.features import Feature, FeatureTestResult
            sage: FeatureTestResult(Feature("SomePresentFeature"), True)  # indirect doctest
            FeatureTestResult('SomePresentFeature', True)
        """
        return "FeatureTestResult({feature!r}, {is_present!r})".format(feature=self.feature.name, is_present=self.is_present)


_cache_package_systems = None


def package_systems():
    """
    Return a list of :class:`~sage.features.pkg_systems.PackageSystem` objects
    representing the available package systems.

    The list is ordered by decreasing preference.

    EXAMPLES::

        sage: from sage.features import package_systems
        sage: package_systems()    # random
        [Feature('homebrew'), Feature('sage_spkg'), Feature('pip')]
    """
    # The current implementation never returns more than one system.
    from subprocess import run, CalledProcessError, PIPE
    global _cache_package_systems
    if _cache_package_systems is None:
        from .pkg_systems import PackageSystem, SagePackageSystem, PipPackageSystem
        _cache_package_systems = []
        # Try to use scripts from SAGE_ROOT (or an installation of sage_bootstrap)
        # to obtain system package advice.
        try:
            proc = run('sage-guess-package-system', shell=True, capture_output=True, text=True, check=True)
            system_name = proc.stdout.strip()
            if system_name != 'unknown':
                _cache_package_systems = [PackageSystem(system_name)]
        except CalledProcessError:
            pass
        more_package_systems = [SagePackageSystem(), PipPackageSystem()]
        _cache_package_systems += [ps for ps in more_package_systems if ps.is_present()]

    return _cache_package_systems


class FileFeature(Feature):
    r"""
    Base class for features that describe a file or directory in the file system.

    A subclass should implement a method :meth:`absolute_filename`.

    EXAMPLES:

    Two direct concrete subclasses of :class:`FileFeature` are defined::

        sage: from sage.features import StaticFile, Executable, FileFeature
        sage: issubclass(StaticFile, FileFeature)
        True
        sage: issubclass(Executable, FileFeature)
        True

    To work with the file described by the feature, use the method :meth:`absolute_filename`.
    A :exc:`FeatureNotPresentError` is raised if the file cannot be found::

        sage: Executable(name='does-not-exist', executable='does-not-exist-xxxxyxyyxyy').absolute_filename()
        Traceback (most recent call last):
        ...
        sage.features.FeatureNotPresentError: does-not-exist is not available.
        Executable 'does-not-exist-xxxxyxyyxyy' not found on PATH.

    A :class:`FileFeature` also provides the :meth:`is_present` method to test for
    the presence of the file at run time. This is inherited from the base class
    :class:`Feature`::

        sage: Executable(name='sh', executable='sh').is_present()
        FeatureTestResult('sh', True)
    """
    def _is_present(self):
        r"""
        Whether the file is present.

        EXAMPLES::

           sage: from sage.features import StaticFile
           sage: StaticFile(name='no_such_file', filename='KaT1aihu', spkg='some_spkg', url='http://rand.om').is_present()
           FeatureTestResult('no_such_file', False)
        """
        try:
            abspath = self.absolute_filename()
            return FeatureTestResult(self, True, reason="Found at `{abspath}`.".format(abspath=abspath))
        except FeatureNotPresentError as e:
            return FeatureTestResult(self, False, reason=e.reason, resolution=e.resolution)

    def absolute_filename(self) -> str:
        r"""
        The absolute path of the file as a string.

        Concrete subclasses must override this abstract method.

        TESTS::

            sage: from sage.features import FileFeature
            sage: FileFeature(name='abstract_file').absolute_filename()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # We do not use sage.misc.abstract_method here because that is provided by
        # the distribution sagemath-objects, which is not an install-requires of
        # the distribution sagemath-environment.
        raise NotImplementedError


class Executable(FileFeature):
    r"""
    A feature describing an executable in the ``PATH``.

    In an installation of Sage with ``SAGE_LOCAL`` different from ``SAGE_VENV``, the
    executable is searched first in ``SAGE_VENV/bin``, then in ``SAGE_LOCAL/bin``,
    then in ``PATH``.

    .. NOTE::

        Overwrite :meth:`is_functional` if you also want to check whether
        the executable shows proper behaviour.

        Calls to :meth:`is_present` are cached. You might want to cache the
        :class:`Executable` object to prevent unnecessary calls to the
        executable.

    EXAMPLES::

        sage: from sage.features import Executable
        sage: Executable(name='sh', executable='sh').is_present()
        FeatureTestResult('sh', True)
        sage: Executable(name='does-not-exist', executable='does-not-exist-xxxxyxyyxyy').is_present()
        FeatureTestResult('does-not-exist', False)
    """
    def __init__(self, name, executable, **kwds):
        r"""
        TESTS::

            sage: from sage.features import Executable
            sage: isinstance(Executable(name='sh', executable='sh'), Executable)
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
            sage: Executable(name='sh', executable='sh').is_present()
            FeatureTestResult('sh', True)
        """
        result = FileFeature._is_present(self)
        if not result:
            return result
        return self.is_functional()

    def is_functional(self):
        r"""
        Return whether an executable in the path is functional.

        This method is used internally and can be overridden in subclasses
        in order to implement a feature test. It should not be called directly.
        Use :meth:`Feature.is_present` instead.

        EXAMPLES:

        The function returns ``True`` unless explicitly overwritten::

            sage: from sage.features import Executable
            sage: Executable(name='sh', executable='sh').is_functional()
            FeatureTestResult('sh', True)
        """
        return FeatureTestResult(self, True)

    def absolute_filename(self) -> str:
        r"""
        The absolute path of the executable as a string.

        EXAMPLES::

            sage: from sage.features import Executable
            sage: Executable(name='sh', executable='sh').absolute_filename()
            '/...bin/sh'

        A :exc:`FeatureNotPresentError` is raised if the file cannot be found::

            sage: Executable(name='does-not-exist', executable='does-not-exist-xxxxyxyyxyy').absolute_filename()
            Traceback (most recent call last):
            ...
            sage.features.FeatureNotPresentError: does-not-exist is not available.
            Executable 'does-not-exist-xxxxyxyyxyy' not found on PATH.
        """
        if SAGE_LOCAL:
            if Path(SAGE_VENV).resolve() != Path(SAGE_LOCAL).resolve():
                # As sage.env currently gives SAGE_LOCAL a fallback value from SAGE_VENV,
                # SAGE_LOCAL is never unset.  So we only use it if it differs from SAGE_VENV.
                search_path = ':'.join([os.path.join(SAGE_VENV, 'bin'),
                                        os.path.join(SAGE_LOCAL, 'bin')])
                path = shutil.which(self.executable, path=search_path)
                if path is not None:
                    return path
        # Now look up in the regular PATH.
        path = shutil.which(self.executable)
        if path is not None:
            return path
        raise FeatureNotPresentError(self,
                                     reason="Executable {executable!r} not found on PATH.".format(executable=self.executable),
                                     resolution=self.resolution())


class StaticFile(FileFeature):
    r"""
    A :class:`Feature` which describes the presence of a certain file such as a
    database.

    EXAMPLES::

        sage: from sage.features import StaticFile
        sage: StaticFile(name='no_such_file', filename='KaT1aihu',              # optional - sage_spkg
        ....:            search_path='/', spkg='some_spkg',
        ....:            url='http://rand.om').require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: no_such_file is not available.
        'KaT1aihu' not found in any of ['/']...
        To install no_such_file...you can try to run...sage -i some_spkg...
        Further installation instructions might be available at http://rand.om.
    """
    def __init__(self, name, filename, *, search_path=None, type='optional', **kwds):
        r"""
        TESTS::

            sage: from sage.features import StaticFile
            sage: StaticFile(name='null', filename='null', search_path='/dev')
            Feature('null')
            sage: sh = StaticFile(name='shell', filename='sh',
            ....:                 search_path=("/dev", "/bin", "/usr"))
            sage: sh
            Feature('shell')
            sage: sh.absolute_filename()
            '/bin/sh'
        """
        Feature.__init__(self, name, type=type, **kwds)
        self.filename = filename
        if search_path is None:
            self.search_path = [SAGE_SHARE]
        elif isinstance(search_path, str):
            self.search_path = [search_path]
        else:
            self.search_path = list(search_path)

    def absolute_filename(self) -> str:
        r"""
        The absolute path of the file as a string.

        EXAMPLES::

            sage: from sage.features import StaticFile
            sage: from sage.misc.temporary_file import tmp_dir
            sage: dir_with_file = tmp_dir()
            sage: file_path = os.path.join(dir_with_file, "file.txt")
            sage: open(file_path, 'a').close() # make sure the file exists
            sage: search_path = ( '/foo/bar', dir_with_file ) # file is somewhere in the search path
            sage: feature = StaticFile(name='file', filename='file.txt', search_path=search_path)
            sage: feature.absolute_filename() == file_path
            True

        A :exc:`FeatureNotPresentError` is raised if the file cannot be found::

            sage: from sage.features import StaticFile
            sage: StaticFile(name='no_such_file', filename='KaT1aihu',\
                             search_path=(), spkg='some_spkg',\
                             url='http://rand.om').absolute_filename()  # optional - sage_spkg
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: no_such_file is not available.
            'KaT1aihu' not found in any of []...
            To install no_such_file...you can try to run...sage -i some_spkg...
            Further installation instructions might be available at http://rand.om.
        """
        for directory in self.search_path:
            path = os.path.join(directory, self.filename)
            if os.path.isfile(path) or os.path.isdir(path):
                return os.path.abspath(path)
        reason = "{filename!r} not found in any of {search_path}".format(filename=self.filename, search_path=self.search_path)
        raise FeatureNotPresentError(self, reason=reason, resolution=self.resolution())


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
        sage: fabs = CythonFeature("fabs", test_code=fabs_test_code,                    # needs sage.misc.cython
        ....:                      spkg='gcc', url='https://gnu.org',
        ....:                      type='standard')
        sage: fabs.is_present()                                                         # needs sage.misc.cython
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
            sage: from sage.features.bliss import BlissLibrary
            sage: isinstance(BlissLibrary(), CythonFeature)  # indirect doctest
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
            sage: empty.is_present()                                                    # needs sage.misc.cython
            FeatureTestResult('empty', True)
        """
        from sage.misc.temporary_file import tmp_filename
        try:
            # Available since https://setuptools.pypa.io/en/latest/history.html#v59-0-0
            from setuptools.errors import CCompilerError
        except ImportError:
            try:
                from distutils.errors import CCompilerError
            except ImportError:
                CCompilerError = ()
        with open(tmp_filename(ext='.pyx'), 'w') as pyx:
            pyx.write(self.test_code)
        try:
            from sage.misc.cython import cython_import
        except ImportError:
            return FeatureTestResult(self, False, reason="sage.misc.cython is not available")
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
            sage: from sage.features.databases import DatabaseKnotInfo
            sage: isinstance(DatabaseKnotInfo(), PythonModule)  # indirect doctest
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
        except ImportError as exception:
            return FeatureTestResult(self, False, reason=f"Failed to import `{self.name}`: {exception}")
        return FeatureTestResult(self, True, reason=f"Successfully imported `{self.name}`.")
