# -*- coding: utf-8 -*-
"""
Sage Packages
"""

# ****************************************************************************
#       Copyright (C) 2015-2016 Volker Braun <vbraun.name@gmail.com>
#                     2018      Jeroen Demeyer
#                     2020-2024 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import logging
import os
import re

from sage_bootstrap.env import SAGE_ROOT

log = logging.getLogger()


class Package(object):

    def __new__(cls, package_name):
        if package_name.startswith("pypi/") or package_name.startswith("generic/"):
            package_name = "pkg:" + package_name
        if package_name.startswith("pkg:"):
            package_name = package_name.replace('_', '-').lower()
            if package_name.startswith("pkg:generic/"):  # fast path
                try:
                    pkg = cls(package_name[len("pkg:generic/"):].replace('-', '_'))
                    if pkg.purl == package_name:
                        return pkg  # assume unique
                except Exception:
                    pass
            elif package_name.startswith("pkg:pypi/"):  # fast path
                try:
                    pkg = cls(package_name[len("pkg:pypi/"):].replace('-', '_'))
                    if pkg.purl == package_name:
                        return pkg  # assume unique
                except Exception:
                    pass
            for pkg in cls.all():
                if pkg.purl == package_name:
                    return pkg  # assume unique
            raise ValueError('no package for PURL {0}'.format(package_name))
        self = object.__new__(cls)
        self.__init__(package_name)
        return self

    def __init__(self, package_name):
        """
        Sage Package

        A package is defined by a subdirectory of
        ``SAGE_ROOT/build/pkgs/``. The name of the package is the name
        of the subdirectory; The metadata of the package is contained
        in various files in the package directory. This class provides
        an abstraction to the metadata, you should never need to
        access the package directory directly.

        INPUT:

        -- ``package_name`` -- string. Name of the package. The Sage
           convention is that all package names are lower case.
        """
        if any(package_name.startswith(prefix)
               for prefix in ["pkg:", "pypi/", "generic"]):
            # Already initialized
            return
        if package_name != package_name.lower():
            raise ValueError('package names should be lowercase, got {0}'.format(package_name))
        if '-' in package_name:
            raise ValueError('package names use underscores, not dashes, got {0}'.format(package_name))

        self.__name = package_name
        self.__tarball = None
        self._init_checksum()
        self._init_version()
        self._init_type()
        self._init_version_requirements()
        self._init_requirements()
        self._init_dependencies()
        self._init_trees()

    def __repr__(self):
        return 'Package {0}'.format(self.name)

    @property
    def name(self):
        """
        Return the package name

         A package is defined by a subdirectory of
        ``SAGE_ROOT/build/pkgs/``. The name of the package is the name
        of the subdirectory.

        OUTPUT:

        String.
        """
        return self.__name

    @property
    def sha1(self):
        """
        Return the SHA1 checksum

        OUTPUT:

        String.
        """
        return self.__sha1

    @property
    def sha256(self):
        """
        Return the SHA256 checksum

        OUTPUT:

        String.
        """
        return self.__sha256

    @property
    def tarball(self):
        """
        Return the (primary) tarball

        If there are multiple tarballs (currently unsupported), this
        property returns the one that is unpacked automatically.

        OUTPUT:

        Instance of :class:`sage_bootstrap.tarball.Tarball`
        """
        if self.__tarball is None:
            from sage_bootstrap.tarball import Tarball
            self.__tarball = Tarball(self.tarball_filename, package=self)
        return self.__tarball

    def _substitute_variables_once(self, pattern):
        """
        Substitute (at most) one occurrence of variables in ``pattern`` by the values.

        These variables are ``VERSION``, ``VERSION_MAJOR``, ``VERSION_MINOR``,
        ``VERSION_MICRO``, either appearing like this or in the form ``${VERSION_MAJOR}``
        etc.

        Return a tuple:
        - the string with the substitution done or the original string
        - whether a substitution was done
        """
        for var in ('VERSION_MAJOR', 'VERSION_MINOR', 'VERSION_MICRO', 'VERSION'):
            # As VERSION is a substring of the other three, it needs to be tested last.
            dollar_brace_var = '${' + var + '}'
            if dollar_brace_var in pattern:
                value = getattr(self, var.lower())
                return pattern.replace(dollar_brace_var, value, 1), True
            elif var in pattern:
                value = getattr(self, var.lower())
                return pattern.replace(var, value, 1), True
        return pattern, False

    def _substitute_variables(self, pattern):
        """
        Substitute all occurrences of ``VERSION`` in ``pattern`` by the actual version.

        Likewise for ``VERSION_MAJOR``, ``VERSION_MINOR``, ``VERSION_MICRO``,
        either appearing like this or in the form ``${VERSION}``, ``${VERSION_MAJOR}``,
        etc.
        """
        not_done = True
        while not_done:
            pattern, not_done = self._substitute_variables_once(pattern)
        return pattern

    @property
    def tarball_pattern(self):
        """
        Return the (primary) tarball file pattern

        If there are multiple tarballs (currently unsupported), this
        property returns the one that is unpacked automatically.

        OUTPUT:

        String. The full-qualified tarball filename, but with
        ``VERSION`` instead of the actual tarball filename.
        """
        return self.__tarball_pattern

    @property
    def tarball_filename(self):
        """
        Return the (primary) tarball filename

        If there are multiple tarballs (currently unsupported), this
        property returns the one that is unpacked automatically.

        OUTPUT:

        String. The full-qualified tarball filename.
        """
        pattern = self.tarball_pattern
        if pattern:
            return self._substitute_variables(pattern)
        else:
            return None

    @property
    def tarball_upstream_url_pattern(self):
        """
        Return the tarball upstream URL pattern

        OUTPUT:

        String. The tarball upstream URL, but with the placeholder
        ``VERSION``.
        """
        return self.__tarball_upstream_url_pattern

    @property
    def tarball_upstream_url(self):
        """
        Return the tarball upstream URL or ``None`` if none is recorded

        OUTPUT:

        String. The URL.
        """
        pattern = self.tarball_upstream_url_pattern
        if pattern:
            return self._substitute_variables(pattern)
        else:
            return None

    @property
    def tarballs_info(self):
        """
        Return information about all tarballs for this package.
        
        This supports packages with multiple platform-specific wheels.
        
        OUTPUT:
        
        List of dictionaries, each containing:
        - 'tarball': tarball filename pattern
        - 'sha256': SHA256 checksum
        - 'sha1': SHA1 checksum (optional)
        - 'upstream_url': upstream URL pattern
        """
        return self.__tarballs_info
    
    def find_tarball_for_platform(self, python_version=None):
        """
        Find the appropriate tarball for the current platform.
        
        For packages with multiple platform-specific wheels, this selects
        the one matching the current platform and Python version using
        the packaging.tags module to ensure compatibility.
        
        Properly handles wheel ABI tags:
        - cp313-cp313 (CPython 3.13 specific)
        - cp313-cp313t (CPython 3.13 free-threaded/nogil)
        - cp313-abi3 (stable ABI, forward compatible)
        - py3-none-any (universal pure Python wheel)
        - pp39-pypy39_pp73 (PyPy wheels)
        
        INPUT:
        
        - ``python_version`` -- Python version string (e.g., '3.11'), or None to auto-detect
        
        OUTPUT:
        
        Dictionary with tarball info, or None if no suitable tarball found.
        The dictionary contains the same fields as tarballs_info entries.
        """
        import subprocess
        
        if not self.__tarballs_info:
            return None
        
        # If only one tarball, return it
        if len(self.__tarballs_info) == 1:
            return self.__tarballs_info[0]
        
        # Get compatible tags from Sage's Python using packaging.tags
        from sage_bootstrap.env import SAGE_ROOT
        sage_script = os.path.join(SAGE_ROOT, 'sage')
        if not os.path.exists(sage_script):
            raise RuntimeError('Sage script not found at: {0}'.format(sage_script))
        
        try:
            # Get all compatible tags from Sage's Python
            result = subprocess.run(
                [sage_script, '-python', '-c', 
                 'import packaging.tags; import json; tags = [str(t) for t in packaging.tags.sys_tags()]; print(json.dumps(tags))'],
                capture_output=True,
                text=True,
                timeout=10,
                cwd=SAGE_ROOT
            )
            if result.returncode != 0:
                raise RuntimeError('Failed to get compatible tags from sage -python: {0}'.format(result.stderr))
            
            import json
            compatible_tags = json.loads(result.stdout.strip())
            
        except subprocess.TimeoutExpired:
            raise RuntimeError('Timeout while querying compatible tags via ./sage -python')
        except Exception as e:
            raise RuntimeError('Error querying compatible tags via ./sage -python: {0}'.format(str(e)))
        
        # Convert tags list to a set for fast lookup with priority
        # Lower index = higher priority
        tag_priority = {tag: idx for idx, tag in enumerate(compatible_tags)}
        
        # Import packaging utilities for parsing wheel filenames
        try:
            import packaging.utils
        except ImportError:
            raise RuntimeError('packaging module not available')
        
        # Find the best matching tarball
        best_match = None
        best_priority = float('inf')
        
        for tarball_info in self.__tarballs_info:
            tarball = tarball_info['tarball']
            
            # Handle non-wheel tarballs (source distributions)
            if not tarball.endswith('.whl'):
                # Source distributions have lowest priority
                if best_priority > len(compatible_tags) + 1000:
                    best_match = tarball_info
                    best_priority = len(compatible_tags) + 1000
                continue
            
            # Parse wheel filename using packaging.utils
            # This properly handles multi-platform wheels like:
            # rpds_py-0.28.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
            try:
                _, _, _, wheel_tags = packaging.utils.parse_wheel_filename(tarball)
            except Exception as e:
                log.warning(f'Could not parse wheel filename {tarball}: {e}')
                continue
            
            # Check each tag in the wheel (multi-platform wheels have multiple tags)
            for wheel_tag in wheel_tags:
                wheel_tag_str = str(wheel_tag)
                if wheel_tag_str in tag_priority:
                    priority = tag_priority[wheel_tag_str]
                    if priority < best_priority:
                        best_match = tarball_info
                        best_priority = priority
                        log.debug(f'Found compatible wheel: {tarball} with tag {wheel_tag_str} (priority: {priority})')
                    break  # Found a match, no need to check other tags for this wheel
        
        if best_match:
            log.debug(f'Selected {best_match["tarball"]} with priority {best_priority}')
            return best_match
        
        # If no match found, return the first one (backward compatibility)
        log.warning(f'No exact platform match found for {self.name}, using first tarball')
        return self.__tarballs_info[0]

    @property
    def tarball_package(self):
        """
        Return the canonical package for the tarball

        This is almost always equal to ``self`` except if the package
        or the ``checksums.ini`` file is a symbolic link. In that case,
        the package of the symbolic link is returned.

        OUTPUT:

        A ``Package`` instance
        """
        n = self.__tarball_package_name
        if n == self.name:
            return self
        else:
            return type(self)(n)

    @property
    def version(self):
        """
        Return the version

        OUTPUT:

        String. The package version. Excludes the Sage-specific
        patchlevel.
        """
        return self.__version

    @property
    def version_major(self):
        """
        Return the major version

        OUTPUT:

        String. The package's major version.
        """
        return self.version.split('.')[0]

    @property
    def version_minor(self):
        """
        Return the minor version

        OUTPUT:

        String. The package's minor version.
        """
        return self.version.split('.')[1]

    @property
    def version_micro(self):
        """
        Return the micro version

        OUTPUT:

        String. The package's micro version.
        """
        return self.version.split('.')[2]

    @property
    def patchlevel(self):
        """
        Return the patchlevel

        OUTPUT:

        Integer. The patchlevel of the package. Excludes the "p"
        prefix.
        """
        return self.__patchlevel

    @property
    def version_with_patchlevel(self):
        """
        Return the version, including the Sage-specific patchlevel

        OUTPUT:

        String.
        """
        v = self.version
        if v is None:
            return v
        p = self.patchlevel
        if p < 0:
            return v
        return "{0}.p{1}".format(v, p)

    @property
    def type(self):
        """
        Return the package type
        """
        return self.__type

    @property
    def source(self):
        """
        Return the package source type
        """
        if self.__requirements is not None:
            return 'pip'
        if self.tarball_filename:
            if self.tarball_filename.endswith('.whl'):
                return 'wheel'
            return 'normal'
        if self.has_file('spkg-install') or self.has_file('spkg-install.in'):
            return 'script'
        return 'none'

    def is_platform_specific_wheel(self):
        """
        Check if this package uses a platform-specific wheel.
        
        Platform-specific wheels need special handling during download
        as they contain platform tags in the filename.
        
        Returns True if this is a wheel package with platform-specific binaries,
        False otherwise.
        """
        if self.source != 'wheel':
            return False
        if not self.tarball_filename:
            return False
        # Platform-independent wheels end with -none-any.whl
        return not self.tarball_filename.endswith('-none-any.whl')

    @property
    def trees(self):
        """
        Return the installation trees for the package

        OUTPUT:

        A white-space-separated string of environment variable names
        """
        if self.__trees is not None:
            return self.__trees
        if self.__version_requirements is not None:
            return 'SAGE_VENV'
        if self.__requirements is not None:
            return 'SAGE_VENV'
        return 'SAGE_LOCAL'

    @property
    def purl(self):
        """
        Return a PURL (Package URL) for the package

        OUTPUT:

        A string in the format ``SCHEME:TYPE/NAMESPACE/NAME``,
        i.e., without components for version, qualifiers, and subpath.
        See https://github.com/package-url/purl-spec/blob/master/PURL-SPECIFICATION.rst#package-url-specification-v10x
        for details
        """
        dist = self.distribution_name
        if dist:
            return 'pkg:pypi/' + dist.lower().replace('_', '-')
        return 'pkg:generic/' + self.name.replace('_', '-')

    @property
    def distribution_name(self):
        """
        Return the Python distribution name or ``None`` for non-Python packages
        """
        if self.__requirements is not None:
            for line in self.__requirements.split('\n'):
                line = line.strip()
                if line.startswith('#'):
                    continue
                for part in line.split():
                    return part
        if self.__version_requirements is None:
            return None
        for line in self.__version_requirements.split('\n'):
            line = line.strip()
            if line.startswith('#'):
                continue
            for part in line.split():
                return part
        return None

    @property
    def dependencies(self):
        """
        Return a list of strings, the package names of the (ordinary) dependencies
        """
        # after a '|', we have order-only dependencies
        return self.__dependencies.partition('|')[0].strip().split()

    @property
    def dependencies_order_only(self):
        """
        Return a list of strings, the package names of the order-only dependencies
        """
        return self.__dependencies.partition('|')[2].strip().split() + self.__dependencies_order_only.strip().split()

    @property
    def dependencies_optional(self):
        """
        Return a list of strings, the package names of the optional build dependencies
        """
        return self.__dependencies_optional.strip().split()

    @property
    def dependencies_runtime(self):
        """
        Return a list of strings, the package names of the runtime dependencies
        """
        # after a '|', we have order-only build dependencies
        return self.__dependencies.partition('|')[0].strip().split()

    @property
    def dependencies_check(self):
        """
        Return a list of strings, the package names of the check dependencies
        """
        return self.__dependencies_order_only.strip().split()

    def __eq__(self, other):
        return self.tarball == other.tarball

    @classmethod
    def all(cls):
        """
        Return all packages
        """
        base = os.path.join(SAGE_ROOT, 'build', 'pkgs')
        for subdir in os.listdir(base):
            path = os.path.join(base, subdir)
            if not os.path.isfile(os.path.join(path, "type")):
                log.debug('%s has no type', subdir)
                continue
            try:
                yield cls(subdir)
            except Exception:
                log.error('Failed to open %s', subdir)
                raise

    @property
    def path(self):
        """
        Return the package directory
        """
        return os.path.join(SAGE_ROOT, 'build', 'pkgs', self.name)

    def has_file(self, filename):
        """
        Return whether the file exists in the package directory
        """
        return os.path.exists(os.path.join(self.path, filename))

    def line_count_file(self, filename):
        """
        Return the number of lines of the file

        Directories are traversed recursively.

        OUTPUT:

        integer; 0 if the file cannot be read, 1 if it is a symlink
        """
        path = os.path.join(self.path, filename)
        if os.path.islink(path):
            return 1
        if os.path.isdir(path):
            return sum(self.line_count_file(os.path.join(filename, entry))
                       for entry in os.listdir(path))
        try:
            with open(path, "rb") as f:
                return len(list(f))
        except OSError:
            return 0

    def _init_checksum(self):
        """
        Load the checksums from the appropriate ``checksums.ini`` file
        
        Supports multiple tarballs with format:
        tarball=package-VERSION-cp311-cp311-manylinux_2_17_x86_64.whl
        sha256=abc123...
        upstream_url=https://...
        tarball=package-VERSION-cp311-cp311-macosx_11_0_arm64.whl
        sha256=def456...
        upstream_url=https://...
        """
        checksums_ini = os.path.join(self.path, 'checksums.ini')
        assignment = re.compile('(?P<var>[a-zA-Z0-9_]*)=(?P<value>.*)')
        
        # Store all entries, supporting multiple values for tarball, sha256, upstream_url
        tarballs = []
        sha256s = []
        sha1s = []
        upstream_urls = []
        
        try:
            with open(checksums_ini, 'rt') as f:
                for line in f.readlines():
                    match = assignment.match(line)
                    if match is None:
                        continue
                    var, value = match.groups()
                    
                    # Collect multiple entries
                    if var == 'tarball':
                        tarballs.append(value)
                    elif var == 'sha256':
                        sha256s.append(value)
                    elif var == 'sha1':
                        sha1s.append(value)
                    elif var == 'upstream_url':
                        upstream_urls.append(value)
        except IOError:
            pass
        
        # Store all tarballs info
        self.__tarballs_info = []
        for i, tarball in enumerate(tarballs):
            info = {
                'tarball': tarball,
                'sha256': sha256s[i] if i < len(sha256s) else None,
                'sha1': sha1s[i] if i < len(sha1s) else None,
                'upstream_url': upstream_urls[i] if i < len(upstream_urls) else None,
            }
            self.__tarballs_info.append(info)
        
        # For backward compatibility, set the first tarball as primary
        self.__sha1 = sha1s[0] if sha1s else None
        self.__sha256 = sha256s[0] if sha256s else None
        self.__tarball_pattern = tarballs[0] if tarballs else None
        self.__tarball_upstream_url_pattern = upstream_urls[0] if upstream_urls else None
        
        # Name of the directory containing the checksums.ini file
        self.__tarball_package_name = os.path.realpath(checksums_ini).split(os.sep)[-2]

    VERSION_PATCHLEVEL = re.compile(r'(?P<version>.*)\.p(?P<patchlevel>[0-9]+)')

    def _init_version(self):
        try:
            with open(os.path.join(self.path, 'package-version.txt')) as f:
                package_version = f.read().strip()
        except IOError:
            self.__version = None
            self.__patchlevel = None
        else:
            match = self.VERSION_PATCHLEVEL.match(package_version)
            if match is None:
                self.__version = package_version
                self.__patchlevel = -1
            else:
                self.__version = match.group('version')
                self.__patchlevel = int(match.group('patchlevel'))

    def _init_type(self):
        with open(os.path.join(self.path, 'type')) as f:
            package_type = f.read().strip()
        assert package_type in [
            'base', 'standard', 'optional', 'experimental'
        ]
        self.__type = package_type

    def _init_version_requirements(self):
        try:
            with open(os.path.join(self.path, 'version_requirements.txt')) as f:
                self.__version_requirements = f.read().strip()
        except IOError:
            self.__version_requirements = None

    def _init_requirements(self):
        try:
            with open(os.path.join(self.path, 'requirements.txt')) as f:
                self.__requirements = f.read().strip()
        except IOError:
            self.__requirements = None

    def _init_dependencies(self):
        try:
            with open(os.path.join(self.path, 'dependencies')) as f:
                self.__dependencies = f.readline().partition('#')[0].strip()
        except IOError:
            self.__dependencies = ''
        try:
            with open(os.path.join(self.path, 'dependencies_check')) as f:
                self.__dependencies_check = f.readline().partition('#')[0].strip()
        except IOError:
            self.__dependencies_check = ''
        try:
            with open(os.path.join(self.path, 'dependencies_optional')) as f:
                self.__dependencies_optional = f.readline().partition('#')[0].strip()
        except IOError:
            self.__dependencies_optional = ''
        try:
            with open(os.path.join(self.path, 'dependencies_order_only')) as f:
                self.__dependencies_order_only = f.readline().partition('#')[0].strip()
        except IOError:
            self.__dependencies_order_only = ''

    def _init_trees(self):
        try:
            with open(os.path.join(self.path, 'trees.txt')) as f:
                self.__trees = f.readline().partition('#')[0].strip()
        except IOError:
            self.__trees = None
