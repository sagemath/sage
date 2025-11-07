# -*- coding: utf-8 -*-
"""
Third-Party Tarballs
"""

# ****************************************************************************
#       Copyright (C) 2014-2015 Volker Braun <vbraun.name@gmail.com>
#                     2017      Jeroen Demeyer
#                     2020      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import sys
import subprocess
import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.download import Download, MirrorList
from sage_bootstrap.package import Package


class ChecksumError(Exception):
    """
    Exception raised when the checksum of the tarball does not match
    """
    pass


class FileNotMirroredError(Exception):
    """
    Exception raised when the tarball cannot be downloaded from the mirrors
    """
    pass


class Tarball(object):

    def __init__(self, tarball_name, package=None, tarball_info=None):
        """
        A (third-party downloadable) tarball

        Note that the tarball might also be a different kind of
        archive format that is supported, it does not necessarily have
        to be tar.

        INPUT:

        - ``tarball_name`` - string. The full filename (``foo-1.3.tar.bz2``)
          of a tarball on the Sage mirror network.
        - ``package`` - Package object, or None to auto-detect
        - ``tarball_info`` - dict with tarball info (for multi-tarball packages)
          containing sha256, sha1, upstream_url
        """
        self.__filename = tarball_name
        self.__tarball_info = tarball_info
        
        if package is None:
            self.__package = None
            for pkg in Package.all():
                if pkg.tarball_filename == tarball_name:
                    self.__package = pkg.tarball_package
            if self.package is None:
                error = 'tarball {0} is not referenced by any Sage package'.format(tarball_name)
                log.error(error)
                raise ValueError(error)
        else:
            self.__package = package
            # For multi-tarball packages (with tarball_info), skip the filename check
            # since we're selecting a platform-specific wheel
            if tarball_info is None and package.tarball_filename != tarball_name:
                error = 'tarball {0} is not referenced by the {1} package'.format(tarball_name, package.name)
                log.error(error)
                raise ValueError(error)

    def __repr__(self):
        return 'Tarball {0}'.format(self.filename)

    @property
    def filename(self):
        """
        Return the tarball filename

        OUTPUT:

        String. The full filename (``foo-1.3.tar.bz2``) of the
        tarball.
        """
        return self.__filename

    @property
    def package(self):
        """
        Return the package that the tarball belongs to

        OUTPUT:

        Instance of :class:`sage_bootstrap.package.Package`
        """
        return self.__package

    @property
    def upstream_fqn(self):
        """
        The fully-qualified (including directory) file name in the upstream directory.
        """
        return os.path.join(SAGE_DISTFILES, self.filename)

    def __eq__(self, other):
        return self.filename == other.filename

    def _compute_hash(self, algorithm):
        with open(self.upstream_fqn, 'rb') as f:
            while True:
                buf = f.read(0x100000)
                if not buf:
                    break
                algorithm.update(buf)
        return algorithm.hexdigest()

    def _compute_sha1(self):
        import hashlib
        return self._compute_hash(hashlib.sha1())

    def _compute_sha256(self):
        import hashlib
        return self._compute_hash(hashlib.sha256())

    def checksum_verifies(self, force_sha256=False):
        """
        Test whether the checksum of the downloaded file is correct.
        
        Uses tarball_info if available (for multi-tarball packages),
        otherwise falls back to package-level checksums.
        """
        # Use tarball_info if available (for multi-tarball packages)
        if self.__tarball_info:
            sha256_expected = self.__tarball_info.get('sha256')
            sha1_expected = self.__tarball_info.get('sha1')
        else:
            sha256_expected = self.package.sha256
            sha1_expected = self.package.sha1
        
        if sha256_expected:
            sha256 = self._compute_sha256()
            if sha256 != sha256_expected:
                log.error(f'SHA256 mismatch for {self.filename}')
                log.error(f'Expected: {sha256_expected}')
                log.error(f'Got:      {sha256}')
                return False
        elif force_sha256:
            log.warning('sha256 not available for {0}'.format(self.package.name))
            return False
        else:
            log.warning('sha256 not available for {0}, using sha1'.format(self.package.name))
            if sha1_expected:
                sha1 = self._compute_sha1()
                if sha1 != sha1_expected:
                    log.error(f'SHA1 mismatch for {self.filename}')
                    log.error(f'Expected: {sha1_expected}')
                    log.error(f'Got:      {sha1}')
                    return False
                return True
            else:
                log.warning('No checksum available for {0}'.format(self.package.name))
                return False
        
        return True

    def is_distributable(self):
        return 'do-not-distribute' not in self.filename

    def is_platform_specific_wheel(self):
        """
        Check if this is a platform-specific wheel.
        
        Platform-specific wheels have platform tags like:
        - manylinux_2_17_x86_64
        - macosx_11_0_arm64
        - win_amd64
        
        Platform-independent wheels end with -none-any.whl
        """
        if not self.filename or not self.filename.endswith('.whl'):
            return False
        return not self.filename.endswith('-none-any.whl')

    def _find_cached_wheel_for_platform(self):
        """
        Find a cached wheel file that matches the current platform.
        
        Looks in SAGE_DISTFILES for any wheel matching the package name,
        version, and current platform/Python version.
        
        Returns the path to the cached wheel if found, None otherwise.
        """
        import platform
        
        # Get platform info
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        py_tag = f"cp{python_version.replace('.', '')}"
        
        system = platform.system().lower()
        machine = platform.machine().lower()
        
        # Build possible platform tags
        platform_patterns = []
        if system == 'linux':
            if machine == 'x86_64':
                platform_patterns = ['manylinux', 'linux_x86_64']
            elif machine in ['aarch64', 'arm64']:
                platform_patterns = ['manylinux', 'linux_aarch64']
        elif system == 'darwin':
            if machine in ['arm64', 'aarch64']:
                platform_patterns = ['macosx', 'arm64']
            else:
                platform_patterns = ['macosx', 'x86_64']
        
        # Look for matching wheel files in upstream directory
        # Note: Wheel filenames use underscores, not dashes
        pkg_name_wheel = self.package.name.replace('-', '_').lower()
        pkg_name_pypi = self.package.name.replace('_', '-').lower()
        pkg_version = self.package.version
        
        log.debug(f'Looking for cached wheel: pkg_name_wheel={pkg_name_wheel}, pkg_name_pypi={pkg_name_pypi}, version={pkg_version}, py_tag={py_tag}')
        log.debug(f'Platform patterns: {platform_patterns}')
        
        try:
            for filename in os.listdir(SAGE_DISTFILES):
                if not filename.endswith('.whl'):
                    continue
                
                filename_lower = filename.lower()
                
                log.debug(f'Checking file: {filename}')
                
                # Check if matches package name and version
                # Wheel filenames use underscores, not dashes
                if not (filename_lower.startswith(pkg_name_wheel) or filename_lower.startswith(pkg_name_pypi)):
                    log.debug(f'  Does not start with {pkg_name_wheel} or {pkg_name_pypi}')
                    continue
                if pkg_version not in filename:
                    log.debug(f'  Does not contain version {pkg_version}')
                    continue
                
                # Check if matches current platform
                matches_platform = any(p in filename_lower for p in platform_patterns)
                matches_python = py_tag in filename
                
                log.debug(f'  matches_platform={matches_platform}, matches_python={matches_python}')
                
                if matches_platform and matches_python:
                    wheel_path = os.path.join(SAGE_DISTFILES, filename)
                    log.info(f'Found cached wheel for platform: {filename}')
                    return wheel_path
        except OSError as e:
            log.warning(f'Error listing directory {SAGE_DISTFILES}: {e}')
            pass
        
        log.warning(f'No cached wheel found for {pkg_name_wheel}/{pkg_name_pypi}-{pkg_version} (py_tag={py_tag}, platform={platform_patterns})')
        return None

    def download(self, allow_upstream=False):
        """
        Download the tarball to the upstream directory.

        If allow_upstream is False and the package cannot be found
        on the sage mirrors, fall back to downloading it from
        the upstream URL if the package has one.
        
        For platform-specific wheels, this method:
        1. Checks for a cached wheel matching the current platform
        2. If cached and checksum valid, uses it
        3. Otherwise, uses pip download to get the correct wheel
        4. Falls back to traditional download if pip fails
        """
        if not self.filename:
            raise ValueError('non-normal package does define a tarball, so cannot download')
        
        destination = self.upstream_fqn
        
        # Check if package has multiple tarballs (multi-platform wheels)
        has_multiple_tarballs = len(self.package.tarballs_info) > 1
        
        if has_multiple_tarballs:
            log.info(f'Package {self.package.name} has {len(self.package.tarballs_info)} platform-specific tarballs')
            return self._download_multiple_wheels(allow_upstream)
        
        # Single tarball case - existing logic
        # Check if file already exists and is valid
        if os.path.isfile(destination):
            if self.checksum_verifies():
                log.info('Using cached file {destination}'.format(destination=destination))
                return
            else:
                # Garbage in the upstream directory? Ignore it.
                # Don't delete it because maybe somebody just forgot to
                # update the checksum (Issue #23972).
                log.warning('Invalid checksum; ignoring cached file {destination}'
                            .format(destination=destination))
                
        # Traditional download logic for tarballs and platform-independent wheels
        successful_download = False
        log.info('Attempting to download package {0} from mirrors'.format(self.filename))
        for mirror in MirrorList():
            url = mirror.replace('${SPKG}', self.package.name)
            if not url.endswith('/'):
                url += '/'
            url += self.filename
            log.info(url)
            try:
                Download(url, destination).run()
                successful_download = True
                break
            except IOError:
                log.debug('File not on mirror')
        if not successful_download:
            url = self.package.tarball_upstream_url
            if allow_upstream and url:
                log.info('Attempting to download from {}'.format(url))
                try:
                    Download(url, destination).run()
                except IOError:
                    raise FileNotMirroredError('tarball does not exist on mirror network and neither at the upstream URL')
            else:
                raise FileNotMirroredError('tarball does not exist on mirror network')
        if not self.checksum_verifies():
            raise ChecksumError('checksum does not match')

    def _download_multiple_wheels(self, allow_upstream=False):
        """
        Handle download for packages with multiple platform-specific wheels.
        
        Strategy:
        1. Check if already cached with valid checksum
        2. If not cached, try pip download (auto-detects platform/Python version)
        3. If pip fails (offline/error), fall back to traditional download from mirrors/upstream
        """
        # Find the appropriate tarball for this platform from checksums.ini
        tarball_info = self.package.find_tarball_for_platform()
        
        if not tarball_info:
            raise ValueError(f'No suitable tarball found for {self.package.name} on current platform')
        
        # Get the actual filename with version substituted
        tarball_pattern = tarball_info['tarball']
        tarball_filename = self.package._substitute_variables(tarball_pattern)
        
        log.info(f'Selected tarball for platform: {tarball_filename}')
        
        # Step 1: Check cache first
        cached_wheel = self._find_cached_wheel_for_platform()
        if cached_wheel:
            # Verify checksum of cached wheel
            try:
                cached_tarball = Tarball(os.path.basename(cached_wheel), 
                                        package=self.package, 
                                        tarball_info=tarball_info)
                if cached_tarball.checksum_verifies():
                    log.info(f'Using cached wheel with valid checksum: {os.path.basename(cached_wheel)}')
                    # Update self to point to the cached wheel
                    self.__filename = os.path.basename(cached_wheel)
                    return
                else:
                    log.warning('Cached wheel has invalid checksum, will re-download')
            except Exception as e:
                log.warning(f'Error checking cached wheel: {e}')
        
        dest_dir = SAGE_DISTFILES

        # Step 2: Try to download from mirrors/upstream
        log.info(f'Downloading {tarball_filename} using traditional method (mirrors/upstream)...')
        destination = os.path.join(dest_dir, tarball_filename)
        upstream_url_pattern = tarball_info.get('upstream_url')
        
        if upstream_url_pattern:
            upstream_url = self.package._substitute_variables(upstream_url_pattern)
        else:
            upstream_url = None
        
        successful_download = False
        
        # Try mirrors first
        for mirror in MirrorList():
            url = mirror.replace('${SPKG}', self.package.name)
            if not url.endswith('/'):
                url += '/'
            url += tarball_filename
            log.debug(f'Trying mirror: {url}')
            try:
                Download(url, destination).run()
                successful_download = True
                log.info(f'Downloaded from mirror: {url}')
                break
            except IOError:
                log.debug('File not on mirror')
        
        # Try upstream if mirrors failed
        if not successful_download and upstream_url and allow_upstream:
            log.info(f'Trying upstream: {upstream_url}')
            try:
                Download(upstream_url, destination).run()
                successful_download = True
                log.info(f'Downloaded from upstream: {upstream_url}')
            except IOError:
                log.debug('File not at upstream URL')
        
        if not successful_download:
            raise FileNotMirroredError(f'Could not download {tarball_filename} from pip, mirrors, or upstream')
        
        # Verify checksum of traditionally-downloaded file
        try:
            downloaded_tarball = Tarball(tarball_filename,
                                        package=self.package,
                                        tarball_info=tarball_info)
            if not downloaded_tarball.checksum_verifies():
                raise ChecksumError(f'Checksum verification failed for {tarball_filename}')
            log.info(f'Successfully downloaded and verified: {tarball_filename}')
            # Update self to point to the actually downloaded file
            self.__filename = tarball_filename
        except Exception as e:
            log.error(f'Error verifying downloaded tarball: {e}')
            raise

    def save_as(self, destination):
        """
        Save the tarball as a new file
        """
        import shutil
        shutil.copy(self.upstream_fqn, destination)
