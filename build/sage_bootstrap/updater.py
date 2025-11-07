# -*- coding: utf-8 -*-
"""
Package Updater
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os

import logging
log = logging.getLogger()

from sage_bootstrap.package import Package
from sage_bootstrap.download import Download


class ChecksumUpdater(object):

    def __init__(self, package_name):
        self.__package = None
        self.package_name = package_name

    @property
    def package(self):
        if self.__package is None:
            self.__package = Package(self.package_name)
        return self.__package

    def fix_checksum(self):
        checksums_ini = os.path.join(self.package.path, 'checksums.ini')
        s = self.checksums_ini()
        with open(checksums_ini, 'w') as f:
            f.write(s)

    def checksums_ini(self):
        tarball = self.package.tarball
        result = [
            'tarball=' + self.package.tarball_pattern,
            'sha1=' + tarball._compute_sha1(),
            'sha256=' + tarball._compute_sha256()
        ]
        if self.package.tarball_upstream_url_pattern:
            result.append('upstream_url=' + self.package.tarball_upstream_url_pattern)
        result.append('')   # newline at end
        return '\n'.join(result)


class PackageUpdater(ChecksumUpdater):

    def __init__(self, package_name, new_version):
        super(PackageUpdater, self).__init__(package_name)
        self._update_version(new_version)

    def _update_version(self, new_version):
        old = Package(self.package_name)
        package_version_txt = os.path.join(old.path, 'package-version.txt')
        with open(package_version_txt, 'w') as f:
            f.write(new_version.strip() + '\n')

    def download_upstream(self, download_url=None):
        tarball = self.package.tarball
        if download_url is None:
            pattern = self.package.tarball_upstream_url_pattern
            if pattern and 'VERSION' not in pattern:
                print('Warning: upstream_url pattern does not use the VERSION variable')
            download_url = self.package.tarball_upstream_url
        if download_url is None:
            raise ValueError("package has no default upstream_url pattern, download_url needed")
        print('Downloading tarball from {0} to {1}'.format(download_url, tarball.upstream_fqn))
        Download(download_url, tarball.upstream_fqn).run()
