# -*- coding: utf-8 -*-
"""
PyPi Version Information
"""


# ****************************************************************************
#       Copyright (C) 2016      Volker Braun <vbraun.name@gmail.com>
#                     2020-2023 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import logging
log = logging.getLogger()

import json

from sage_bootstrap.package import Package
from sage_bootstrap.updater import PackageUpdater
from sage_bootstrap.compat import urllib


class PyPiNotFound(Exception):
    pass


class PyPiError(Exception):
    pass


class PyPiVersion(object):

    def __init__(self, package_name, source='normal'):
        self.name = package_name
        self.json = self._get_json()
        # Replace provided name with the canonical name
        self.name = self.json['info']['name']
        if source == 'wheel':
            self.python_version = 'py3'
        else:
            self.python_version = 'source'

    def _get_json(self):
        response = urllib.urlopen(self.json_url)
        if response.getcode() != 200:
            raise PyPiNotFound('%s not on pypi', self.name)
        data = response.read()
        text = data.decode('utf-8')
        return json.loads(text)

    @property
    def json_url(self):
        return 'https://pypi.python.org/pypi/{0}/json'.format(self.name)

    @property
    def version(self):
        """
        Return the current version
        """
        return self.json['info']['version']

    @property
    def url(self):
        """
        Return the source url
        """
        for download in self.json['urls']:
            if self.python_version in download['python_version']:
                self.python_version = download['python_version']
                return download['url']
        raise PyPiError('No %s url for %s found', self.python_version, self.name)

    @property
    def tarball(self):
        """
        Return the source tarball name
        """
        for download in self.json['urls']:
            if self.python_version in download['python_version']:
                self.python_version = download['python_version']
                return download['filename']
        raise PyPiError('No %s url for %s found', self.python_version, self.name)

    @property
    def package_url(self):
        """
        Return the package URL
        """
        return self.json['info']['package_url']

    @property
    def license(self):
        """
        Return the package license
        """
        return self.json['info']['license']

    @property
    def summary(self):
        """
        Return the package summary
        """
        return self.json['info']['summary']

    @property
    def requires_dist(self):
        """
        Return the dependencies
        """
        return self.json['info']['requires_dist']

    @property
    def requires_python(self):
        """
        Return the requires_python attribute
        """
        return self.json['info']['requires_python']

    def update(self, package=None):
        if package is None:
            package = Package(self.name)
        if package.version == self.version:
            log.info('%s is already at the latest version', self.name)
            return
        log.info('Updating %s: %s -> %s', package.name, package.version, self.version)
        update = PackageUpdater(package.name, self.version)
        update.download_upstream(self.url)
        update.fix_checksum()
