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

        If the license field contains overly long text (the full license text),
        try to extract a short license identifier from the classifiers instead.
        """
        license_text = self.json['info']['license']
        license_expression = self.json['info'].get('license_expression')

        # If there's a license expression (PEP 639), prefer that
        if license_expression:
            return license_expression

        # If the license text is short enough, use it directly
        if license_text and len(license_text) <= 100:
            return license_text

        # Try to extract license from classifiers
        classifiers = self.json['info'].get('classifiers', [])
        license_classifiers = []
        for classifier in classifiers:
            if classifier.startswith('License :: '):
                # Extract just the license name from the classifier
                # e.g., "License :: OSI Approved :: BSD License" -> "BSD License"
                parts = classifier.split(' :: ')
                if len(parts) >= 3:
                    license_classifiers.append(parts[-1])
                elif len(parts) == 2:
                    license_classifiers.append(parts[-1])

        if license_classifiers:
            return ', '.join(license_classifiers)

        # If we have a long license text but no classifiers, try to extract
        # a short identifier from the first line or truncate
        if license_text:
            first_line = license_text.split('\n')[0].strip()
            # Check if first line looks like a license name (short and descriptive)
            if len(first_line) <= 100:
                return first_line
            # Otherwise, check for common license patterns in the text
            license_patterns = [
                ('BSD 3-Clause', ['BSD 3-Clause', 'BSD-3-Clause', 'three conditions']),
                ('BSD 2-Clause', ['BSD 2-Clause', 'BSD-2-Clause', 'two conditions', 'Simplified BSD']),
                ('MIT', ['MIT License', 'Permission is hereby granted, free of charge']),
                ('Apache 2.0', ['Apache License', 'Version 2.0']),
                ('GPL', ['GNU General Public License']),
                ('LGPL', ['GNU Lesser General Public License']),
            ]
            for short_name, patterns in license_patterns:
                if all(p in license_text for p in patterns[:1]) or \
                   (len(patterns) > 1 and patterns[1] in license_text):
                    return short_name

        return license_text

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
