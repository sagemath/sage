# -*- coding: utf-8 -*-
"""
Controller for the commandline actions
"""


# ****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
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
from sage_bootstrap.tarball import Tarball, FileNotMirroredError
from sage_bootstrap.updater import ChecksumUpdater, PackageUpdater
from sage_bootstrap.creator import PackageCreator
from sage_bootstrap.pypi import PyPiVersion, PyPiNotFound, PyPiError
from sage_bootstrap.fileserver import FileServer
from sage_bootstrap.expand_class import PackageClass


class Application(object):

    def config(self):
        """
        Print the configuration

        $ sage --package config
        Configuration:
          * log = info
          * interactive = True
        """
        log.debug('Printing configuration')
        from sage_bootstrap.config import Configuration
        print(Configuration())

    def list_cls(self, *package_classes, **filters):
        """
        Print a list of all available packages

        $ sage --package list | sort
        4ti2
        arb
        autotools
        [...]
        zn_poly

        $ sage -package list --has-file=spkg-configure.m4 :experimental:
        perl_term_readline_gnu

        $ sage -package list --has-file=spkg-configure.m4 --has-file=distros/debian.txt | sort
        arb
        boost_cropped
        brial
        [...]
        zn_poly
        """
        log.debug('Listing packages')
        pc = PackageClass(*package_classes, **filters)
        for pkg_name in pc.names:
            print(pkg_name)

    def name(self, tarball_filename):
        """
        Find the package name given a tarball filename
    
        $ sage --package name pari-2.8-1564-gdeac36e.tar.gz
        pari
        """
        log.debug('Looking up package name for %s', tarball_filename)
        tarball = Tarball(os.path.basename(tarball_filename))
        print(tarball.package.name)

    def tarball(self, package_name):
        """
        Find the tarball filename given a package name
    
        $ sage --package tarball pari
        pari-2.8-1564-gdeac36e.tar.gz
        """
        log.debug('Looking up tarball name for %s', package_name)
        package = Package(package_name)
        print(package.tarball.filename)

    def apropos(self, incorrect_name):
        """
        Find up to 5 package names that are close to the given name

        $ sage --package apropos python
        Did you mean: cython, ipython, python2, python3, patch?
        """
        log.debug('Apropos for %s', incorrect_name)
        from sage_bootstrap.levenshtein import Levenshtein, DistanceExceeded
        levenshtein = Levenshtein(5)
        names = []
        for pkg in Package.all():
            try:
                names.append([levenshtein(pkg.name, incorrect_name), pkg.name])
            except DistanceExceeded:
                pass
        if names:
            names = sorted(names)[:5]
            print('Did you mean: {0}?'.format(', '.join(name[1] for name in names)))
        else:
            print('There is no package similar to {0}'.format(incorrect_name))
            print('You can find further packages at http://files.sagemath.org/spkg/')

    def update(self, package_name, new_version, url=None):
        """
        Update a package. This modifies the Sage sources. 
    
        $ sage --package update pari 2015 --url=http://localhost/pari/tarball.tgz
        """
        log.debug('Updating %s to %s', package_name, new_version)
        update = PackageUpdater(package_name, new_version)
        if url is not None or update.package.tarball_upstream_url:
            log.debug('Downloading %s', url)
            update.download_upstream(url)
        update.fix_checksum()

    def update_latest(self, package_name):
        """
        Update a package to the latest version. This modifies the Sage sources. 
        """
        try:
            pypi = PyPiVersion(package_name)
        except PyPiNotFound:
            log.debug('%s is not a pypi package', package_name)
            return
        else:
            pypi.update(Package(package_name))

    def update_latest_cls(self, package_name_or_class):
        exclude = [
            'cypari'   # Name conflict
        ]
        # Restrict to normal Python packages
        pc = PackageClass(package_name_or_class, has_files=['checksums.ini', 'install-requires.txt'])
        if not pc.names:
            log.warn('nothing to do (does not name a normal Python package)')
        for package_name in sorted(pc.names):
            if package_name in exclude:
                log.debug('skipping %s because of pypi name collision', package_name)
                continue
            try:
                self.update_latest(package_name)
            except PyPiError as e:
                log.warn('updating %s failed: %s', package_name, e)

    def download(self, package_name, allow_upstream=False):
        """
        Download a package

        $ sage --package download pari
        Using cached file /home/vbraun/Code/sage.git/upstream/pari-2.8-2044-g89b0f1e.tar.gz
        /home/vbraun/Code/sage.git/upstream/pari-2.8-2044-g89b0f1e.tar.gz
        """
        log.debug('Downloading %s', package_name)
        package = Package(package_name)
        package.tarball.download(allow_upstream=allow_upstream)
        print(package.tarball.upstream_fqn)

    def download_cls(self, package_name_or_class, allow_upstream=False, on_error='stop'):
        """
        Download a package or a class of packages
        """
        pc = PackageClass(package_name_or_class, has_files=['checksums.ini'])
        def download_with_args(package):
            try:
                self.download(package, allow_upstream=allow_upstream)
            except FileNotMirroredError:
                if on_error == 'stop':
                    raise
                elif on_error == 'warn':
                    log.warn('Unable to download tarball of %s', package)
                else:
                    raise ValueError('on_error must be one of "stop" and "warn"')
        pc.apply(download_with_args)

    def upload(self, package_name):
        """
        Upload a package to the Sage mirror network

        $ sage --package upload pari
        Uploading /home/vbraun/Code/sage.git/upstream/pari-2.8-2044-g89b0f1e.tar.gz
        """
        package = Package(package_name)
        if not os.path.exists(package.tarball.upstream_fqn):
            log.debug('Skipping %s because there is no local tarball', package_name)
            return
        if not package.tarball.is_distributable():
            log.info('Skipping %s because the tarball is marked as not distributable',
                     package_name)
            return
        log.info('Uploading %s', package.tarball.upstream_fqn)
        fs = FileServer()
        fs.upload(package)

    def upload_cls(self, package_name_or_class):
        pc = PackageClass(package_name_or_class)
        pc.apply(self.upload)
        fs = FileServer()
        log.info('Publishing')
        fs.publish()
        
    def fix_checksum_cls(self, *package_classes):
        """
        Fix the checksum of packages

        $ sage --package fix-checksum
        """
        pc = PackageClass(*package_classes, has_files=['checksums.ini'])
        pc.apply(self.fix_checksum)

    def fix_checksum(self, package_name):
        """
        Fix the checksum of a package

        $ sage --package fix-checksum pari
        Updating checksum of pari-2.8-2044-g89b0f1e.tar.gz
        """
        log.debug('Correcting the checksum of %s', package_name)
        update = ChecksumUpdater(package_name)
        pkg = update.package
        if not pkg.tarball_filename:
            log.info('Ignoring {0} because it is not a normal package'.format(package_name))
            return
        if not os.path.exists(pkg.tarball.upstream_fqn):
            log.info('Ignoring {0} because tarball is not cached'.format(package_name))
            return
        if pkg.tarball.checksum_verifies():
            log.info('Checksum of {0} (tarball {1}) unchanged'.format(package_name, pkg.tarball_filename))
        else:
            log.info('Updating checksum of {0} (tarball {1})'.format(package_name, pkg.tarball_filename))
            update.fix_checksum()

    def create(self, package_name, version=None, tarball=None, pkg_type=None, upstream_url=None,
               description=None, license=None, upstream_contact=None, pypi=False, source='normal'):
        """
        Create a normal package
        """
        if '-' in package_name:
            raise ValueError('package names must not contain dashes, use underscore instead')
        if pypi:
            pypi_version = PyPiVersion(package_name)
            if source == 'normal':
                if not tarball:
                    # Guess the general format of the tarball name.
                    tarball = pypi_version.tarball.replace(pypi_version.version, 'VERSION')
                if not version:
                    version = pypi_version.version
                # Use a URL from pypi.io instead of the specific URL received from the PyPI query
                # because it follows a simple pattern.
                upstream_url = 'https://pypi.io/packages/source/{0:1.1}/{0}/{1}'.format(package_name, tarball)
            if not description:
                description = pypi_version.summary
            if not license:
                license = pypi_version.license
            if not upstream_contact:
                upstream_contact = pypi_version.package_url
        if tarball and not pkg_type:
            # If we set a tarball, also make sure to create a "type" file,
            # so that subsequent operations (downloading of tarballs) work.
            pkg_type = 'optional'
        log.debug('Creating %s: %s, %s, %s', package_name, version, tarball, pkg_type)
        creator = PackageCreator(package_name)
        if version:
            creator.set_version(version)
        if pkg_type:
            creator.set_type(pkg_type)
        if description or license or upstream_contact:
            creator.set_description(description, license, upstream_contact)
        if pypi or source == 'pip':
            creator.set_python_data_and_scripts(pypi_package_name=pypi_version.name, source=source)
        if tarball:
            creator.set_tarball(tarball, upstream_url)
            if upstream_url and version:
                update = PackageUpdater(package_name, version)
                update.download_upstream()
            else:
                update = ChecksumUpdater(package_name)
            update.fix_checksum()
            
