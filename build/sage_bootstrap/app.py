# -*- coding: utf-8 -*-
"""
Controller for the commandline actions

AUTHORS:

    - Volker Braun (2016): initial version
    - Thierry Monteil (2022): clean option to remove outdated source tarballs
"""


# ****************************************************************************
#       Copyright (C) 2016      Volker Braun <vbraun.name@gmail.com>
#                     2020-2024 Matthias Koeppe
#                     2022      Thierry Monteil
#                     2024      Marc Culler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import os
import re
import logging
log = logging.getLogger()

from collections import defaultdict

from sage_bootstrap.package import Package
from sage_bootstrap.tarball import Tarball, FileNotMirroredError
from sage_bootstrap.updater import ChecksumUpdater, PackageUpdater
from sage_bootstrap.creator import PackageCreator
from sage_bootstrap.pypi import PyPiVersion, PyPiNotFound, PyPiError
from sage_bootstrap.fileserver import FileServer
from sage_bootstrap.expand_class import PackageClass
from sage_bootstrap.env import SAGE_DISTFILES


# Approximation of https://peps.python.org/pep-0508/#names dependency specification
dep_re = re.compile('^ *([-A-Z0-9._]+)', re.IGNORECASE)


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

        $ sage --package list
        4ti2
        _bootstrap
        _develop
        [...]
        zlib

        $ sage -package list --has-file=spkg-configure.m4 :experimental:
        perl_term_readline_gnu

        $ sage -package list --has-file=spkg-configure.m4 --has-file=distros/debian.txt
        4ti2
        _develop
        _prereq
        [...]
        zlib
        """
        log.debug('Listing packages')
        pc = PackageClass(*package_classes, **filters)
        for pkg_name in pc.names:
            print(pkg_name)

    def properties(self, *package_classes, **kwds):
        """
        Show the properties of given packages

        $ sage --package properties --format shell maxima
        path_maxima='........./build/pkgs/maxima'
        version_with_patchlevel_maxima='5.46.0'
        type_maxima='standard'
        source_maxima='normal'
        trees_maxima='SAGE_LOCAL'
        """
        props = kwds.pop('props', ['path', 'version_with_patchlevel', 'type', 'source', 'trees', 'purl'])
        format = kwds.pop('format', 'plain')
        log.debug('Looking up properties')
        pc = PackageClass(*package_classes)
        for package_name in pc.names:
            package = Package(package_name)
            if len(pc.names) > 1:
                if format == 'plain':
                    print("{0}:".format(package_name))
            for p in props:
                value = getattr(package, p)
                if value is None:
                    if p.startswith('version'):
                        value = 'none'
                    else:
                        value = ''
                if format == 'plain':
                    print("        {0:28} {1}".format(p + ":", value))
                else:
                    print("{0}_{1}='{2}'".format(p, package_name, value))

    def dependencies(self, *package_classes, **kwds):
        """
        Find the dependencies given package names

        $ sage --package dependencies maxima --runtime --order-only --format=shell
        order_only_deps_maxima='info'
        runtime_deps_maxima='ecl'
        """
        types = kwds.pop('types', None)
        format = kwds.pop('format', 'plain')
        log.debug('Looking up dependencies')
        pc = PackageClass(*package_classes)
        if format in ['plain', 'rst']:
            if types is None:
                typesets = [['order_only', 'runtime']]
            else:
                typesets = [[t] for t in types]
        elif format == 'shell':
            if types is None:
                types = ['order_only', 'optional', 'runtime', 'check']
            typesets = [[t] for t in types]
        else:
            raise ValueError('format must be one of "plain", "rst", and "shell"')

        for package_name in pc.names:
            package = Package(package_name)
            if len(pc.names) > 1:
                if format == 'plain':
                    print("{0}:".format(package_name))
                    indent1 = "        "
                elif format == 'rst':
                    print("\n{0}\n{1}\n".format(package_name, "~" * len(package_name)))
                    indent1 = ""
            else:
                indent1 = ""

            for typeset in typesets:
                if len(typesets) > 1:
                    if format == 'plain':
                        print(indent1 + "{0}: ".format('/'.join(typeset)))
                        indent2 = indent1 + "        "
                    elif format == 'rst':
                        print("\n" + indent1 + ".. tab:: {0}\n".format('/'.join(typeset)))
                        indent2 = indent1 + "    "
                else:
                    indent2 = indent1

                deps = []
                for t in typeset:
                    deps.extend(getattr(package, 'dependencies_' + t))
                deps = sorted(set(deps))

                if format in ['plain', 'rst']:
                    for dep in deps:
                        if '/' in dep:
                            # Suppress dependencies on source files, e.g. of the form $(SAGE_ROOT)/..., $(SAGE_SRC)/...
                            continue
                        if dep == 'FORCE':
                            # Suppress FORCE
                            continue
                        if dep.startswith('$('):
                            # Dependencies like $(BLAS)
                            print(indent2 + "- {0}".format(dep))
                        elif format == 'rst' and Package(dep).has_file('SPKG.rst'):
                            # This RST label is set in src/doc/bootstrap
                            print(indent2 + "- :ref:`spkg_{0}`".format(dep))
                        else:
                            print(indent2 + "- {0}".format(dep))
                elif format == 'shell':
                    # We single-quote the values because dependencies
                    # may contain Makefile variable substitutions
                    print("{0}_deps_{1}='{2}'".format(t, package_name, ' '.join(deps)))

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

    def commit(self, package_name, message=None):
        """
        Commit the changes to the Sage source tree for the given package
        """
        package = Package(package_name)
        if message is None:
            message = 'build/pkgs/{0}: Update to {1}'.format(package_name, package.version)
        os.system('git commit -m "{0}" {1}'.format(message, package.path))

    def update(self, package_name, new_version, url=None, commit=False):
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
        if commit:
            self.commit(package_name)

    def update_latest(self, package_name, commit=False):
        """
        Update a package to the latest version. This modifies the Sage sources.
        """
        pkg = Package(package_name)
        if pkg.source not in ['normal', 'wheel']:
            log.debug('update_latest can only update normal and wheel packages; %s is a %s package' % (pkg, pkg.source))
            return
        dist_name = pkg.distribution_name
        if dist_name is None:
            log.debug('%s does not have Python distribution info in version_requirements.txt' % pkg)
            return
        if pkg.tarball_pattern.endswith('.whl'):
            source = 'wheel'
        else:
            source = 'pypi'
        try:
            pypi = PyPiVersion(dist_name, source=source)
        except PyPiNotFound:
            log.debug('%s is not a pypi package', dist_name)
            return
        else:
            pypi.update(pkg)
        if commit:
            self.commit(package_name)

    def update_latest_cls(self, package_name_or_class, commit=False):
        exclude = [
            'cypari'   # Name conflict
        ]
        # Restrict to normal Python packages
        pc = PackageClass(package_name_or_class, has_files=['checksums.ini', 'version_requirements.txt'])
        if not pc.names:
            log.warn('nothing to do (does not name a normal Python package)')
        for package_name in sorted(pc.names):
            if package_name in exclude:
                log.debug('skipping %s because of pypi name collision', package_name)
                continue
            try:
                self.update_latest(package_name, commit=commit)
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

    def download_cls(self, *package_classes, **kwds):
        """
        Download a package or a class of packages
        """
        allow_upstream = kwds.pop('allow_upstream', False)
        on_error = kwds.pop('on_error', 'stop')
        has_files = list(kwds.pop('has_files', []))
        pc = PackageClass(*package_classes, has_files=has_files + ['checksums.ini'], **kwds)

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
        if pkg.tarball.checksum_verifies(force_sha256=True):
            log.info('Checksum of {0} (tarball {1}) unchanged'.format(package_name, pkg.tarball_filename))
        else:
            log.info('Updating checksum of {0} (tarball {1})'.format(package_name, pkg.tarball_filename))
            update.fix_checksum()

    def create(self, package_name, version=None, tarball=None, pkg_type=None, upstream_url=None,
               description=None, license=None, upstream_contact=None, pypi=False, source=None,
               dependencies=None):
        """
        Create a package

        $ sage --package create foo --version 1.3 --tarball FoO-VERSION.tar.gz --type experimental

        $ sage --package create scikit_spatial --pypi --type optional

        $ sage --package create torch --pypi --source pip --type optional

        $ sage --package create jupyterlab_markup --pypi --source wheel --type optional
        """
        if package_name.startswith('pypi/'):
            package_name = 'pkg:' + package_name
        if package_name.startswith('pkg:pypi/'):
            pypi = True
            package_name = package_name[len('pkg:pypi/'):].lower().replace('-', '_').replace('.', '_')
        elif '-' in package_name:
            raise ValueError('package names must not contain dashes, use underscore instead')
        if pypi:
            if source is None:
                try:
                    if PyPiVersion(package_name, source='wheel').tarball.endswith('-none-any.whl'):
                        source = 'wheel'
                    else:
                        source = 'normal'
                except PyPiError:
                    source = 'normal'
            pypi_version = PyPiVersion(package_name, source=source)
            if source == 'normal':
                if not tarball:
                    # Guess the general format of the tarball name.
                    tarball = pypi_version.tarball.replace(pypi_version.version, 'VERSION')
                if not version:
                    version = pypi_version.version
                # Use a URL from pypi.io instead of the specific URL received from the PyPI query
                # because it follows a simple pattern.
                upstream_url = 'https://pypi.io/packages/source/{0:1.1}/{0}/{1}'.format(package_name, tarball)
            elif source == 'wheel':
                if not tarball:
                    tarball = pypi_version.tarball.replace(pypi_version.version, 'VERSION')
                if not tarball.endswith('-none-any.whl'):
                    raise ValueError('Only platform-independent wheels can be used for wheel packages, got {0}'.format(tarball))
                if not version:
                    version = pypi_version.version
                if dependencies is None:
                    log.info('Requires-Python: {0}'.format(pypi_version.requires_python))
                    requires_dist = pypi_version.requires_dist
                    if requires_dist:
                        dependencies = []
                        for item in requires_dist:
                            if "extra ==" in item:
                                continue
                            try:
                                dep = dep_re.match(item).groups()[0].strip()
                            except Exception:
                                continue
                            dep = 'pkg:pypi/' + dep
                            try:
                                dep = Package(dep).name
                            except ValueError:
                                self.create(dep, pkg_type=pkg_type)
                                dep = Package(dep).name
                            dependencies.append(dep)
                upstream_url = 'https://pypi.io/packages/{2}/{0:1.1}/{0}/{1}'.format(package_name, tarball, pypi_version.python_version)
            if not description:
                description = pypi_version.summary
            if not license:
                license = pypi_version.license
            if not upstream_contact:
                upstream_contact = pypi_version.package_url
        if upstream_url and not tarball:
            tarball = upstream_url.rpartition('/')[2]
        if tarball and source is None:
            source = 'normal'
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
            creator.set_python_data_and_scripts(pypi_package_name=pypi_version.name, source=source,
                                                dependencies=dependencies)
        if tarball:
            creator.set_tarball(tarball, upstream_url)
            if upstream_url and version:
                update = PackageUpdater(package_name, version)
                update.download_upstream()
            else:
                update = ChecksumUpdater(package_name)
            update.fix_checksum()

    def clean(self):
        """
        Remove outdated source tarballs from the upstream/ directory

        $ sage --package clean
        42 files were removed from the .../upstream directory
        """
        log.debug('Cleaning upstream/ directory')
        package_names = PackageClass(':all:').names
        keep = [Package(package_name).tarball.filename for package_name in package_names]
        count = 0
        for filename in os.listdir(SAGE_DISTFILES):
            if filename not in keep:
                filepath = os.path.join(SAGE_DISTFILES, filename)
                if os.path.isfile(filepath):
                    log.debug('Removing file {}'.format(filepath))
                    os.remove(filepath)
                    count += 1
        print('{} files were removed from the {} directory'.format(count, SAGE_DISTFILES))

    def metrics_cls(self, *package_classes):
        """
        Show the metrics of given packages

        $ sage --package metrics :standard:
        has_file_distros_arch_txt=131
        has_file_distros_conda_txt=216
        has_file_distros_debian_txt=125
        has_file_distros_fedora_txt=138
        has_file_distros_gentoo_txt=181
        has_file_distros_homebrew_txt=61
        has_file_distros_macports_txt=129
        has_file_distros_nix_txt=51
        has_file_distros_opensuse_txt=146
        has_file_distros_slackware_txt=25
        has_file_distros_void_txt=184
        has_file_patches=35
        has_file_spkg_check=59
        has_file_spkg_configure_m4=222
        has_file_spkg_install=198
        has_tarball_upstream_url=231
        line_count_file_patches=22561
        line_count_file_spkg_check=402
        line_count_file_spkg_configure_m4=2792
        line_count_file_spkg_install=2960
        packages=272
        type_standard=272
        """
        log.debug('Computing metrics')
        metrics = defaultdict(int)
        pc = PackageClass(*package_classes)
        for package_name in pc.names:
            package = Package(package_name)
            metrics['packages'] += 1
            metrics['type_' + package.type] += 1
            for filenames in [['spkg-configure.m4'],
                              ['spkg-install', 'spkg-install.in'],
                              ['spkg-check', 'spkg-check.in'],
                              ['distros/arch.txt'],
                              ['distros/conda.txt'],
                              ['distros/debian.txt'],
                              ['distros/fedora.txt'],
                              ['distros/gentoo.txt'],
                              ['distros/homebrew.txt'],
                              ['distros/macports.txt'],
                              ['distros/nix.txt'],
                              ['distros/opensuse.txt'],
                              ['distros/slackware.txt'],
                              ['distros/void.txt'],
                              ['patches']]:
                key = filenames[0].replace('.', '_').replace('-', '_').replace('/', '_')
                metrics['has_file_' + key] += int(any(package.has_file(filename)
                                                      for filename in filenames))
                if not key.startswith('distros_'):
                    metrics['line_count_file_' + key] += sum(package.line_count_file(filename)
                                                             for filename in filenames)
            metrics['has_tarball_upstream_url'] += int(bool(package.tarball_upstream_url))
        for key, value in sorted(metrics.items()):
            print('{0}={1}'.format(key, value))
