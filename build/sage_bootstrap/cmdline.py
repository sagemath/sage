# -*- coding: utf-8 -*-
"""
View for the Commandline UI

This module handles the main "sage-package" commandline utility, which
is also exposed as "sage --package".

AUTHORS:

    - Volker Braun (2016): initial version
    - Thierry Monteil (2022): clean option to remove outdated source tarballs
"""

# ****************************************************************************
#       Copyright (C) 2015-2016 Volker Braun <vbraun.name@gmail.com>
#                     2020-2024 Matthias Koeppe
#                     2022      Thierry Monteil
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys
import logging
log = logging.getLogger()


import argparse

from sage_bootstrap.app import Application


description = \
"""
SageMath Bootstrap Library

Provides scripts to manage the packages of Sage-the-distribution,
including SageMath's database of equivalent system packages,
and to download and upload tarballs from/to SageMath servers.
"""


epilog = \
"""
The individual subcommands have their own detailed help, for example
run "sage --package config -h" to see the help on the config option.
"""


epilog_config = \
"""
Print the configuration

EXAMPLE:

    $ sage --package config
    Configuration:
    * log = info
    * interactive = True
"""


epilog_list = \
"""
Print a list of packages known to Sage (sorted alphabetically)

EXAMPLE:

    $ sage --package list
    4ti2
    [...]
    zlib

    $ sage --package list :standard:
    _prereq
    [...]
    zlib
"""


epilog_properties = \
"""
Print properties of given package.

EXAMPLE:

    $ sage --package properties maxima
    maxima:
            path:                        /.../build/pkgs/maxima
            version_with_patchlevel:     5.46.0
            type:                        standard
            source:                      normal
            trees:                       SAGE_LOCAL
"""


epilog_dependencies = \
"""
Print the list of packages that are dependencies of given package.
By default, list a summary of the build, order-only, and runtime
dependencies.

EXAMPLE:

    $ sage --package dependencies maxima openblas
    maxima:
            - ecl
            - info
    openblas:
            - gfortran
    $ sage --package dependencies maxima --runtime
    - ecl

    $ sage --package dependencies maxima openblas --runtime --order-only
    maxima:
            order_only:
                    - info
            runtime:
                    - ecl
    openblas:
            order_only:
            runtime:
                    - gfortran
"""


epilog_name = \
"""
Find the package name given a tarball filename

EXAMPLE:

    $ sage --package name pari-2.8-1564-gdeac36e.tar.gz
    pari
"""


epilog_tarball = \
"""
Find the tarball filename given a package name

EXAMPLE:

    $ sage --package tarball pari
    pari-2.8-1564-gdeac36e.tar.gz
"""


epilog_apropos = \
"""
Find up to 5 package names that are close to the given name

EXAMPLE:

    $ sage --package apropos python
    Did you mean: cython, ipython, python2, python3, patch?
"""


epilog_update = \
"""
Update a package. This modifies the Sage sources.

EXAMPLE:

    $ sage --package update pari 2015 --url=http://localhost/pari/tarball.tgz
"""


epilog_update_latest = \
"""
Update a package to the latest version. This modifies the Sage sources.

EXAMPLE:

    $ sage --package update-latest ipython
"""


epilog_download = \
"""
Download the tarball for a package and print the filename to stdout

EXAMPLE:

    $ sage --package download pari
    Using cached file /home/vbraun/Code/sage.git/upstream/pari-2.8-2044-g89b0f1e.tar.gz
    /home/vbraun/Code/sage.git/upstream/pari-2.8-2044-g89b0f1e.tar.gz
"""


epilog_upload = \
"""
Upload the tarball to the Sage mirror network (requires ssh key authentication)

EXAMPLE:

    $ sage --package upload pari
    Uploading /home/vbraun/Code/sage.git/upstream/pari-2.8-2044-g89b0f1e.tar.gz
"""


epilog_fix_checksum = \
"""
Fix the checksum of a package

EXAMPLE:

    $ sage --package fix-checksum pari
    Updating checksum of pari (tarball pari-2.8-2044-g89b0f1e.tar.gz)
"""

epilog_create = \
"""
Create new package, or overwrite existing package

EXAMPLE:

    $ sage --package create foo --version=3.14 --tarball=Foo-VERSION.tar.bz2 --type=standard
    Creating new package "foo"
"""

epilog_clean = \
"""
Remove outdated source tarballs from the upstream/ directory

EXAMPLE:

    $ sage --package clean
    42 files were removed from the .../upstream directory
"""

epilog_metrics = \
"""
Print metrics of given package.

EXAMPLE:

    $ sage --package metrics :standard:
    has_file_distros_arch_txt=131
    ...
    has_file_distros_void_txt=184
    has_file_patches=35
    has_file_spkg_check=59
    has_file_spkg_configure_m4=222
    has_file_spkg_install=198
    has_tarball_upstream_url=231
    line_count_file_patches=22561
    ...
    packages=272
    type_standard=272
"""


def make_parser():
    """
    The main commandline argument parser
    """
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        prog='sage --package',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--log', dest='log', default=None,
                        help='one of [DEBUG, INFO, ERROR, WARNING, CRITICAL]')
    subparsers = parser.add_subparsers(dest='subcommand')

    parser_config = subparsers.add_parser(
        'config', epilog=epilog_config,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='print the configuration')

    parser_list = subparsers.add_parser(
        'list', epilog=epilog_list,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='print a list of packages known to Sage')
    parser_list.add_argument(
        'package_class', metavar='[PACKAGE_NAME|pkg:pypi/DISTRIBUTION-NAME|:PACKAGE_TYPE:]',
        type=str, default=[':all-or-nothing:'], nargs='*',
        help=('package name, pkg:pypi/ followed by a distribution name, '
              'or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, and :experimental:); '
              'default: :all: (or nothing when --include-dependencies or --exclude-dependencies is given)'))
    parser_list.add_argument(
        '--has-file', action='append', default=[], metavar='FILENAME', dest='has_files',
        help=('only include packages that have this file in their metadata directory '
              '(examples: SPKG.rst, spkg-configure.m4, distros/debian.txt, spkg-install|spkg-install.in)'))
    parser_list.add_argument(
        '--no-file', action='append', default=[], metavar='FILENAME', dest='no_files',
        help=('only include packages that do not have this file in their metadata directory '
              '(examples: huge, patches, huge|has_nonfree_dependencies)'))
    parser_list.add_argument(
        '--exclude', nargs='*', action='append', default=[], metavar='PACKAGE_NAME',
        help='exclude package from list')
    parser_list.add_argument(
        '--include-dependencies', action='store_true',
        help='include (ordinary) dependencies of the packages recursively')
    parser_list.add_argument(
        '--exclude-dependencies', action='store_true',
        help='exclude (ordinary) dependencies of the packages recursively')

    parser_properties = subparsers.add_parser(
        'properties', epilog=epilog_properties,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='print properties of given packages')
    parser_properties.add_argument(
        'package_class', metavar='[PACKAGE_NAME|pkg:pypi/DISTRIBUTION-NAME|:PACKAGE_TYPE:]',
        type=str, nargs='+',
        help=('package name, pkg:pypi/ followed by a distribution name, '
              'or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, and :experimental:)'))
    parser_properties.add_argument(
        '--format', type=str, default='plain',
        help='output format (one of plain and shell; default: plain)')

    parser_dependencies = subparsers.add_parser(
        'dependencies', epilog=epilog_dependencies,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='print the list of packages that are dependencies of given packages')
    parser_dependencies.add_argument(
        'package_class', metavar='[package_name|:package_type:]',
        type=str, nargs='+',
        help=('package name or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, and :experimental:)'))
    parser_dependencies.add_argument(
        '--order-only', action='store_true',
        help='list the order-only build dependencies')
    parser_dependencies.add_argument(
        '--optional', action='store_true',
        help='list the optional build dependencies')
    parser_dependencies.add_argument(
        '--runtime', action='store_true',
        help='list the runtime dependencies')
    parser_dependencies.add_argument(
        '--check', action='store_true',
        help='list the check dependencies')
    parser_dependencies.add_argument(
        '--format', type=str, default='plain',
        help='output format (one of plain, rst, and shell; default: plain)')

    parser_name = subparsers.add_parser(
        'name', epilog=epilog_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='find the package name given a tarball filename')
    parser_name.add_argument('tarball_filename', type=str, help='tarball filename')

    parser_tarball = subparsers.add_parser(
        'tarball', epilog=epilog_tarball,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='find the tarball filename given a package name')
    parser_tarball.add_argument('package_name', type=str, help='package name')

    parser_apropos = subparsers.add_parser(
        'apropos', epilog=epilog_apropos,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='find up to 5 package names that are close to the given name')
    parser_apropos.add_argument(
        'incorrect_name', type=str,
        help='fuzzy name to search for')

    parser_update = subparsers.add_parser(
        'update', epilog=epilog_update,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='update a package, modifying the Sage sources')
    parser_update.add_argument(
        'package_name', type=str, help='package name')
    parser_update.add_argument(
        'new_version', type=str, help='new version')
    parser_update.add_argument(
        '--url', type=str, default=None, help='download URL')
    parser_update.add_argument(
        '--commit', action="store_true",
        help='whether to run "git commit"')

    parser_update_latest = subparsers.add_parser(
        'update-latest', epilog=epilog_update_latest,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='update a package to the latest version, modifying the Sage sources')
    parser_update_latest.add_argument(
        'package_name', type=str, help='package name (:all: for all packages)')
    parser_update_latest.add_argument(
        '--commit', action="store_true",
        help='whether to run "git commit"')

    parser_download = subparsers.add_parser(
        'download', epilog=epilog_download,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='download tarball')
    parser_download.add_argument(
        'package_class', metavar='[package_name|:package_type:]',
        type=str, nargs='+',
        help=('package name or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, and :experimental:)'))
    parser_download.add_argument(
        '--has-file', action='append', default=[], metavar='FILENAME', dest='has_files',
        help=('only include packages that have this file in their metadata directory '
              '(examples: SPKG.rst, spkg-configure.m4, distros/debian.txt, spkg-install|spkg-install.in)'))
    parser_download.add_argument(
        '--no-file', action='append', default=[], metavar='FILENAME', dest='no_files',
        help=('only include packages that do not have this file in their metadata directory '
              '(examples: huge, patches, huge|has_nonfree_dependencies)'))
    parser_download.add_argument(
        '--exclude', nargs='*', action='append', default=[], metavar='PACKAGE_NAME',
        help='exclude package from list')
    parser_download.add_argument(
        '--allow-upstream', action="store_true",
        help='whether to fall back to downloading from the upstream URL')
    parser_download.add_argument(
        '--on-error', choices=['stop', 'warn'], default='stop',
        help='what to do if the tarball cannot be downloaded')
    parser_download.add_argument(
        '--no-check-certificate', action='store_true',
        help='do not check SSL certificates for https connections')
    parser_download.add_argument(
        '--tags', nargs='*', default=[],
        help='ordered list of platform compatibility tags to try (for platform-dependent wheels only)')

    parser_upload = subparsers.add_parser(
        'upload', epilog=epilog_upload,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='upload tarball to Sage mirrors')
    parser_upload.add_argument(
        'package_name', type=str, help='package name or :type:')

    parser_fix_checksum = subparsers.add_parser(
        'fix-checksum', epilog=epilog_fix_checksum,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='fix the checksum of normal packages')
    parser_fix_checksum.add_argument(
        'package_class', metavar='[PACKAGE_NAME|pkg:pypi/DISTRIBUTION-NAME|:PACKAGE_TYPE:]',
        type=str, default=[':all:'], nargs='*',
        help=('package name, pkg:pypi/ followed by a distribution name, '
              'or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, and :experimental:; default: :all:)'))

    parser_create = subparsers.add_parser(
        'create', epilog=epilog_create,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='create or overwrite a package')
    parser_create.add_argument(
        'package_name', default=None, type=str,
        help='package name')
    parser_create.add_argument(
        '--source', type=str, default=None, help='package source (one of normal, wheel, script, pip); default depends on provided arguments')
    parser_create.add_argument(
        '--version', type=str, default=None, help='package version')
    parser_create.add_argument(
        '--tarball', type=str, default=None, help='tarball filename pattern, e.g. Foo-VERSION.tar.bz2')
    parser_create.add_argument(
        '--type', type=str, default=None, help='package type')
    parser_create.add_argument(
        '--url', type=str, default=None, help='download URL pattern, e.g. http://example.org/Foo-VERSION.tar.bz2')
    parser_create.add_argument(
        '--description', type=str, default=None, help='short description of the package (for SPKG.rst)')
    parser_create.add_argument(
        '--license', type=str, default=None, help='license of the package (for SPKG.rst)')
    parser_create.add_argument(
        '--upstream-contact', type=str, default=None, help='upstream contact (for SPKG.rst)')
    parser_create.add_argument(
        '--pypi', action="store_true",
        help='create a package for a Python package available on PyPI')

    parser_clean = subparsers.add_parser(
        'clean', epilog=epilog_clean,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='remove outdated source tarballs from the upstream/ directory')

    parser_metrics = subparsers.add_parser(
        'metrics', epilog=epilog_metrics,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='print metrics of given packages')
    parser_metrics.add_argument(
        'package_class', metavar='[PACKAGE_NAME|pkg:pypi/DISTRIBUTION-NAME|:PACKAGE_TYPE:]',
        type=str, nargs='*', default=[':all:'],
        help=('package name, pkg:pypi/ followed by a distribution name, '
              'or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, and :experimental:; default: :all:)'))

    return parser


def run():
    parser = make_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args(sys.argv[1:])
    if args.log is not None:
        level = getattr(logging, args.log.upper())
        log.setLevel(level=level)
    log.debug('Commandline arguments: %s', args)
    app = Application()
    if args.subcommand == 'config':
        app.config()
    elif args.subcommand == 'list':
        if args.package_class == [':all-or-nothing:']:
            if args.include_dependencies or args.exclude_dependencies:
                args.package_class = []
            else:
                args.package_class = [':all:']
        app.list_cls(*args.package_class,
                     has_files=args.has_files, no_files=args.no_files,
                     exclude=args.exclude,
                     include_dependencies=args.include_dependencies,
                     exclude_dependencies=args.exclude_dependencies)
    elif args.subcommand == 'properties':
        app.properties(*args.package_class, format=args.format)
    elif args.subcommand == 'dependencies':
        types = []
        if args.order_only:
            types.append('order_only')
        if args.optional:
            types.append('optional')
        if args.runtime:
            types.append('runtime')
        if args.check:
            types.append('check')
        if not types:
            types = None
        app.dependencies(*args.package_class, types=types, format=args.format)
    elif args.subcommand == 'name':
        app.name(args.tarball_filename)
    elif args.subcommand == 'tarball':
        app.tarball(args.package_name)
    elif args.subcommand == 'apropos':
        app.apropos(args.incorrect_name)
    elif args.subcommand == 'update':
        app.update(args.package_name, args.new_version, url=args.url, commit=args.commit)
    elif args.subcommand == 'update-latest':
        app.update_latest_cls(args.package_name, commit=args.commit)
    elif args.subcommand == 'download':
        if args.no_check_certificate:
            try:
                import ssl
                ssl._create_default_https_context = ssl._create_unverified_context
            except ImportError:
                pass
        app.download_cls(*args.package_class,
                         has_files=args.has_files, no_files=args.no_files,
                         exclude=args.exclude,
                         allow_upstream=args.allow_upstream,
                         on_error=args.on_error,
                         tags=args.tags)
    elif args.subcommand == 'create':
        app.create(args.package_name, args.version, args.tarball, args.type, args.url,
                   args.description, args.license, args.upstream_contact,
                   pypi=args.pypi, source=args.source)
    elif args.subcommand == 'upload':
        app.upload_cls(args.package_name)
    elif args.subcommand == 'fix-checksum':
        app.fix_checksum_cls(*args.package_class)
    elif args.subcommand == 'clean':
        app.clean()
    elif args.subcommand == 'metrics':
        app.metrics_cls(*args.package_class)
    else:
        raise RuntimeError('unknown subcommand: {0}'.format(args))


if __name__ == '__main__':
    run()
