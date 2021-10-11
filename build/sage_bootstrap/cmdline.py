# -*- coding: utf-8 -*-
"""
View for the Commandline UI

This module handles the main "sage-package" commandline utility, which
is also exposed as "sage --package".
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
Print a list of packages known to Sage

EXAMPLE:

    $ sage --package list | sort
    4ti2
    arb
    autotools
    [...]
    zn_poly

    $ sage --package list :standard: | sort
    arb
    backports_ssl_match_hostname
    [...]
    zn_poly
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
        help='Print the configuration')

    parser_list = subparsers.add_parser(
        'list', epilog=epilog_list,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Print a list of packages known to Sage')
    parser_list.add_argument(
        'package_class', metavar='[package_name|:package_type:]',
        type=str, default=[':all:'], nargs='*',
        help=('package name or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, :experimental:, and :huge:); '
              'default: :all:'))
    parser_list.add_argument(
        '--has-file', action='append', default=[], metavar='FILENAME', dest='has_files',
        help=('only include packages that have this file in their metadata directory'
              '(examples: SPKG.rst, spkg-configure.m4, distros/debian.txt)'))
    parser_name = subparsers.add_parser(
        'name', epilog=epilog_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Find the package name given a tarball filename')
    parser_name.add_argument('tarball_filename', type=str, help='Tarball filename')

    parser_tarball = subparsers.add_parser(
        'tarball', epilog=epilog_tarball,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Find the tarball filename given a package name')
    parser_tarball.add_argument('package_name', type=str, help='Package name')
    
    parser_apropos = subparsers.add_parser(
        'apropos', epilog=epilog_apropos,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Find up to 5 package names that are close to the given name')
    parser_apropos.add_argument(
        'incorrect_name', type=str, 
        help='Fuzzy name to search for')

    parser_update = subparsers.add_parser(
        'update', epilog=epilog_update,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Update a package. This modifies the Sage sources.')
    parser_update.add_argument(
        'package_name', type=str, help='Package name')
    parser_update.add_argument(
        'new_version', type=str, help='New version')
    parser_update.add_argument(
        '--url', type=str, default=None, help='Download URL')

    parser_update_latest = subparsers.add_parser(
        'update-latest', epilog=epilog_update_latest,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Update a package to the latest version. This modifies the Sage sources.')
    parser_update_latest.add_argument(
        'package_name', type=str, help='Package name (:all: for all packages)')

    parser_download = subparsers.add_parser(
        'download', epilog=epilog_download,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Download tarball')
    parser_download.add_argument(
        'package_name', type=str, help='Package name or :type:')
    parser_download.add_argument(
        '--allow-upstream', action="store_true",
        help='Whether to fall back to downloading from the upstream URL')
    parser_download.add_argument(
        '--on-error', choices=['stop', 'warn'], default='stop',
        help='What to do if the tarball cannot be downloaded')

    parser_upload = subparsers.add_parser(
        'upload', epilog=epilog_upload,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Upload tarball to Sage mirrors')
    parser_upload.add_argument(
        'package_name', type=str, help='Package name or :type:')
    
    parser_fix_checksum = subparsers.add_parser(
        'fix-checksum', epilog=epilog_fix_checksum,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Fix the checksum of normal packages.')
    parser_fix_checksum.add_argument(
        'package_class', metavar='[package_name|:package_type:]',
        type=str, default=[':all:'], nargs='*',
        help=('package name or designator for all packages of a given type '
              '(one of :all:, :standard:, :optional:, :experimental:, and :huge:); '
              'default: :all:'))

    parser_create = subparsers.add_parser(
        'create', epilog=epilog_create,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Create or overwrite package.')
    parser_create.add_argument(
        'package_name', default=None, type=str,
        help='Package name.')
    parser_create.add_argument(
        '--source', type=str, default='normal', help='Package source (one of normal, script, pip)')
    parser_create.add_argument(
        '--version', type=str, default=None, help='Package version')
    parser_create.add_argument(
        '--tarball', type=str, default=None, help='Tarball filename pattern, e.g. Foo-VERSION.tar.bz2')
    parser_create.add_argument(
        '--type', type=str, default=None, help='Package type')
    parser_create.add_argument(
        '--url', type=str, default=None, help='Download URL pattern, e.g. http://example.org/Foo-VERSION.tar.bz2')
    parser_create.add_argument(
        '--description', type=str, default=None, help='Short description of the package (for SPKG.rst)')
    parser_create.add_argument(
        '--license', type=str, default=None, help='License of the package (for SPKG.rst)')
    parser_create.add_argument(
        '--upstream-contact', type=str, default=None, help='Upstream contact (for SPKG.rst)')
    parser_create.add_argument(
        '--pypi', action="store_true",
        help='Create a package for a Python package available on PyPI')

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
        app.list_cls(*args.package_class, has_files=args.has_files)
    elif args.subcommand == 'name':
        app.name(args.tarball_filename)
    elif args.subcommand == 'tarball':
        app.tarball(args.package_name)
    elif args.subcommand == 'apropos':
        app.apropos(args.incorrect_name)
    elif args.subcommand == 'update':
        app.update(args.package_name, args.new_version, url=args.url)
    elif args.subcommand == 'update-latest':
        app.update_latest_cls(args.package_name)
    elif args.subcommand == 'download':
        app.download_cls(args.package_name,
                         allow_upstream=args.allow_upstream,
                         on_error=args.on_error)
    elif args.subcommand == 'create':
        app.create(args.package_name, args.version, args.tarball, args.type, args.url,
                   args.description, args.license, args.upstream_contact,
                   pypi=args.pypi, source=args.source)
    elif args.subcommand == 'upload':
        app.upload_cls(args.package_name)
    elif args.subcommand == 'fix-checksum':
        app.fix_checksum_cls(*args.package_class)
    else:
        raise RuntimeError('unknown subcommand: {0}'.format(args))

        
if __name__ == '__main__':
    run()
