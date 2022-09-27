"""
Command-line script for checking an installed SPKG in an installation tree ($SAGE_LOCAL, $SAGE_VENV).
"""

# ****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.m.bray@gmail.com>
#                     2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function

import glob
import json
import os
import shutil
import subprocess
import sys
import argparse
import warnings

from .env import SAGE_ROOT

pth = os.path
PKGS = pth.join(SAGE_ROOT, 'build', 'pkgs')
"""Directory where all spkg sources are found."""


def check_lib_auditwheel(f, verbose=False):
    from auditwheel.lddtree import lddtree
    for lib, info in lddtree(f)["libs"].items():
        if verbose:
            print('- {0}: {1}'.format(lib, info["realpath"]), file=sys.stderr)
        if info["realpath"] is None:
            raise RuntimeError('Shared library {0} needed by {1} is not found'
                               .format(lib, f))


def installcheck(spkg_name, sage_local, verbose=False):
    """
    Given a package name and path to an installation tree (SAGE_LOCAL or SAGE_VENV),
    check the installation of the package in that tree.
    """

    # The default path to this directory; however its value should be read
    # from the environment if possible
    spkg_inst = pth.join(sage_local, 'var', 'lib', 'sage', 'installed')

    # Find all stamp files for the package; there should be only one, but if
    # there is somehow more than one we'll work with the most recent one.
    pattern = pth.join(spkg_inst, '{0}-*'.format(spkg_name))
    stamp_files = sorted(glob.glob(pattern), key=pth.getmtime)

    if stamp_files:
        stamp_file = stamp_files[-1]
    else:
        stamp_file = None

    spkg_meta = {}
    if stamp_file:
        try:
            with open(stamp_file) as f:
                spkg_meta = json.load(f)
        except (OSError, ValueError):
            pass

    if 'files' not in spkg_meta:
        if stamp_file:
            print("Old-style or corrupt stamp file '{0}'"
                  .format(stamp_file), file=sys.stderr)
        else:
            print("Package '{0}' is currently not installed in '{1}'"
                  .format(spkg_name, sage_local), file=sys.stderr)
    else:
        files = spkg_meta['files']

        for f in files:
            f = os.path.join(sage_local, f)
            if f.endswith(('.so', '.dylib')):
                if verbose:
                    print("Checking shared library file '{0}'"
                          .format(f), file=sys.stderr)
                if sys.platform == 'darwin':
                    try:
                        from delocate.libsana import _tree_libs_from_libraries, _filter_system_libs
                    except ImportError:
                        warnings.warn('delocate is not available, so nothing is actually checked')
                    else:
                        _tree_libs_from_libraries([f],
                                                  lib_filt_func=_filter_system_libs,
                                                  copy_filt_func=lambda path: True)
                else:
                    try:
                        check_lib_auditwheel(f, verbose=False)
                    except ImportError:
                        warnings.warn('auditwheel is not available, so nothing is actually checked')
            elif f.endswith('-any.whl'):
                # pure Python wheel, nothing to check
                pass
            elif f.endswith('.whl'):
                if verbose:
                    print("Checking wheel file '{0}'"
                          .format(f), file=sys.stderr)
                if sys.platform == 'darwin':
                    try:
                        from delocate import wheel_libs
                    except ImportError:
                        warnings.warn('delocate is not available, so nothing is actually checked')
                    else:
                        wheel_libs(f)
                else:
                    try:
                        from delocate.tmpdirs import TemporaryDirectory
                        from delocate.tools import zip2dir
                    except ImportError:
                        warnings.warn('delocate is not available, so nothing is actually checked')
                    else:
                        try:
                            with TemporaryDirectory() as tmpdir:
                                zip2dir(f, tmpdir)
                                for dirpath, dirnames, basenames in os.walk(tmpdir):
                                    for base in basenames:
                                        if base.endswith('.so'):
                                            depending_path = os.path.realpath(os.path.join(dirpath, base))
                                            check_lib_auditwheel(depending_path, verbose=False)
                        except ImportError:
                            warnings.warn('auditwheel is not available, so nothing is actually checked')


def dir_type(path):
    """
    A custom argument 'type' for directory paths.
    """

    if path and not pth.isdir(path):
        raise argparse.ArgumentTypeError(
            "'{0}' is not a directory".format(path))

    return path


def spkg_type(pkg):
    """
    A custom argument 'type' for spkgs--checks whether the given package name
    is a known spkg.
    """
    pkgbase = pth.join(PKGS, pkg)

    if not pth.isdir(pkgbase):
        raise argparse.ArgumentTypeError(
            "'{0}' is not an spkg listed in '{1}'".format(pkg, PKGS))

    return pkg


def make_parser():
    """Returns the command-line argument parser for sage-spkg-installcheck."""

    doc_lines = __doc__.strip().splitlines()

    parser = argparse.ArgumentParser(
        description=doc_lines[0],
        epilog='\n'.join(doc_lines[1:]).strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('spkg', type=spkg_type, help='the spkg to check')
    parser.add_argument('sage_local', type=dir_type, nargs='?',
                        default=os.environ.get('SAGE_LOCAL'),
                        help='the path of the installation tree (default: the $SAGE_LOCAL '
                             'environment variable if set)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output showing all files removed')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    return parser


def run(argv=None):
    parser = make_parser()

    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    if args.sage_local is None:
        print('Error: An installation tree must be specified either at the command '
              'line or in the $SAGE_LOCAL environment variable',
              file=sys.stderr)
        sys.exit(1)

    try:
        installcheck(args.spkg, args.sage_local,
                  verbose=args.verbose)
    except Exception as exc:
        print("Error during installcheck of '{0}': {1}".format(
            args.spkg, exc), file=sys.stderr)

        if args.debug:
            raise

        sys.exit(1)


if __name__ == '__main__':
    run()
