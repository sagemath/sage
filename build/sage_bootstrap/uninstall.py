"""
Command-line script for uninstalling an existing SPKG from an installation tree ($SAGE_LOCAL, $SAGE_VENV).

This performs two types of uninstallation:

    1) Old-style uninstallation: This is close to what existed before this
       script, where *some* packages had uninstall steps (mostly consisting of
       some broad `rm -rf` commands) that were run before installing new
       versions of the packages.  This convention was applied inconsistently,
       but for those packages that did have old-style uninstall steps, those
       steps should be in a script called `spkg-legacy-uninstall` under the
       spkg directory (build/pkgs/<pkg_name>).  If this script is found, it is
       run for backwards-compatibility support.

    2) New-style uninstallation: More recently installed packages that were
       installed with staged installation have a record of all files installed
       by that package.  That file is stored in the directory
       $SAGE_LOCAL/var/lib/sage/installed or $SAGE_VENV/var/lib/sage/installed
       and is created when the spkg is installed.
       This is a JSON file containing some meta-data about
       the package, including the list of all files associated with the
       package.  This script removes all these files, including the record
       file.  Any directories that are empty after files are removed from them
       are also removed.
"""
# ****************************************************************************
#       Copyright (C) 2017-2018 Erik M. Bray <erik.m.bray@gmail.com>
#                     2019      Jeroen Demeyer
#                     2021-2022 Matthias Koeppe
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

from .env import SAGE_ROOT

pth = os.path
PKGS = pth.join(SAGE_ROOT, 'build', 'pkgs')
"""Directory where all spkg sources are found."""


def uninstall(spkg_name, sage_local, keep_files=False,
              verbose=False, log_directory=None):
    """
    Given a package name and path to an installation tree (SAGE_LOCAL or SAGE_VENV),
    uninstall that package from that tree if it is currently installed.
    """

    # The path to the installation records.
    # See SPKG_INST_RELDIR in SAGE_ROOT/build/make/Makefile.in
    spkg_inst = pth.join(sage_local, 'var', 'lib', 'sage', 'installed')

    # Find all stamp files for the package; there should be only one, but if
    # there is somehow more than one we'll work with the most recent and delete
    # the rest
    pattern = pth.join(spkg_inst, '{0}-*'.format(spkg_name))
    stamp_files = sorted(glob.glob(pattern), key=pth.getmtime)

    if keep_files:
        print("Removing stamp file but keeping package files")
        remove_stamp_files(stamp_files)
        return

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
            print("Old-style or corrupt stamp file {0}"
                  .format(stamp_file), file=sys.stderr)
        else:
            print("No stamp file for package '{0}' in {1}"
                  .format(spkg_name, spkg_inst), file=sys.stderr)

        # Run legacy uninstaller even if there is no stamp file: the
        # package may be partially installed without a stamp file
        legacy_uninstall(spkg_name,
                         verbose=verbose, log_directory=log_directory)
    else:
        files = spkg_meta['files']
        if not files:
            print("Warning: No files to uninstall for "
                  "'{0}'".format(spkg_name), file=sys.stderr)

        modern_uninstall(spkg_name, sage_local, files,
                         verbose=verbose, log_directory=log_directory)

    remove_stamp_files(stamp_files, verbose=verbose)


def legacy_uninstall(spkg_name,
                     verbose=False, log_directory=None):
    """
    Run the spkg's legacy uninstall script, if one exists; otherwise do
    nothing.
    """
    spkg_dir = pth.join(PKGS, spkg_name)

    # Any errors from this, including a non-zero return code will
    # bubble up and exit the uninstaller
    run_spkg_script(spkg_name, spkg_dir, 'legacy-uninstall',
                    if_does_not_exist='log', log_directory=log_directory)


def modern_uninstall(spkg_name, sage_local, files,
                     verbose=False, log_directory=None):
    """
    Remove all listed files from the given installation tree (SAGE_LOCAL or SAGE_VENV).

    All file paths should be assumed relative to the installation tree.

    This is otherwise (currently) agnostic about what package is actually
    being uninstalled--all it cares about is removing a list of files.

    If the directory containing any of the listed files is empty after all
    files are removed then the directory is removed as well.
    """

    spkg_scripts = pth.join(sage_local, 'var', 'lib', 'sage', 'scripts')
    spkg_scripts = os.environ.get('SAGE_SPKG_SCRIPTS', spkg_scripts)
    spkg_scripts = pth.join(spkg_scripts, spkg_name)

    # Sort the given files first by the highest directory depth, then by name,
    # so that we can easily remove a directory once it's been emptied
    files.sort(key=lambda f: (-f.count(os.sep), f))

    if verbose:
        print("Uninstalling existing '{0}'".format(spkg_name))

    # Run the package's prerm script, if it exists.
    # If an error occurs here we abort the uninstallation for now.
    # This means a prerm script actually has the ability to abort an
    # uninstallation, for example, if some manual intervention is needed
    # to proceed.
    try:
        run_spkg_script(spkg_name, spkg_scripts, 'prerm',
                        log_directory=log_directory)
    except Exception as exc:
        script_path = pth.join(spkg_scripts, 'spkg-prerm')
        print("Error: The pre-uninstall script for '{0}' failed; the "
              "package will not be uninstalled, and some manual intervention "
              "may be needed to repair the package's state before "
              "uninstallation can proceed.  Check further up in this log "
              "for more details, or the pre-uninstall script itself at "
              "{1}.".format(spkg_name, script_path), file=sys.stderr)
        if isinstance(exc, subprocess.CalledProcessError):
            sys.exit(exc.returncode)
        else:
            sys.exit(1)

    # Run the package's piprm script, if it exists.
    # Since #36452 the spkg-requirements.txt file appears in the installation
    # manifest, so this step has to happen before removing the files.
    try:
        run_spkg_script(spkg_name, spkg_scripts, 'piprm',
                        log_directory=log_directory)
    except Exception:
        print("Warning: Error running the pip-uninstall script for "
              "'{0}'; uninstallation may have left behind some files".format(
              spkg_name), file=sys.stderr)

    def rmdir(dirname):
        if pth.isdir(dirname):
            if not os.listdir(dirname):
                if verbose:
                    print('rmdir "{}"'.format(dirname))
                os.rmdir(dirname)
        else:
            print("Warning: Directory {0} not found".format(
                dirname), file=sys.stderr)

    # Remove the files; if a directory is empty after removing a file
    # from it, remove the directory too.
    for filename in files:
        # Just in case: use lstrip to remove leading "/" from
        # filename. See https://github.com/sagemath/sage/issues/26013.
        filename = pth.join(sage_local, filename.lstrip(os.sep))
        dirname = pth.dirname(filename)
        if pth.lexists(filename):
            if verbose:
                print('rm "{}"'.format(filename))
            os.remove(filename)
        else:
            print("Warning: File {0} not found".format(filename),
                  file=sys.stderr)

        # Remove file's directory if it is now empty
        rmdir(dirname)

    # Run the package's postrm script, if it exists.
    # If an error occurs here print a warning, but complete the
    # uninstallation; otherwise we leave the package in a broken
    # state--looking as though it's still 'installed', but with all its
    # files removed.
    try:
        run_spkg_script(spkg_name, spkg_scripts, 'postrm',
                        log_directory=log_directory)
    except Exception:
        print("Warning: Error running the post-uninstall script for "
              "'{0}'; the package will still be uninstalled, but "
              "may have left behind some files or settings".format(
              spkg_name), file=sys.stderr)

    try:
        shutil.rmtree(spkg_scripts)
    except Exception:
        pass


def remove_stamp_files(stamp_files, verbose=False):
    # Finally, if all went well, delete all the stamp files.
    for stamp_file in stamp_files:
        print("Removing stamp file {0}".format(stamp_file))
        os.remove(stamp_file)


def run_spkg_script(spkg_name, path, script_name,
                    if_does_not_exist='ignore', log_directory=None):
    """
    Runs the specified ``spkg-<foo>`` script under the given ``path``,
    if it exists.
    """
    script_name = 'spkg-{0}'.format(script_name)
    script = pth.join(path, script_name)
    if pth.exists(script):
        if log_directory:
            log_file = pth.join(log_directory, script + '.log')
            subprocess.check_call(['sage-logger', '-p', script, log_file])
        else:
            subprocess.check_call(['sage-logger', '-P', script_name, script])
    elif if_does_not_exist == 'ignore':
        pass
    elif if_does_not_exist == 'log':
        print('No {0} script; nothing to do'.format(script_name), file=sys.stderr)
    else:
        raise ValueError('unknown if_does_not_exist value: {0}'.format(if_does_not_exist))


def dir_type(path):
    """
    A custom argument 'type' for directory paths.
    """

    if path and not pth.isdir(path):
        raise argparse.ArgumentTypeError(
            "{0} is not a directory".format(path))

    return path


def make_parser():
    """Returns the command-line argument parser for sage-spkg-uninstall."""

    doc_lines = __doc__.strip().splitlines()

    parser = argparse.ArgumentParser(
        description=doc_lines[0],
        epilog='\n'.join(doc_lines[1:]).strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('spkg', type=str, help='the spkg to uninstall')
    parser.add_argument('sage_local', type=dir_type, nargs='?',
                        default=os.environ.get('SAGE_LOCAL'),
                        help='the path of the installation tree (default: the $SAGE_LOCAL '
                             'environment variable if set)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output showing all files removed')
    parser.add_argument('-k', '--keep-files', action='store_true',
                        help="only delete the package's installation record, "
                             "but do not remove files installed by the "
                             "package")
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--log-directory', type=str,
                        help="directory where to create log files (default: none)")

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
        uninstall(args.spkg, args.sage_local, keep_files=args.keep_files,
                  verbose=args.verbose, log_directory=args.log_directory)
    except Exception as exc:
        print("Error during uninstallation of '{0}': {1}".format(
            args.spkg, exc), file=sys.stderr)

        if args.debug:
            raise

        sys.exit(1)


if __name__ == '__main__':
    run()
