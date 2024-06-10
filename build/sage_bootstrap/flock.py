#!/usr/bin/env sage-bootstrap-python
# vim: set filetype=python:
"""
This script runs the given command under a file lock (similar to the flock
command on some systems).
"""

# This is originally motivated by pip, but has since been generalized.  We
# should avoid running pip while uninstalling a package because that is prone
# to race conditions. This script runs pip under a lock.  For details, see
# https://github.com/sagemath/sage/issues/21672

import fcntl
import os
import pipes
import sys
import argparse


class FileType(argparse.FileType):
    """
    Version of argparse.FileType with the option to ensure that the full path
    to the file exists.
    """

    def __init__(self, mode='r', makedirs=False):
        # Note, the base class __init__ takes other arguments too depending on
        # the Python version but we don't care about them for this purpose
        super(FileType, self).__init__(mode=mode)
        self._makedirs = makedirs

    def __call__(self, string):
        if self._makedirs and string != '-':
            dirname = os.path.dirname(string)
            try:
                os.makedirs(dirname)
            except OSError as exc:
                if not os.path.isdir(dirname):
                    raise argparse.ArgumentTypeError(
                        "can't create '{0}': {1}".format(dirname, exc))

        return super(FileType, self).__call__(string)


class IntOrFileType(FileType):
    """
    Like FileType but also accepts an int (e.g. for a file descriptor).
    """

    def __call__(self, string):
        try:
            return int(string)
        except ValueError:
            return super(IntOrFileType, self).__call__(string)


def run(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--shared', action='store_true',
                       help='create a shared lock')

    # Note: A exclusive lock is created by default if no other flags are given,
    # but supplying the --exclusive flag explicitly may help clarity
    group.add_argument('-x', '--exclusive', action='store_true',
                       help='create an exclusive lock (the default)')
    group.add_argument('-u', '--unlock', action='store_true',
                       help='remove an existing lock')
    parser.add_argument('lock', metavar='LOCK',
                        type=IntOrFileType('w+', makedirs=True),
                        help='filename of the lock an integer file descriptor')
    parser.add_argument('command', metavar='COMMAND', nargs=argparse.REMAINDER,
                        help='command to run with the lock including any '
                             'arguments to that command')

    args = parser.parse_args(argv)

    if args.shared:
        locktype = fcntl.LOCK_SH
    elif args.unlock:
        locktype = fcntl.LOCK_UN
    else:
        locktype = fcntl.LOCK_EX

    lock = args.lock
    command = args.command

    if isinstance(lock, int) and command:
        parser.error('sage-flock does not accept a command when passed '
                     'a file descriptor number')

    # First try a non-blocking lock such that we can give an informative
    # message while the user is waiting.
    try:
        fcntl.flock(lock, locktype | fcntl.LOCK_NB)
    except IOError as exc:
        if locktype == fcntl.LOCK_SH:
            kind = "shared"
        elif locktype == fcntl.LOCK_UN:
            # This shouldn't happen
            sys.stderr.write(
                "Unexpected error trying to unlock fd: {0}\n".format(exc))
            return 1
        else:
            kind = "exclusive"

        sys.stderr.write("Waiting for {0} lock to run {1} ... ".format(
            kind, ' '.join(pipes.quote(arg) for arg in command)))
        fcntl.flock(lock, locktype)
        sys.stderr.write("ok\n")

    if not (args.unlock or isinstance(lock, int)):
        os.execvp(command[0], command)


if __name__ == '__main__':
    sys.exit(run())
