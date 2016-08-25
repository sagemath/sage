# -*- coding: utf-8 -*-
"""
View for the Download Commandline UI

This module handles the main "sage-download-file" commandline utility.
"""


#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import sys
import logging
log = logging.getLogger()

# Note that argparse is not part of Python 2.6, so we bundle it
try:
    import argparse
except ImportError:
    from sage_bootstrap.compat import argparse

from sage_bootstrap.download.app import Application
from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.util import is_url


description = \
"""
Download files from a given URL or from the Sage mirror network.
"""



def make_parser():
    """
    The commandline argument parser for sage-download-file
    """
    parser = argparse.ArgumentParser(
        description=description,
    )
    parser.add_argument('--log', dest='log', default=None,
                        help='one of [DEBUG, INFO, ERROR, WARNING, CRITICAL]')

    parser.add_argument(
        '--print-fastest-mirror', action='store_true', 
        help='Print out the fastest mirror. All other arguments are ignored in that case.')

    parser.add_argument(
        '--quiet', action='store_true',
        help='Hide progress bar')

    parser.add_argument(
        '--timeout', type=float, default=None,
        help='Timeout for network operations')

    parser.add_argument(
        'url_or_tarball', type=str, nargs='?', default=None,
        help="""A http:// url or a tarball filename. In the latter case, the
        tarball is downloaded from the mirror network and its checksum
        is verified.""")

    parser.add_argument(
        'destination', type=str, nargs='?', default=None,
        help="""Where to write the file. If the destination is not specified, a url
        will be downloaded and the content written to stdout and a
        tarball will be saved under {SAGE_DISTFILES}""".format(SAGE_DISTFILES=SAGE_DISTFILES))
    
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
    app = Application(timeout=args.timeout, quiet=args.quiet)
    if (not args.print_fastest_mirror) and (args.url_or_tarball is None):
        parser.print_help()
        print('')
        print('error: either --print-fastest-mirror or url_or_tarball is required')
        sys.exit(2)
    if args.print_fastest_mirror:
        app.print_fastest_mirror()
    elif is_url(args.url_or_tarball):
        app.download_url(args.url_or_tarball, args.destination)
    else:
        app.download_tarball(args.url_or_tarball, args.destination)


def format_error(message):
    stars = '*' * 72 + '\n'
    sys.stderr.write(stars)
    try:
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.stderr.write(stars)
    except:
        pass
    sys.stderr.write(message)
    sys.stderr.write(stars)

                
def run_safe():
    try:
        run()
    except StandardError as error:
        try:
            format_error(error)
        finally:
            sys.exit(1)

                
        
if __name__ == '__main__':
    run_safe()
