"""
Commandline handling for sage-uncompress-spkg
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


from __future__ import print_function

import os
import sys
import argparse

from sage_bootstrap.uncompress.action import (
    open_archive, unpack_archive
)
                                              

    
def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='dir', metavar='DIR',
                        help='directory to extract archive contents into')
    parser.add_argument('pkg', nargs=1, metavar='PKG',
                        help='the archive to extract')
    parser.add_argument('file', nargs='?', metavar='FILE',
                        help='(deprecated) print the contents of the given '
                             'archive member to stdout')
    return parser


def run():
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    filename = args.pkg[0]
    dirname = args.dir

    try:
        archive = open_archive(filename)
    except ValueError:
        print('Error: Unknown file type: {}'.format(filename),
              file=sys.stderr)
        return 1
        
    if args.file:
        contents = archive.extractbytes(args.file)
        if contents:
            print(contents, end='')
            return 0
        else:
            return 1

    if dirname and os.path.exists(dirname):
        print('Error: Directory {} already exists'.format(dirname),
              file=sys.stderr)
        return 1

    unpack_archive(archive, dirname)
    return 0


if __name__ == '__main__':
    sys.exit(run())
