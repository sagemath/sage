"""
Filtering out OS-specific files
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


def filter_os_files(filenames):
    """
    Given a list of filenames, returns a filtered list with OS-specific
    special files removed.

    Currently removes OSX .DS_Store files and AppleDouble format ._ files.
    """

    files_set = set(filenames)

    def is_os_file(path):
        dirname, name = os.path.split(path)

        if name == '.DS_Store':
            return True

        if name.startswith('._'):
            name = os.path.join(dirname, name[2:])
            # These files store extended attributes on OSX
            # In principle this could be a false positive but it's
            # unlikely, and to be really sure we'd have to extract the file
            # (or at least the first four bytes to check for the magic number
            # documented in
            # http://kaiser-edv.de/documents/AppleSingle_AppleDouble.pdf)
            if name in files_set or os.path.normpath(name) in files_set:
                return True

        return False

    filenames = filter(lambda f: not is_os_file(f), filenames)

    if sys.version_info[0] == 2:
        return filenames
    else:
        # Make sure to return a list on Python >= 3
        return list(filenames)
