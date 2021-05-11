"""
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

from sage_bootstrap.uncompress.tar_file import SageTarFile, SageTarXZFile
from sage_bootstrap.uncompress.zip_file import SageZipFile
from sage_bootstrap.util import retry

ARCHIVE_TYPES = [SageTarFile, SageZipFile, SageTarXZFile]



def open_archive(filename):
    """
    Automatically detect archive type
    """
    for cls in ARCHIVE_TYPES:
        if cls.can_read(filename):
            break
    else:
        raise ValueError

    # For now ZipFile and TarFile both have default open modes that are
    # acceptable
    return cls(filename)


def unpack_archive(archive, dirname=None):
    """
    Unpack archive
    """
    top_level = None

    if dirname:
        top_levels = set()
        for member in archive.names:
            # Zip and tar files all use forward slashes as separators
            # internally
            top_levels.add(member.split('/', 1)[0])

        if len(top_levels) == 1:
            top_level = top_levels.pop()
        else:
            os.makedirs(dirname)

    prev_cwd = os.getcwd()

    if dirname and not top_level:
        # We want to extract content into dirname, but there is not
        # a single top-level directory for the tarball, so we cd into
        # the extraction target first
        os.chdir(dirname)

    try:
        archive.extractall(members=archive.names)
        if dirname and top_level:
            # On Windows os.rename can fail unexpectedly with a permission
            # error if a virus scanner or other background process is
            # inspecting the newly extracted files
            rename = lambda: os.rename(top_level, dirname)
            retry(rename, OSError, tries=len(archive.names))

            # Apply typical umask to the top-level directory in case it wasn't
            # already; see https://trac.sagemath.org/ticket/24567
            # and later https://trac.sagemath.org/ticket/28596
            os.chmod(dirname, os.stat(dirname).st_mode & ~0o022)
    finally:
        os.chdir(prev_cwd)
