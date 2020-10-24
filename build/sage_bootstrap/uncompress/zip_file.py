"""
Zip file support
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

import zipfile
from sage_bootstrap.uncompress.filter_os_files import filter_os_files


class SageZipFile(zipfile.ZipFile):
    """
    Wrapper for zipfile.ZipFile to provide better API fidelity with
    SageTarFile insofar as it's used by this script.
    """

    @classmethod
    def can_read(cls, filename):
        """
        Given an archive filename, returns True if this class can read and
        process the archive format of that file.
        """

        return zipfile.is_zipfile(filename)

    @property
    def names(self):
        """
        List of filenames in the archive.

        Filters out names of OS-related files that shouldn't be in the
        archive (.DS_Store, etc.)
        """

        return filter_os_files(self.namelist())

    def extractbytes(self, member):
        """
        Return the contents of the specified archive member as bytes.

        If the member does not exist, returns None.
        """

        if member in self.namelist():
            return self.read(member)


