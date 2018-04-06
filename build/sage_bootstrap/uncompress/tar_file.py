"""
Tar file support
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
import copy
import tarfile
import stat
import subprocess
from io import BytesIO

from sage_bootstrap.uncompress.filter_os_files import filter_os_files


class SageBaseTarFile(tarfile.TarFile):
    """
    Sage as tarfile.TarFile, but applies the user's current umask to the
    permissions of all extracted files and directories.

    This mimics the default behavior of the ``tar`` utility.

    See http://trac.sagemath.org/ticket/20218#comment:16 for more background.
    """
    def __init__(self, *args, **kwargs):
        # Unfortunately the only way to get the current umask is to set it
        # and then restore it
        super(SageBaseTarFile, self).__init__(*args, **kwargs)
        self.umask = os.umask(0o777)
        os.umask(self.umask)

    @property
    def names(self):
        """
        List of filenames in the archive.

        Filters out names of OS-related files that shouldn't be in the
        archive (.DS_Store, etc.)
        """

        return filter_os_files(self.getnames())

    def chmod(self, tarinfo, target):
        tarinfo = copy.copy(tarinfo)
        tarinfo.mode &= ~self.umask
        tarinfo.mode &= ~(stat.S_ISUID | stat.S_ISGID)
        return super(SageBaseTarFile, self).chmod(tarinfo, target)

    def extractall(self, path='.', members=None):
        """
        Same as tarfile.TarFile.extractall but allows filenames for
        the members argument (like zipfile.ZipFile).
        """
        if members:
            name_to_member = dict([member.name, member] for member in self.getmembers())
            members = [m if isinstance(m, tarfile.TarInfo)
                       else name_to_member[m]
                       for m in members]
        return super(SageBaseTarFile, self).extractall(path=path, members=members)

    def extractbytes(self, member):
        """
        Return the contents of the specified archive member as bytes.

        If the member does not exist, returns None.
        """

        if member in self.getnames():
            reader = self.extractfile(member)
            return reader.read()


class SageTarFile(SageBaseTarFile):
    """
    A wrapper around SageBaseTarFile such that SageTarFile(filename) is
    essentially equivalent to TarFile.open(filename) which is more
    flexible than the basic TarFile.__init__
    """
    def __new__(cls, filename):
        return SageBaseTarFile.open(filename)

    @staticmethod
    def can_read(filename):
        """
        Given an archive filename, returns True if this class can read and
        process the archive format of that file.
        """
        return tarfile.is_tarfile(filename)


class SageTarXZFile(SageBaseTarFile):
    """
    A ``.tar.xz`` file which is uncompressed in memory.
    """
    def __new__(cls, filename):
        # Read uncompressed data through a pipe
        proc = subprocess.Popen(["xz", "-d", "-c", filename], stdout=subprocess.PIPE)
        data, _ = proc.communicate()
        return SageBaseTarFile(mode="r", fileobj=BytesIO(data))

    @staticmethod
    def can_read(filename):
        """
        Given an archive filename, returns True if this class can read and
        process the archive format of that file.
        """
        devnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(["xz", "-l", filename], stdout=devnull, stderr=devnull)
        except Exception:
            return False
        return True
