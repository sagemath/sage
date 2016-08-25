# -*- coding: utf-8 -*-
"""
Third-Party Tarballs
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.download import Download, MirrorList
from sage_bootstrap.package import Package


class ChecksumError(Exception):
    """
    Exception raised when the checksum of the tarball does not match
    """
    pass


class FileNotMirroredError(Exception):
    """
    Exception raised when the tarball cannot be downloaded from the mirrors
    """
    pass


class Tarball(object):
    
    def __init__(self, tarball_name, package=None):
        """
        A (third-party downloadable) tarball

        Note that the tarball might also be a different kind of
        archive format that is supported, it does not necessarily have
        to be tar.

        INPUT:

        - ``name`` - string. The full filename (``foo-1.3.tar.bz2``)
          of a tarball on the Sage mirror network.
        """
        self.__filename = tarball_name
        if package is None:
            self.__package = None
            for pkg in Package.all():
                if pkg.tarball_filename == tarball_name:
                    self.__package = pkg
            if self.package is None:
                error = 'tarball {0} is not referenced by any Sage package'.format(tarball_name)
                log.error(error)
                raise ValueError(error)
        else:
            self.__package = package
            if package.tarball_filename != tarball_name:
                error = 'tarball {0} is not referenced by the {1} package'.format(tarball_name, package.name)
                log.error(error)
                raise ValueError(error)

    def __repr__(self):
        return 'Tarball {0}'.format(self.filename)
            
    @property
    def filename(self):
        """
        Return the tarball filename

        OUTPUT:

        String. The full filename (``foo-1.3.tar.bz2``) of the
        tarball.
        """
        return self.__filename
        
    @property
    def package(self):
        """
        Return the package that the tarball belongs to

        OUTPUT:

        Instance of :class:`sage_bootstrap.package.Package`
        """
        return self.__package
        
    @property
    def upstream_fqn(self):
        """
        The fully-qualified (including directory) file name in the upstream directory.
        """
        return os.path.join(SAGE_DISTFILES, self.filename)

    def __eq__(self, other):
        return self.filename == other.filename
        
    def _compute_hash(self, algorithm):
        with open(self.upstream_fqn, 'rb') as f:
            while True:
                buf = f.read(0x100000)
                if not buf:
                    break
                algorithm.update(buf)
        return algorithm.hexdigest()

    def _compute_sha1(self):
        import hashlib
        return self._compute_hash(hashlib.sha1())

    def _compute_md5(self):
        import hashlib
        return self._compute_hash(hashlib.md5())
    
    def _compute_cksum(self):
        from sage_bootstrap.cksum import CksumAlgorithm
        return self._compute_hash(CksumAlgorithm())
    
    def checksum_verifies(self):
        """
        Test whether the checksum of the downloaded file is correct.
        """
        sha1 = self._compute_sha1()
        return sha1 == self.package.sha1

    def download(self):
        """
        Download the tarball to the upstream directory.
        """
        destination = os.path.join(SAGE_DISTFILES, self.filename)
        if os.path.isfile(destination):
            if self.checksum_verifies():
                log.info('Using cached file {destination}'.format(destination=destination))
                return
            else:
                # Garbage in the upstream directory? Delete and re-download
                log.info('Invalid checksum for cached file {destination}, deleting'
                         .format(destination=destination))
                os.remove(destination)
        successful_download = False
        log.info('Attempting to download package {0} from mirrors'.format(self.filename))
        for mirror in MirrorList():
            url = mirror + '/'.join(['spkg', 'upstream', self.package.name, self.filename])
            log.info(url)
            try:
                Download(url, self.upstream_fqn).run()
                successful_download = True
                break
            except IOError:
                log.debug('File not on mirror')
        if not successful_download:
            raise FileNotMirroredError('tarball does not exist on mirror network')
        if not self.checksum_verifies():
            raise ChecksumError('checksum does not match')

    def save_as(self, destination):
        """
        Save the tarball as a new file
        """
        import shutil
        shutil.copy(self.upstream_fqn, destination)

