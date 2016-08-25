# -*- coding: utf-8 -*-
"""
Sage Packages
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

import re
import os

import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_ROOT, SAGE_DISTFILES



class Package(object):

    def __init__(self, package_name):
        """
        Sage Package

        A package is defined by a subdirectory of
        ``SAGE_ROOT/build/pkgs/``. The name of the package is the name
        of the subdirectory; The metadata of the package is contained
        in various files in the package directory. This class provides
        an abstraction to the metadata, you should never need to
        access the package directory directly.

        INPUT: 

        -- ``package_name`` -- string. Name of the package. The Sage
           convention is that all package names are lower case.
        """
        self.__name = package_name
        self.__tarball = None
        self._init_checksum()
        self._init_version()

    def __repr__(self):
        return 'Package {0}'.format(self.name)
            
    @property
    def name(self):
        """
        Return the package name

         A package is defined by a subdirectory of
        ``SAGE_ROOT/build/pkgs/``. The name of the package is the name
        of the subdirectory.

        OUTPUT:

        String.
        """
        return self.__name

    @property
    def md5(self):
        """
        Return the MD5 checksum
        
        Do not use, this is ancient! Use :meth:`sha1` instead.

        OUTPUT:

        String.
        """
        return self.__md5

    @property
    def sha1(self):
        """
        Return the SHA1 checksum

        OUTPUT:

        String.
        """
        return self.__sha1

    @property
    def cksum(self):
        """
        Return the Ck sum checksum
        
        Do not use, this is ancient! Use :meth:`sha1` instead.

        OUTPUT:

        String.
        """
        return self.__cksum

    @property
    def tarball(self):
        """
        Return the (primary) tarball

        If there are multiple tarballs (currently unsupported), this
        property returns the one that is unpacked automatically.

        OUTPUT:

        Instance of :class:`sage_bootstrap.tarball.Tarball`
        """
        if self.__tarball is None:
            from sage_bootstrap.tarball import Tarball
            self.__tarball = Tarball(self.tarball_filename, package=self)
        return self.__tarball

    @property
    def tarball_pattern(self):
        """
        Return the (primary) tarball file pattern

        If there are multiple tarballs (currently unsupported), this
        property returns the one that is unpacked automatically.

        OUTPUT:

        String. The full-qualified tarball filename, but with
        ``VERSION`` instead of the actual tarball filename.
        """
        return self.__tarball_pattern

    @property
    def tarball_filename(self):
        """
        Return the (primary) tarball filename

        If there are multiple tarballs (currently unsupported), this
        property returns the one that is unpacked automatically.

        OUTPUT:

        String. The full-qualified tarball filename.
        """
        return self.tarball_pattern.replace('VERSION', self.version)

    @property
    def version(self):
        """
        Return the version

        OUTPUT:

        String. The package version. Excludes the Sage-specific
        patchlevel.
        """
        return self.__version

    @property
    def patchlevel(self):
        """
        Return the patchlevel

        OUTPUT:

        Integer. The patchlevel of the package. Excludes the "p"
        prefix.
        """
        return self.__patchlevel

    def __eq__(self, other):
        return self.tarball == other.tarball
        
    @classmethod
    def all(cls):
        """
        Return all packages
        """
        base = os.path.join(SAGE_ROOT, 'build', 'pkgs')
        for subdir in os.listdir(base):
            path = os.path.join(base, subdir) 
            if not os.path.isfile(os.path.join(path, "checksums.ini")):
                continue
            yield cls(subdir)

    @property
    def path(self):
        """
        Return the package directory
        """
        return os.path.join(SAGE_ROOT, 'build', 'pkgs', self.name)
            
    def _init_checksum(self):
        """
        Load the checksums from the appropriate ``checksums.ini`` file
        """
        checksums_ini = os.path.join(self.path, 'checksums.ini')
        assignment = re.compile('(?P<var>[a-zA-Z0-9]*)=(?P<value>.*)')
        result = dict()
        with open(checksums_ini, 'rt') as f:
            for line in f.readlines():
                match = assignment.match(line)
                if match is None:
                    continue
                var, value = match.groups()
                result[var] = value
        self.__md5 = result.get('md5', None)
        self.__sha1 = result.get('sha1', None)
        self.__cksum = result.get('cksum', None)
        self.__tarball_pattern = result['tarball']
        
    VERSION_PATCHLEVEL = re.compile('(?P<version>.*)\.p(?P<patchlevel>[0-9]+)')
    
    def _init_version(self):
        with open(os.path.join(self.path, 'package-version.txt')) as f:
            package_version = f.read().strip()
        match = self.VERSION_PATCHLEVEL.match(package_version)
        if match is None:
            self.__version = package_version
            self.__patchlevel = -1
        else:
            self.__version = match.group('version')
            self.__patchlevel = int(match.group('patchlevel'))
        
        
