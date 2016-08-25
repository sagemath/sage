# -*- coding: utf-8 -*-
"""
Package Creator
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

import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_ROOT



class PackageCreator(object):

    def __init__(self, package_name):
        self.package_name = package_name
        self.path = os.path.join(SAGE_ROOT, 'build', 'pkgs', package_name)
        try:
            os.mkdir(self.path)
        except OSError:
            pass

    def set_version(self, version):
        """
        Write the version to ``package-version.txt``
        """
        with open(os.path.join(self.path, 'package-version.txt'), 'w+') as f:
            f.write(version)
            f.write('\n')
            
    def set_type(self, pkg_type):
        """
        Write the package type to ``type``
        """
        with open(os.path.join(self.path, 'type'), 'w+') as f:
            f.write(pkg_type)
            f.write('\n')
            
    def set_tarball(self, tarball):
        """
        Write the tarball name pattern to ``checksums.ini``
        """
        with open(os.path.join(self.path, 'checksums.ini'), 'w+') as f:
            f.write('tarball={0}'.format(tarball))
            f.write('\n')
            
        
