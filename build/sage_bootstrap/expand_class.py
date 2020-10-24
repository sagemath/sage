# -*- coding: utf-8 -*-
"""
Utility to specify classes of packages like `:all:`
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

from sage_bootstrap.package import Package


class PackageClass(object):

    def __init__(self, package_name_or_class):
        if package_name_or_class == ':all:':
            self._init_all()
        elif package_name_or_class == ':standard:':
            self._init_standard()
        elif package_name_or_class == ':optional:':
            self._init_optional()
        elif package_name_or_class == ':experimental:':
            self._init_experimental()
        elif package_name_or_class == ':huge:':
            self._init_huge()
        else:
            if package_name_or_class.startswith(':'):
                raise ValueError('Package name cannot start with ":", got %s', package_name_or_class)
            if package_name_or_class.endswith(':'):
                raise ValueError('Package name cannot end with ":", got %s', package_name_or_class)
            self.names = [package_name_or_class]

    def _init_all(self):
        self.names = [pkg.name for pkg in Package.all()]
            
    def _init_standard(self):
        self.names = [pkg.name for pkg in Package.all() if pkg.type == 'standard']
            
    def _init_optional(self):
        self.names = [pkg.name for pkg in Package.all() if pkg.type == 'optional']
            
    def _init_experimental(self):
        self.names = [pkg.name for pkg in Package.all() if pkg.type == 'experimental']
            
    def _init_huge(self):
        self.names = [pkg.name for pkg in Package.all() if pkg.type == 'huge']
            
    def apply(self, function, *args, **kwds):
        for package_name in self.names:
            function(package_name, *args, **kwds)
