# -*- coding: utf-8 -*-
"""
Test Sage Package Handling
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import unittest
from sage_bootstrap.package import Package
from sage_bootstrap.tarball import Tarball


class PackageTestCase(unittest.TestCase):

    maxDiff = None

    def test_package(self):
        pkg = Package('pari')
        self.assertTrue(pkg.name, 'pari')
        self.assertTrue(pkg.path.endswith('build/pkgs/pari'))
        self.assertEqual(pkg.tarball_pattern, 'pari-VERSION.tar.gz')
        self.assertEqual(pkg.tarball_filename, pkg.tarball.filename)
        self.assertTrue(pkg.tarball.filename.startswith('pari-') and
                        pkg.tarball.filename.endswith('.tar.gz'))
        self.assertTrue(pkg.tarball.filename.startswith('pari-') and
                        pkg.tarball.filename.endswith('.tar.gz'))
        self.assertTrue(isinstance(pkg.tarball, Tarball))

    def test_all(self):
        pari = Package('pari')
        self.assertTrue(pari in Package.all())
