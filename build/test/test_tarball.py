# -*- coding: utf-8 -*-
"""
Test Sage Third-Party Tarball Handling
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
import logging

from sage_bootstrap.package import Package
from sage_bootstrap.tarball import Tarball

from .capture import CapturedLog, CapturedOutput


log = logging.getLogger()


class TarballTestCase(unittest.TestCase):

    def test_tarball(self):
        pkg = Package('configure')
        tarball = Tarball(pkg.tarball_filename)
        self.assertEqual(tarball, pkg.tarball)
        self.assertEqual(pkg, tarball.package)
        with CapturedOutput() as (stdout, stderr):
            with CapturedLog() as _:
                tarball.download()
        self.assertEqual(stdout.getvalue(), '')
        self.assertTrue(tarball.checksum_verifies())

    def test_checksum(self):
        pkg = Package('configure')
        tarball = pkg.tarball
        with CapturedOutput() as (stdout, stderr):
            with CapturedLog() as log:
                tarball.download()
        self.assertTrue(tarball.checksum_verifies())
        with open(tarball.upstream_fqn, 'w') as f:
            f.write('foobar')
        self.assertFalse(tarball.checksum_verifies())
        with CapturedOutput() as (stdout, stderr):
            with CapturedLog() as log:
                tarball.download()
        msg = log.messages()
        self.assertTrue(
            ('INFO', 'Attempting to download package {0} from mirrors'.format(pkg.tarball_filename)) in msg)
        self.assertEqual(stdout.getvalue(), '')
        self.assertEqual(stderr.getvalue(),
            '[......................................................................]\n')
        self.assertTrue(tarball.checksum_verifies())
