# -*- coding: utf-8 -*-
"""
Test sage-download-file commandline utility
"""

# ****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import unittest
import subprocess
import logging

from sage_bootstrap.download.mirror_list import MIRRORLIST_FILENAME
from sage_bootstrap.util import is_url


log = logging.getLogger()


EXECUTABLE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'bin',
    'sage-download-file',
)


class SageDownloadFileTestCase(unittest.TestCase):

    maxDiff = None

    def test_print_mirror_list_no_network(self):
        """
        Subsequent runs of sage-download-file
        """
        try:
            os.remove(MIRRORLIST_FILENAME)
        except OSError:
            pass
        env = dict(os.environ)
        env['http_proxy'] = 'http://192.0.2.0:5187/'
        env['https_proxy'] = 'http://192.0.2.0:5187/'
        env['ftp_proxy'] = 'http://192.0.2.0:5187/'
        env['rsync_proxy'] = 'http://192.0.2.0:5187/'
        proc = subprocess.Popen(
            [EXECUTABLE, '--print-fastest-mirror', '--timeout=0.001'],
            stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=env,
        )
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        rc = proc.returncode
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints single url (the default) to stdout
        self.assertTrue(is_url(stdout))
        # Prints error to stderr
        self.assertTrue(stderr.find('Downloading the mirror list failed') >= 0)

    def test_print_mirror_list_timing(self):
        """
        The first run of sage-download-file
        """
        try:
            os.remove(MIRRORLIST_FILENAME)
        except OSError:
            pass
        proc = subprocess.Popen(
            [EXECUTABLE, '--print-fastest-mirror'],
            stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        rc = proc.returncode
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints single url to stdout
        self.assertTrue(is_url(stdout))
        # Prints mirrors to stderr
        self.assertTrue(len(stderr.strip().splitlines()) > 3)

    def test_print_mirror_list_cached(self):
        """
        Subsequent runs of sage-download-file
        """
        proc = subprocess.Popen(
            [EXECUTABLE, '--print-fastest-mirror'],
            stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        rc = proc.returncode
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints single url to stdout
        self.assertTrue(is_url(stdout))
        # May or may not print to stderr depending on whether cache was saved
