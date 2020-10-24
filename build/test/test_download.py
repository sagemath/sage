# -*- coding: utf-8 -*-
"""
Test downloading files
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
from __future__ import print_function, absolute_import

import unittest
import tempfile


from .capture import CapturedLog
from sage_bootstrap.download import Download, MirrorList
from sage_bootstrap.compat import StringIO


class DownloadTestCase(unittest.TestCase):

    def test_download_mirror_list(self):
        tmp = tempfile.NamedTemporaryFile()
        tmp.close()
        progress = StringIO()
        Download(MirrorList.URL, tmp.name, progress=progress).run()
        self.assertEqual(progress.getvalue(),
            '[......................................................................]\n')
        with open(tmp.name, 'r') as f:
            content = f.read()
            self.assertTrue(content.startswith('# Sage Mirror List'))

    def test_error(self):
        URL = 'http://files.sagemath.org/sage_bootstrap/this_url_does_not_exist'
        progress = StringIO()
        download = Download(URL, progress=progress)
        log = CapturedLog()

        def action():
            with log:
                download.run()
        self.assertRaises(IOError, action)
        self.assertIsNotFoundError(log.messages())
        self.assertEqual(progress.getvalue(),
            '[xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx]\n')

    def test_ignore_errors(self):
        URL = 'http://files.sagemath.org/sage_bootstrap/this_url_does_not_exist'
        with CapturedLog() as log:
            Download(URL, progress=False, ignore_errors=True).run()
        self.assertIsNotFoundError(log.messages())

    def assertIsNotFoundError(self, messages):
        self.assertEqual(len(messages), 1)
        self.assertEqual(messages[0][0], 'ERROR')
        self.assertTrue(messages[0][1].startswith('[Errno'))
        self.assertTrue(messages[0][1].endswith(
            "[Errno 404] Not Found: '//files.sagemath.org/sage_bootstrap/this_url_does_not_exist'"))
