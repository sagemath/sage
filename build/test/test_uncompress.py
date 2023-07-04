# -*- coding: utf-8 -*-
"""
Test sage extraction of tarball / zip files
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import print_function, absolute_import

import os
import unittest
import shutil
import tempfile
import subprocess

from sage_bootstrap.uncompress.zip_file import SageZipFile
from sage_bootstrap.uncompress.tar_file import SageTarFile
from sage_bootstrap.uncompress.action import (
    open_archive, unpack_archive
)


class UncompressTarFileTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.filename = os.path.join(self.tmp, 'test.tar.gz')
        self.make_tarfile()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def make_tarfile(self):
        src = os.path.join(self.tmp, 'src')
        os.mkdir(src)
        os.mkdir(os.path.join(src, 'foo'))
        with open(os.path.join(src, 'content'), 'w+') as f:
            f.write('root file test')
        with open(os.path.join(src, 'foo', 'subcontent'), 'w+') as f:
            f.write('subdirectory file test')
        subprocess.check_call([
            'tar', 'czf', self.filename, 'content', 'foo'
        ], cwd=src)

    def test_can_read(self):
        self.assertTrue(SageTarFile.can_read(self.filename))
        self.assertFalse(SageZipFile.can_read(self.filename))

    def test_tarball(self):
        archive = open_archive(self.filename)
        content = archive.extractbytes('content')
        self.assertEqual(content, b'root file test')
        dst = os.path.join(self.tmp, 'dst')
        unpack_archive(archive, dst)
        subprocess.check_call([
            'diff', '-r', 'src', 'dst'
        ], cwd=self.tmp)


class UncompressZipFileTestCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.filename = os.path.join(self.tmp, 'test.zip')
        self.make_zipfile()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def make_zipfile(self):
        src = os.path.join(self.tmp, 'src')
        os.mkdir(src)
        os.mkdir(os.path.join(src, 'foo'))
        with open(os.path.join(src, 'content'), 'w+') as f:
            f.write('root file test')
        with open(os.path.join(src, 'foo', 'subcontent'), 'w+') as f:
            f.write('subdirectory file test')
        subprocess.check_call([
            'zip', '-q', '-r', self.filename, 'content', 'foo'
        ], cwd=src)

    def test_can_read(self):
        self.assertTrue(SageZipFile.can_read(self.filename))
        self.assertFalse(SageTarFile.can_read(self.filename))

    def test_zipfile(self):
        archive = open_archive(self.filename)
        content = archive.extractbytes('content')
        self.assertEqual(content, b'root file test')
        dst = os.path.join(self.tmp, 'dst')
        unpack_archive(archive, dst)
        subprocess.check_call([
            'diff', '-r', 'src', 'dst'
        ], cwd=self.tmp)
